from os import name
from typing import Optional, Dict, Set, List, Iterator, Union, AbstractSet, DefaultDict, Any

import collections

from math import sin, cos, pi

import stim

CLIFFORD_1Q = 'C1'
CLIFFORD_2Q = 'C2'
ANNOTATION = 'info'
MPP = 'MPP'
MEASURE_RESET_1Q = 'MR1'
JUST_MEASURE_1Q = 'M1'
JUST_RESET_1Q = 'R1'
NOISE = '!?'

OP_TYPES = {
    'I': CLIFFORD_1Q,
    'X': CLIFFORD_1Q,
    'Y': CLIFFORD_1Q,
    'Z': CLIFFORD_1Q,
    'C_XYZ': CLIFFORD_1Q,
    'C_ZYX': CLIFFORD_1Q,
    'H': CLIFFORD_1Q,
    'H_XY': CLIFFORD_1Q,
    'H_XZ': CLIFFORD_1Q,
    'H_YZ': CLIFFORD_1Q,
    'S': CLIFFORD_1Q,
    'SQRT_X': CLIFFORD_1Q,
    'SQRT_X_DAG': CLIFFORD_1Q,
    'SQRT_Y': CLIFFORD_1Q,
    'SQRT_Y_DAG': CLIFFORD_1Q,
    'SQRT_Z': CLIFFORD_1Q,
    'SQRT_Z_DAG': CLIFFORD_1Q,
    'S_DAG': CLIFFORD_1Q,

    'CNOT': CLIFFORD_2Q,
    'CX': CLIFFORD_2Q,
    'CY': CLIFFORD_2Q,
    'CZ': CLIFFORD_2Q,
    'ISWAP': CLIFFORD_2Q,
    'ISWAP_DAG': CLIFFORD_2Q,
    'CXSWAP': CLIFFORD_2Q,
    'SWAPCX': CLIFFORD_2Q,
    'SQRT_XX': CLIFFORD_2Q,
    'SQRT_XX_DAG': CLIFFORD_2Q,
    'SQRT_YY': CLIFFORD_2Q,
    'SQRT_YY_DAG': CLIFFORD_2Q,
    'SQRT_ZZ': CLIFFORD_2Q,
    'SQRT_ZZ_DAG': CLIFFORD_2Q,
    'SWAP': CLIFFORD_2Q,
    'XCX': CLIFFORD_2Q,
    'XCY': CLIFFORD_2Q,
    'XCZ': CLIFFORD_2Q,
    'YCX': CLIFFORD_2Q,
    'YCY': CLIFFORD_2Q,
    'YCZ': CLIFFORD_2Q,
    'ZCX': CLIFFORD_2Q,
    'ZCY': CLIFFORD_2Q,
    'ZCZ': CLIFFORD_2Q,

    'MPP': MPP,
    'MR': MEASURE_RESET_1Q,
    'MRX': MEASURE_RESET_1Q,
    'MRY': MEASURE_RESET_1Q,
    'MRZ': MEASURE_RESET_1Q,
    'M': JUST_MEASURE_1Q,
    'MX': JUST_MEASURE_1Q,
    'MY': JUST_MEASURE_1Q,
    'MZ': JUST_MEASURE_1Q,
    'R': JUST_RESET_1Q,
    'RX': JUST_RESET_1Q,
    'RY': JUST_RESET_1Q,
    'RZ': JUST_RESET_1Q,

    'DETECTOR': ANNOTATION,
    'OBSERVABLE_INCLUDE': ANNOTATION,
    'QUBIT_COORDS': ANNOTATION,
    'SHIFT_COORDS': ANNOTATION,
    'TICK': ANNOTATION,
    'E': ANNOTATION,

    'DEPOLARIZE1': NOISE,
    'DEPOLARIZE2': NOISE,
    'PAULI_CHANNEL_1': NOISE,
    'PAULI_CHANNEL_2': NOISE,
    'X_ERROR': NOISE,
    'Y_ERROR': NOISE,
    'Z_ERROR': NOISE,
    # Not supported.
    # 'CORRELATED_ERROR': NOISE,
    # 'E': NOISE,
    # 'ELSE_CORRELATED_ERROR',
}
OP_MEASURE_BASES = {
    'M': 'Z',
    'MX': 'X',
    'MY': 'Y',
    'MZ': 'Z',
    'MPP': '',
}
COLLAPSING_OPS = {op for op, t in OP_TYPES.items() if t == JUST_RESET_1Q or t == JUST_MEASURE_1Q or t == MPP or t == MEASURE_RESET_1Q}


class NoiseRule:
    """Describes how to add noise to an operation."""

    def __init__(self,
                 *,
                 after: Dict[str, float],
                 flip_result: float = 0):
        """
        Args:
            after: A dictionary mapping noise rule names to their probability argument.
                For example, {"DEPOLARIZE2": 0.01, "X_ERROR": 0.02} will add two qubit
                depolarization with parameter 0.01 and also add 2% bit flip noise. These
                noise channels occur after all other operations in the moment and are applied
                to the same targets as the relevant operation.
            flip_result: The probability that a measurement result should be reported incorrectly.
                Only valid when applied to operations that produce measurement results.
        """
        if not (0 <= flip_result <= 1):
            raise ValueError(f'not (0 <= {flip_result=} <= 1)')
        for k, p in after.items():
            if OP_TYPES[k] != NOISE:
                raise ValueError(f'not a noise channel: {k} from {after=}')
            if not (0 <= p <= 1):
                raise ValueError(f'not (0 <= {p} <= 1) from {after=}')
        self.after = after
        self.flip_result = flip_result

    def append_noisy_version_of(self,
                                *,
                                split_op: stim.CircuitInstruction,
                                out_during_moment: stim.Circuit,
                                after_moments: DefaultDict[Any, stim.Circuit],
                                immune_qubits: AbstractSet[int]) -> None:
        targets = split_op.targets_copy()
        if immune_qubits and any((t.is_qubit_target or t.is_x_target or t.is_y_target or t.is_z_target) and t.value in immune_qubits for t in targets):
            out_during_moment.append(split_op)
            return

        args = split_op.gate_args_copy()
        if self.flip_result:
            t = OP_TYPES[split_op.name]
            assert t == MPP or t == JUST_MEASURE_1Q or t == MEASURE_RESET_1Q
            assert len(args) == 0
            args = [self.flip_result]

        out_during_moment.append(split_op.name, targets, args)
        raw_targets = [t.value for t in targets if not t.is_combiner]
        for op_name, arg in self.after.items():
            after_moments[(op_name, arg)].append(op_name, raw_targets, arg)


class NoiseModel:
    def __init__(self,
                 idle_depolarization: float = 0,
                 tick_noise: Optional[NoiseRule] = None,
                 additional_depolarization_waiting_for_m_or_r: float = 0,
                 gate_rules: Optional[Dict[str, NoiseRule]] = None,
                 measure_rules: Optional[Dict[str, NoiseRule]] = None,
                 any_measurement_rule: Optional[NoiseRule] = None,
                 any_clifford_1q_rule: Optional[NoiseRule] = None,
                 any_clifford_2q_rule: Optional[NoiseRule] = None,
                 allow_multiple_uses_of_a_qubit_in_one_tick: bool = False):
        self.idle_depolarization = idle_depolarization
        self.tick_noise = tick_noise
        self.additional_depolarization_waiting_for_m_or_r = additional_depolarization_waiting_for_m_or_r
        self.gate_rules = {} if gate_rules is None else gate_rules
        self.measure_rules = measure_rules
        self.any_measurement_rule = any_measurement_rule
        self.any_clifford_1q_rule = any_clifford_1q_rule
        self.any_clifford_2q_rule = any_clifford_2q_rule
        self.allow_multiple_uses_of_a_qubit_in_one_tick = allow_multiple_uses_of_a_qubit_in_one_tick
        assert self.tick_noise is None or not self.tick_noise.flip_result

    @staticmethod
    def si1000(p: float) -> 'NoiseModel':
        """Superconducting inspired noise.

        As defined in "A Fault-Tolerant Honeycomb Memory" https://arxiv.org/abs/2108.10457

        Small tweak when measurements aren't immediately followed by a reset: the measurement result
        is probabilistically flipped instead of the input qubit. The input qubit is depolarized after
        the measurement.
        """
        return NoiseModel(
            idle_depolarization=p / 10,
            additional_depolarization_waiting_for_m_or_r=2 * p,
            any_clifford_1q_rule=NoiseRule(after={'DEPOLARIZE1': p / 10}),
            any_clifford_2q_rule=NoiseRule(after={'DEPOLARIZE2': p}),
            measure_rules={
                'Z': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p * 5),
                'ZZ': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=p * 5),
            },
            gate_rules={
                'R': NoiseRule(after={'X_ERROR': p * 2}),
            }
        )

    @staticmethod
    def heavy_hex() -> 'NoiseModel':
        """
        Noise model tuned to the heavy-hex device parameters (Table I).

        ┌────────────────────────┬────────┐
        │ parameter              │ value  │
        ├────────────────────────┼────────┤
        │ p₁Q   (1-qubit gate)   │ 0.02 % │
        │ p₂Q   (2-qubit gate)   │ 0.41 % │
        │ p_q.meas (quantum)     │ 1.2 %  │
        │ p_c.meas (classical)   │ 4.2 %  │
        │ p_idle                 │ 1.2 %  │
        │ p_reset                │ 7.5 %  │
        └────────────────────────┴────────┘

        *   Gate errors are modelled as depolarisation **before** each Clifford gate.
        *   Measurement errors are split:
            - quantum part → an `X_ERROR` on the measured qubit with p_q.meas
            - classical part → result flip probability p_c.meas
        *   Data qubits that are idle while other qubits are being measured/reset
            suffer depolarising noise with p_idle.
        *   Reset is followed by an `X_ERROR` with p_reset.
        """
        p_1q = 0.0002      # 0.02 %
        p_2q = 0.0041      # 0.41 %
        p_q_meas = 0.012   # 1.2 %
        p_c_meas = 0.042   # 4.2 %
        p_idle = 0.012     # 1.2 %
        
        p_reset = 0.075    # 7.5 %

        return NoiseModel(
            # Idling depolarisation (applies during measurement / reset periods)
            idle_depolarization=p_idle,

            # Clifford‑gate noise
            any_clifford_1q_rule=NoiseRule(after={'DEPOLARIZE1': p_1q}),
            any_clifford_2q_rule=NoiseRule(after={'DEPOLARIZE2': p_2q}),

            # Measurement noise (quantum part + classical result flip)
            measure_rules={
                'Z':  NoiseRule(after={'X_ERROR': p_q_meas}, flip_result=p_c_meas),
                'X':  NoiseRule(after={'X_ERROR': p_q_meas}, flip_result=p_c_meas),
                'Y':  NoiseRule(after={'X_ERROR': p_q_meas}, flip_result=p_c_meas),
                'ZZ': NoiseRule(after={'DEPOLARIZE2': p_q_meas}, flip_result=p_c_meas),
                'XX': NoiseRule(after={'DEPOLARIZE2': p_q_meas}, flip_result=p_c_meas),
                'YY': NoiseRule(after={'DEPOLARIZE2': p_q_meas}, flip_result=p_c_meas),
            },

            # Reset noise
            gate_rules={
                'R':  NoiseRule(after={'X_ERROR': p_reset}),
                'RX': NoiseRule(after={'Z_ERROR': p_reset}),
                'RY': NoiseRule(after={'X_ERROR': p_reset}),
            }
        )
    @staticmethod
    def IBM1907(p: float) -> 'NoiseModel':
        """Circuit-level depolarizing noise model based on arXiv:1907.09528 and flag fault tolerance assumptions.

        1. Each 1-qubit gate is followed by a random Pauli {X,Y,Z} with probability p.
        2. Each 2-qubit gate is followed by a uniform Pauli error on 2 qubits (not I⊗I) with probability p.
        3. |0⟩ and |+⟩ state prep is faulty with probability 2p/3, flipped to |1⟩ or |−⟩ respectively.
        4. Measurement results are flipped with probability 2p/3.
        5. Idling qubits also suffer random Pauli {X,Y,Z} errors with probability p.
        """
        return NoiseModel(
            idle_depolarization=p,
            any_clifford_1q_rule=NoiseRule(after={'DEPOLARIZE1': p}),
            any_clifford_2q_rule=NoiseRule(after={'DEPOLARIZE2': p}),
            measure_rules={
                'X': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=2 * p / 3),
                'Y': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=2 * p / 3),
                'Z': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=2 * p / 3),
                'XX': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=2 * p / 3),
                'YY': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=2 * p / 3),
                'ZZ': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=2 * p / 3),
            },
            gate_rules={
                'RX': NoiseRule(after={'Z_ERROR': 2 * p / 3}),
                'RY': NoiseRule(after={'X_ERROR': 2 * p / 3}),
                'R':  NoiseRule(after={'X_ERROR': 2 * p / 3}),
            }
        )
    @staticmethod
    def uniform_depolarizing(p: float) -> 'NoiseModel':
        """Near-standard circuit depolarizing noise.

        Everything has the same parameter p.
        Single qubit clifford gates get single qubit depolarization.
        Two qubit clifford gates get single qubit depolarization.
        Dissipative gates have their result probabilistically bit flipped (or phase flipped if appropriate).

        Non-demolition measurement is treated a bit unusually in that it is the result that is flipped instead of
        the input qubit. The input qubit is depolarized.
        """
        return NoiseModel(
            idle_depolarization=p,
            any_clifford_1q_rule=NoiseRule(after={'DEPOLARIZE1': p}),
            any_clifford_2q_rule=NoiseRule(after={'DEPOLARIZE2': p}),
            measure_rules={
                'X': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p),
                'Y': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p),
                'Z': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p),
                'XX': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=p),
                'YY': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=p),
                'ZZ': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=p),
            },
            gate_rules={
                'RX': NoiseRule(after={'Z_ERROR': p}),
                'RY': NoiseRule(after={'X_ERROR': p}),
                'R': NoiseRule(after={'X_ERROR': p}),
            }
        )


    @staticmethod
    def corr_mixing(
        # Global background and sweep
        base_prob: float,                  # default probability for all unspecified slots
        sweep_slot: str,                   # the slot whose probability you sweep to find a threshold
        sweep_prob: float,                 # current value for the sweep slot (varies when searching threshold)

        # Two independently biased (fixed) slots for the X/Y axes of your 3D surface
        bias_slot_1: str,                  # e.g. 'D1' , 'M' or 'D2'
        bias_slot_1_prob: float,           # direct probability for bias_slot_1 (this is your X axis value)
        bias_slot_2: str,                  # e.g. 'D1' , 'M' or 'D2'
        bias_slot_2_prob: float,           # direct probability for bias_slot_2 (this is your Y axis value)
    ) -> 'NoiseModel':
        """
        Construct a noise model with:
        1) a global baseline probability (base_prob) applied to all slots by default,
        2) two *fixed* biased slots set to explicit probabilities (bias_slot_1_prob, bias_slot_2_prob),
        3) one *sweep* slot whose probability (sweep_prob) is varied externally when searching for a threshold.

        D1 = 'DEPOLARIZE1'
        D2 = 'DEPOLARIZE2'
        M  = 'M' (measurement result flip)
        idle = 'idle' (idling depolarization)
        Parameters
        ----------
        base_prob : float
            Probability used for every noise slot not explicitly assigned below.
            Clamped into [0, 1].

        sweep_slot : str
            Name of the noise slot that is tied to `sweep_prob`. This slot is varied
            during threshold searches. Must be distinct from `bias_slot_1` and `bias_slot_2`.

        sweep_prob : float
            Current probability for `sweep_slot`. In a threshold search you will call
            this function repeatedly with different `sweep_prob` values. Clamped into [0, 1].

        bias_slot_1 : str
            Name of the first biased noise slot. This is what you place on the X axis
            when producing a 3D surface.

        bias_slot_1_prob : float
            Explicit probability for `bias_slot_1`. This value is fixed while you
            sweep `sweep_prob`. Clamped into [0, 1].

        bias_slot_2 : str
            Name of the second biased noise slot. This is what you place on the Y axis
            when producing a 3D surface.

        bias_slot_2_prob : float
            Explicit probability for `bias_slot_2`. This value is fixed while you
            sweep `sweep_prob`. Clamped into [0, 1].

        Returns
        -------
        NoiseModel
            A model where each slot probability is determined by the selector below:
            - if name == sweep_slot:        sweep_prob
            - elif name == bias_slot_1:     bias_slot_1_prob
            - elif name == bias_slot_2:     bias_slot_2_prob
            - else:                         base_prob

        Constraints
        -----------
        - All three slots must be distinct: {sweep_slot, bias_slot_1, bias_slot_2} must have size 3.
        - Values are clamped into [0, 1].
        - If two names collide, a ValueError is raised.

        Measure/Reset Semantics
        -----------------------
        - Measurement flip probabilities are taken from sel('M').
        - Reset uses the standard union of independent failure modes:
            p_reset = p_measure + p_idle - p_measure * p_idle
        where p_measure = sel('M') and p_idle = sel('idle').
        - 1q Clifford noise uses sel('DEPOLARIZE1'), 2q Clifford uses sel('DEPOLARIZE2').
        Adjust these keys if your backend uses different slot names.

        Raises
        ------
        ValueError
            If any of the three slot names are not distinct.

        Examples
        --------
        >>> # Build a model at a single grid point (px, py) while the threshold search
        >>> # will vary sweep_prob elsewhere:
        >>> model = corr_mixing_direct(
        ...     base_prob=1e-3,
        ...     sweep_slot='D2',
        ...     sweep_prob=sweep_prob,                 # vary this
        ...     bias_slot_1='D1',      # X axis
        ...     bias_slot_1_prob=0.0025,        # px
        ...     bias_slot_2='M',                # Y axis
        ...     bias_slot_2_prob=0.0010         # py
        ... )
        """
        def clamp(p: float) -> float:
            return max(0.0, min(1.0, float(p)))

        base_prob = clamp(base_prob)
        sweep_prob = clamp(sweep_prob)
        bias_slot_1_prob = clamp(bias_slot_1_prob)
        bias_slot_2_prob = clamp(bias_slot_2_prob)

        if len({sweep_slot, bias_slot_1, bias_slot_2}) != 3:
            raise ValueError(
                f"Slot names must be distinct: sweep_slot='{sweep_slot}', "
                f"bias_slot_1='{bias_slot_1}', bias_slot_2='{bias_slot_2}'."
            )

        def sel(name: str) -> float:
            if name == sweep_slot:
                return sweep_prob
            if name == bias_slot_1:
                return bias_slot_1_prob
            if name == bias_slot_2:
                return bias_slot_2_prob
            return base_prob

        def reset_prob() -> float:
            p_m = sel('M')
            p_i = sel('idle')
            return p_m + p_i - p_m * p_i

        return NoiseModel(
            idle_depolarization=sel('idle'),
            any_clifford_1q_rule=NoiseRule(after={'DEPOLARIZE1': sel('D1')}),
            any_clifford_2q_rule=NoiseRule(after={'DEPOLARIZE2': sel('D2')}),
            measure_rules={
                'X':  NoiseRule(after={'DEPOLARIZE1': sel('D1')}, flip_result=sel('M')),
                'Y':  NoiseRule(after={'DEPOLARIZE1': sel('D1')}, flip_result=sel('M')),
                'Z':  NoiseRule(after={'DEPOLARIZE1': sel('D1')}, flip_result=sel('M')),
                # 'XX': NoiseRule(after={'DEPOLARIZE2': sel('DEPOLARIZE2')}, flip_result=sel('XX')),
                # 'YY': NoiseRule(after={'DEPOLARIZE2': sel('DEPOLARIZE2')}, flip_result=sel('YY')),
                # 'ZZ': NoiseRule(after={'DEPOLARIZE2': sel('DEPOLARIZE2')}, flip_result=sel('ZZ')),
            },
            gate_rules={
                # Reset and Rx/Ry map to specific Pauli error channels in this codebase.
                'RX': NoiseRule(after={'Z_ERROR': reset_prob()}),
                # 'RY': NoiseRule(after={'X_ERROR': sel('X_ERROR')}),
                'R':  NoiseRule(after={'X_ERROR': reset_prob()}),
            }
        )

    @staticmethod
    def depolarizing_two_body_measurement_noise(p: float) -> 'NoiseModel':
        return NoiseModel(
            idle_depolarization=p,
            any_clifford_1q_rule=NoiseRule(after={'DEPOLARIZE1': p}),
            measure_rules={
                'XX': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=p),
                'YY': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=p),
                'ZZ': NoiseRule(after={'DEPOLARIZE2': p}, flip_result=p),
                'X': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p),
                'Y': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p),
                'Z': NoiseRule(after={'DEPOLARIZE1': p}, flip_result=p),
            },
            gate_rules={
                'RX': NoiseRule(after={'Z_ERROR': p}),
                'RY': NoiseRule(after={'X_ERROR': p}),
                'R': NoiseRule(after={'X_ERROR': p}),
            }
        )

    def _noise_rule_for_split_operation(self, *, split_op: stim.CircuitInstruction) -> Optional[NoiseRule]:
        if occurs_in_classical_control_system(split_op):
            return None

        rule = self.gate_rules.get(split_op.name)
        if rule is not None:
            return rule

        t = OP_TYPES[split_op.name]

        if self.any_clifford_1q_rule is not None and t == CLIFFORD_1Q:
            return self.any_clifford_1q_rule
        if self.any_clifford_2q_rule is not None and t == CLIFFORD_2Q:
            return self.any_clifford_2q_rule
        if self.measure_rules is not None:
            rule = self.measure_rules.get(_measure_basis(split_op=split_op))
            if rule is not None:
                return rule
        if self.any_measurement_rule is not None and t in [JUST_MEASURE_1Q, MEASURE_RESET_1Q, MPP]:
            return self.any_measurement_rule

        raise ValueError(f"No noise (or lack of noise) specified for {split_op=}.")

    def _append_idle_error(self,
                           *,
                           moment_split_ops: List[stim.CircuitInstruction],
                           out: stim.Circuit,
                           system_qubits: AbstractSet[int],
                           immune_qubits: AbstractSet[int],
                           ) -> None:
        collapse_qubits = []
        clifford_qubits = []
        for split_op in moment_split_ops:
            if occurs_in_classical_control_system(split_op):
                continue
            if split_op.name in COLLAPSING_OPS:
                qubits_out = collapse_qubits
            else:
                qubits_out = clifford_qubits
            for target in split_op.targets_copy():
                if not target.is_combiner:
                    qubits_out.append(target.value)

        # Safety check for operation collisions.
        usage_counts = collections.Counter(collapse_qubits + clifford_qubits)
        qubits_used_multiple_times = {q for q, c in usage_counts.items() if c != 1}
        if qubits_used_multiple_times and not self.allow_multiple_uses_of_a_qubit_in_one_tick:
            moment = stim.Circuit()
            for op in moment_split_ops:
                moment.append(op)
            raise ValueError(f"Qubits were operated on multiple times without a TICK in between:\n"
                             f"multiple uses: {sorted(qubits_used_multiple_times)}\n"
                             f"moment:\n"
                             f"{moment}")

        collapse_qubits_set = set(collapse_qubits)
        clifford_qubits_set = set(clifford_qubits)
        idle = sorted(system_qubits - collapse_qubits_set - clifford_qubits_set - immune_qubits)
        if idle and self.idle_depolarization:
            out.append('DEPOLARIZE1', idle, self.idle_depolarization)

        waiting_for_mr = sorted(system_qubits - collapse_qubits_set - immune_qubits)
        if collapse_qubits_set and waiting_for_mr and self.additional_depolarization_waiting_for_m_or_r:
            out.append('DEPOLARIZE1', idle, self.additional_depolarization_waiting_for_m_or_r)

        if self.tick_noise is not None:
            for k, p in self.tick_noise.after.items():
                out.append(k, system_qubits, p)

    def _append_noisy_moment(self,
                             *,
                             moment_split_ops: List[stim.CircuitInstruction],
                             out: stim.Circuit,
                             system_qubits: AbstractSet[int],
                             immune_qubits: AbstractSet[int],
                             ) -> None:
        after = collections.defaultdict(stim.Circuit)
        for split_op in moment_split_ops:
            rule = self._noise_rule_for_split_operation(split_op=split_op)
            if rule is None:
                out.append(split_op)
            else:
                rule.append_noisy_version_of(
                    split_op=split_op,
                    out_during_moment=out,
                    after_moments=after,
                    immune_qubits=immune_qubits,
                )
        for k in sorted(after.keys()):
            out += after[k]

        self._append_idle_error(
            moment_split_ops=moment_split_ops,
            out=out,
            system_qubits=system_qubits,
            immune_qubits=immune_qubits,
        )

    def noisy_circuit(self,
                      circuit: stim.Circuit,
                      *,
                      system_qubits: Optional[Set[int]] = None,
                      immune_qubits: Optional[AbstractSet[int]] = None,
                      ) -> stim.Circuit:
        """Returns a noisy version of the given circuit, by applying the receiving noise model.

        Args:
            circuit: The circuit to layer noise over.
            system_qubits: All qubits used by the circuit. These are the qubits eligible for idling noise.
            immune_qubits: Qubits to not apply noise to, even if they are operated on.

        Returns:
            The noisy version of the circuit.
        """
        if system_qubits is None:
            system_qubits = set(range(circuit.num_qubits))
        if immune_qubits is None:
            immune_qubits = set()

        result = stim.Circuit()

        first = True
        for moment_split_ops in _iter_split_op_moments(circuit, immune_qubits=immune_qubits):
            if first:
                first = False
            elif result and isinstance(result[-1], stim.CircuitRepeatBlock):
                pass
            else:
                result.append('TICK')
            if isinstance(moment_split_ops, stim.CircuitRepeatBlock):
                noisy_body = self.noisy_circuit(
                    moment_split_ops.body_copy(),
                    system_qubits=system_qubits,
                    immune_qubits=immune_qubits,
                )
                noisy_body.append('TICK')
                result.append(stim.CircuitRepeatBlock(repeat_count=moment_split_ops.repeat_count, body=noisy_body))
            else:
                self._append_noisy_moment(
                    moment_split_ops=moment_split_ops,
                    out=result,
                    system_qubits=system_qubits,
                    immune_qubits=immune_qubits,
                )

        return result


def occurs_in_classical_control_system(op: stim.CircuitInstruction) -> bool:
    """Determines if an operation is an annotation or a classical control system update."""
    t = OP_TYPES[op.name]
    if t == ANNOTATION:
        return True
    if t == CLIFFORD_2Q:
        targets = op.targets_copy()
        for k in range(0, len(targets), 2):
            a = targets[k]
            b = targets[k + 1]
            classical_0 = a.is_measurement_record_target or a.is_sweep_bit_target
            classical_1 = b.is_measurement_record_target or b.is_sweep_bit_target
            if not (classical_0 or classical_1):
                return False
        return True
    return False


def _split_targets_if_needed(op: stim.CircuitInstruction, immune_qubits: AbstractSet[int]) -> List[stim.CircuitInstruction]:
    """Splits operations into pieces as needed (e.g. MPP into each product, classical control away from quantum ops)."""
    t = OP_TYPES[op.name]
    if t == CLIFFORD_2Q:
        yield from _split_targets_if_needed_clifford_2q(op, immune_qubits)
    elif t == MPP:
        yield from _split_targets_if_needed_m_basis(op, immune_qubits)
    elif t in [NOISE, ANNOTATION]:
        yield op
    else:
        yield from _split_targets_if_needed_clifford_1q(op, immune_qubits)


def _split_targets_if_needed_clifford_1q(op: stim.CircuitInstruction, immune_qubits: AbstractSet[int]) -> List[stim.CircuitInstruction]:
    if immune_qubits:
        args = op.gate_args_copy()
        for t in op.targets_copy():
            yield stim.CircuitInstruction(op.name, [t], args)
    else:
        yield op


def _split_targets_if_needed_clifford_2q(op: stim.CircuitInstruction, immune_qubits: AbstractSet[int]) -> List[stim.CircuitInstruction]:
    """Splits classical control system operations away from things actually happening on the quantum computer."""
    assert OP_TYPES[op.name] == CLIFFORD_2Q
    targets = op.targets_copy()
    if immune_qubits or any(t.is_measurement_record_target for t in targets):
        args = op.gate_args_copy()
        for k in range(0, len(targets), 2):
            yield stim.CircuitInstruction(op.name, targets[k:k+2], args)
    else:
        yield op


def _split_targets_if_needed_m_basis(op: stim.CircuitInstruction, immune_qubits: AbstractSet[int]) -> List[stim.CircuitInstruction]:
    """Splits an MPP operation into one operation for each Pauli product it measures."""
    targets = op.targets_copy()
    args = op.gate_args_copy()
    k = 0
    start = k
    while k < len(targets):
        if k + 1 == len(targets) or not targets[k + 1].is_combiner:
            yield stim.CircuitInstruction(op.name, targets[start:k + 1], args)
            k += 1
            start = k
        else:
            k += 2
    assert k == len(targets)


def _iter_split_op_moments(circuit: stim.Circuit, *, immune_qubits: AbstractSet[int]) -> Iterator[Union[stim.CircuitRepeatBlock, List[stim.CircuitInstruction]]]:
    """Splits a circuit into moments and some operations into pieces.

    Classical control system operations like CX rec[-1] 0 are split from quantum operations like CX 1 0.

    MPP operations are split into one operation per Pauli product.

    Yields:
        Lists of operations corresponding to one moment in the circuit, with any problematic operations
        like MPPs split into pieces.

        (A moment is the time between two TICKs.)
    """
    cur_moment = []

    for op in circuit:
        if isinstance(op, stim.CircuitRepeatBlock):
            if cur_moment:
                yield cur_moment
                cur_moment = []
            yield op
        elif isinstance(op, stim.CircuitInstruction):
            if op.name == 'TICK':
                yield cur_moment
                cur_moment = []
            else:
                cur_moment.extend(_split_targets_if_needed(op, immune_qubits=immune_qubits))
    if cur_moment:
        yield cur_moment


def _measure_basis(*, split_op: stim.CircuitInstruction) -> Optional[str]:
    """Converts an operation into a string describing the Pauli product basis it measures.

    Returns:
        None: This is not a measurement (or not *just* a measurement).
        str: Pauli product string that the operation measures (e.g. "XX" or "Y").
    """
    result = OP_MEASURE_BASES.get(split_op.name)
    targets = split_op.targets_copy()
    if result == '':
        for k in range(0, len(targets), 2):
            t = targets[k]
            if t.is_x_target:
                result += 'X'
            elif t.is_y_target:
                result += 'Y'
            elif t.is_z_target:
                result += 'Z'
            else:
                raise NotImplementedError(f'{targets=}')
    return result
