from typing import Iterable, Dict, Callable, Any, Optional, List, Tuple, TYPE_CHECKING

import dataclasses

import stim

from gen._util import complex_key, sorted_complex

if TYPE_CHECKING:
    from gen._patch import Patch
    from gen._interaction_planner import InteractionPlanner
    from gen._interaction_planner import SingleQubitGatesPlanner


SYMMETRIC_GATES = {
    'CZ',
    'XCX',
    'YCY',
    'ZCZ',
    'SWAP',
    'ISWAP',
    'ISWAP_DAG',
    'SQRT_XX',
    'SQRT_YY',
    'SQRT_ZZ',
    'SQRT_XX_DAG',
    'SQRT_YY_DAG',
    'SQRT_ZZ_DAG',
}


@dataclasses.dataclass(frozen=True)
class AtLayer:
    """A special class that indicates the layer to read a measurement key from."""
    key: Any
    layer: Any


class MeasurementTracker:
    """Tracks measurements and groups of measurements, for producing stim record targets."""
    def __init__(self):
        self.recorded: Dict[Any, Optional[List[int]]] = {}
        self.next_measurement_index = 0

    def copy(self) -> 'MeasurementTracker':
        result = MeasurementTracker()
        result.recorded = {k: list(v) for k, v in self.recorded.items()}
        result.next_measurement_index = self.next_measurement_index
        return result

    def _rec(self, key: Any, value: Optional[List[int]]) -> None:
        if key in self.recorded:
            raise ValueError(f'Measurement key collision: {key=}')
        self.recorded[key] = value

    def record_measurement(self, key: Any) -> None:
        self._rec(key, [self.next_measurement_index])
        self.next_measurement_index += 1

    def make_measurement_group(self, sub_keys: Iterable[Any], *, key: Any) -> None:
        self._rec(key, self.measurement_indices(sub_keys))

    def record_obstacle(self, key: Any) -> None:
        self._rec(key, None)

    def measurement_indices(self, keys: Iterable[Any]) -> List[int]:
        result = set()
        for key in keys:
            if key not in self.recorded:
                raise ValueError(f"No such measurement: {key=}")
            for v in self.recorded[key]:
                if v is None:
                    raise ValueError(f"Obstacle at {key=}")
                if v in result:
                    result.remove(v)
                else:
                    result.add(v)
        return sorted(result)

    def current_measurement_record_targets_for(self, keys: Iterable[Any]) -> List[stim.GateTarget]:
        t0 = self.next_measurement_index
        times = self.measurement_indices(keys)
        return [stim.target_rec(t - t0) for t in sorted(times)]


class Builder:
    """Helper class for building stim circuits.

    Handles qubit indexing (complex -> int conversion).
    Handles measurement tracking (naming results and referring to them by name).
    """

    def __init__(self,
                 *,
                 q2i: Dict[complex, int],
                 circuit: stim.Circuit,
                 tracker: MeasurementTracker):
        self.q2i = q2i
        self.circuit = circuit
        self.tracker = tracker

    def copy(self) -> 'Builder':
        """Returns a Builder with independent copies of this builder's circuit and tracking data."""
        return Builder(q2i=dict(self.q2i), circuit=self.circuit.copy(), tracker=self.tracker.copy())

    def fork(self) -> 'Builder':
        """Returns a Builder with the same underlying tracking but which appends into a different circuit.
        """
        return Builder(q2i=self.q2i, circuit=stim.Circuit(), tracker=self.tracker)

    @staticmethod
    def for_qubits(
            qubits: Iterable[complex],
            *,
            to_circuit_coord_data: Callable[[complex], complex] = lambda e: e) -> 'Builder':
        q2i = {q: i for i, q in enumerate(sorted_complex(set(qubits)))}
        circuit = stim.Circuit()
        for q, i in q2i.items():
            c = to_circuit_coord_data(q)
            circuit.append("QUBIT_COORDS", [i], [c.real, c.imag])
        return Builder(
            q2i=q2i,
            circuit=circuit,
            tracker=MeasurementTracker(),
        )

    def gate(self,
             name: str,
             qubits: Iterable[complex]) -> None:
        assert name not in ['CZ', 'ZCZ', 'XCX', 'YCY', 'ISWAP', 'ISWAP_DAG', 'SWAP', 'M', 'MX', 'MY']
        qubits = sorted_complex(qubits)
        if not qubits:
            return
        self.circuit.append(name, [self.q2i[q] for q in qubits])

    def gate2(self,
              name: str,
              pairs: Iterable[Tuple[complex, complex]]) -> None:
        pairs = sorted(pairs, key=lambda pair: (complex_key(pair[0]), complex_key(pair[1])))
        if name == 'XCZ':
            pairs = [pair[::-1] for pair in pairs]
            name = 'CX'
        if name == 'YCZ':
            pairs = [pair[::-1] for pair in pairs]
            name = 'CY'
        if name == 'SWAPCX':
            pairs = [pair[::-1] for pair in pairs]
            name = 'CXSWAP'
        if name in SYMMETRIC_GATES:
            pairs = [sorted_complex(pair) for pair in pairs]
        if not pairs:
            return
        self.circuit.append(name, [self.q2i[q] for pair in pairs for q in pair])

    def shift_coords(self, *, dp: complex = 0, dt: int):
        self.circuit.append("SHIFT_COORDS", [], [dp.real, dp.imag, dt])

    def measure_patch(self, patch: 'Patch', *, save_layer: Any, cmp_layer: Optional[Any] = None) -> None:
        for tile in patch.tiles:
            self.measure_pauli_product(q2b={
                tile.ordered_data_qubits[k]: tile.bases[k]
                for k in range(len(tile.ordered_data_qubits))
                if tile.ordered_data_qubits[k] is not None
            }, key=AtLayer(tile.measurement_qubit, save_layer))
        if cmp_layer is not None:
            for tile in patch.tiles:
                m = tile.measurement_qubit
                self.detector([AtLayer(m, save_layer), AtLayer(m, cmp_layer)], pos=m)

    def demolition_measure_with_feedback_passthrough(
            self,
            xs: Iterable[complex] = (),
            ys: Iterable[complex] = (),
            zs: Iterable[complex] = (),
            *,
            tracker_key: Callable[[complex], Any] = lambda e: e,
            save_layer: Any) -> None:
        """Performs demolition measurements that look like measurements w.r.t. detectors.

        This is done by adding feedback operations that flip the demolished qubits depending
        on the measurement result. This feedback can then later be removed using
        stim.Circuit.with_inlined_feedback. The benefit is that it can be easier to
        programmatically create the detectors using the passthrough measurements, and
        then they can be automatically converted.
        """
        self.measure(qubits=xs, basis='X', tracker_key=tracker_key, save_layer=save_layer)
        self.measure(qubits=ys, basis='Y', tracker_key=tracker_key, save_layer=save_layer)
        self.measure(qubits=zs, basis='Z', tracker_key=tracker_key, save_layer=save_layer)
        self.tick()
        self.gate('RX', xs)
        self.gate('RY', ys)
        self.gate('RZ', zs)
        for (qs, b) in [(xs, 'Z'), (ys, 'X'), (zs, 'X')]:
            for q in qs:
                self.classical_paulis(control_keys=[AtLayer(tracker_key(q), save_layer)], targets=[q], basis=b)

    def measure(self,
                qubits: Iterable[complex],
                *,
                basis: str = 'Z',
                tracker_key: Callable[[complex], Any] = lambda e: e,
                save_layer: Any) -> None:
        qubits = sorted_complex(qubits)
        if not qubits:
            return
        self.circuit.append(f"M{basis}", [self.q2i[q] for q in qubits])
        for q in qubits:
            self.tracker.record_measurement(AtLayer(tracker_key(q), save_layer))

    def measure_pauli_product(self,
                              *,
                              xs: Iterable[complex] = (),
                              ys: Iterable[complex] = (),
                              zs: Iterable[complex] = (),
                              b2qs: Dict[str, Iterable[complex]] = None,
                              q2b: Dict[complex, str] = None,
                              noise: Optional[float] = None,
                              key: Any):
        """Adds an MPP operation to measure the given qubits. Supports a variety of formats.

        Note that all formats are combined as if multiplying Pauli observables (ignoring phase and
        sign).

        Args:
            xs: A list of qubits to include in the product as X basis terms.
            ys: A list of qubits to include in the product as Y basis terms.
            zs: A list of qubits to include in the product as Z basis terms.
            b2qs: A dictionary mapping Pauli bases ('X', 'Y', 'Z') to iterables of qubits. A mapping from the desired basis to qubits to include in the product in that basis.
            q2b: A dictionary mapping qubits to their corresponding Pauli bases.  A mapping from qubit to basis. Each qubit:basis pair is includedin the product.
            noise: Make the measurement noisy.
            key: Measurement key to track the result under.
        """
        x = set(xs)
        y = set(ys)
        z = set(zs)

        # Adds qubits to the appropriate Pauli basis set based on the b2qs mapping.
        if b2qs is not None:
            for b, bqs in b2qs.items():
                if b == 'X':
                    x |= set(bqs)
                elif b == 'Y':
                    y |= set(bqs)
                elif b == 'Z':
                    z |= set(bqs)
                else:
                    raise NotImplementedError(f'{b=}')
                
        # Adds individual qubit-Pauli pairs to the respective sets based on the q2b mapping.
        if q2b is not None:
            for q, b in q2b.items():
                if b == 'X':
                    x.add(q)
                elif b == 'Y':
                    y.add(q)
                elif b == 'Z':
                    z.add(q)
                else:
                    raise NotImplementedError(f'{b=}')
        xz = x & z  # Qubits in both X and Z
        xy = x & y  # Qubits in both X and Y
        yz = y & z  # Qubits in both Y and Z

        # •	Symmetric Difference Handling:
        # •	If a qubit is assigned both X and Y, it effectively becomes a Z operator, because  X \cdot Y = iZ .
        # •	If assigned both Y and Z, it becomes X ( Y \cdot Z = iX ).
        # •	If assigned both X and Z, it becomes Y ( X \cdot Z = -iY ).
        # •	The function adjusts the sets accordingly to ensure each qubit has a single Pauli operator.
        x -= xz
        x -= xy
        z -= xz
        z -= yz
        y -= xy
        y -= yz
        x |= yz
        y |= xz
        z |= xy

        # Create Mapping of Qubits to Pauli Targets
        # Maps each qubit to its corresponding stim Pauli target, using self.q2i to get the integer index of the qubit.
        vals = {}
        for q in x:
            vals[q] = stim.target_x(self.q2i[q])
        for q in y:
            vals[q] = stim.target_y(self.q2i[q])
        for q in z:
            vals[q] = stim.target_z(self.q2i[q])


        
        
        
        # •	Recording the Measurement:
        # •	The measurement is recorded in self.tracker with the provided key.
        targets = []
        # stim.target_combiner(): Represents the multiplication operator between Pauli operators in the MPP instruction.
        comb = stim.target_combiner() 

        # Building Targets:
        # For each qubit, the function appends the Pauli target and the combiner to the target list.
        # The combiner is appended after each Pauli target except the last one.
        for q in sorted_complex(vals.keys()):
            targets.append(vals[q])
            targets.append(comb)

        # Appending MPP Instruction:
        # If there are targets, it appends the MPP instruction to the circuit.
        if targets:
            targets.pop() # Remove the last combiner
            self.circuit.append('MPP', targets, noise)
            #The measurement is recorded in self.tracker with the provided key.
            self.tracker.record_measurement(key) 
        else:
            #If there are no qubits to measure (i.e., the Pauli product is the identity operator), the function creates an empty measurement group.
            self.tracker.make_measurement_group([], key=key) 

    def detector(self,
                 keys: Iterable[Any],
                 *,
                 pos: Optional[complex],
                 t: float = 0,
                 extra_coords: Iterable[float] = (),
                 mark_as_post_selected: bool = False,
                 ignore_non_existent: bool = False) -> None:
        """
        The detector function adds a DETECTOR instruction to the circuit, which defines 
        a detection event (syndrome bit) based on a set of measurement results. 
        In the context of error correction, detectors are used to identify errors by 
        comparing expected and actual measurement outcomes.

        Args:
            keys: An iterable of measurement keys whose results are used to define the detector.
            pos: An optional complex number specifying the spatial position of the detector for visualization purposes.
            t: A float specifying the time coordinate of the detector (default is 0).
            extra_coords: Additional coordinates for visualization.
            mark_as_post_selected: A boolean indicating if the detector is associated with post-selection (default is False).
            ignore_non_existent: If True, non-existent measurement keys are ignored; otherwise, an error is raised if a key does not exist (default is False).
        """
        if pos is not None: # If a position (pos) is provided, it constructs a list of coordinates for visualization.
            coords = [pos.real, pos.imag, t] + list(extra_coords) # The time coordinate t and any extra_coords are included.
            if mark_as_post_selected: # If mark_as_post_selected is True, an additional coordinate 1 is appended.
                coords.append(1)
        else: # Raises a ValueError if pos is None but extra_coords are provided or if mark_as_post_selected is True.
            if list(extra_coords): 
                raise ValueError('pos is None but extra_coords is not empty')
            if mark_as_post_selected:
                raise ValueError('pos is None and mark_as_post_selected')
            coords = None

        if ignore_non_existent: # Filters out keys that are not present in the measurement tracker if ignore_non_existent is True.
            keys = [k for k in keys if k in self.tracker.recorded]

        # Converts the measurement keys into stim record targets (e.g., rec[-1], rec[-2]), which refer to prior measurement outcomes.
        targets = self.tracker.current_measurement_record_targets_for(keys)
        # Adds the DETECTOR instruction to the circuit with the specified targets and coordinates.
        self.circuit.append('DETECTOR', targets, coords)

    def obs_include(self,
                    keys: Iterable[Any],
                    *,
                    obs_index: int) -> None:
        """
        The obs_include function adds an observable include instruction to the circuit. It specifies which 
        measurements contribute to a logical observable (e.g., a logical qubits Pauli operator) in the context 
        of quantum error correction.

        args:
            keys: An iterable of measurement keys that identify the measurement results to include in the observable.
		    obs_index: An integer index representing the specific observable (e.g., the logical qubit index).
        """
        # The function calls self.tracker.current_measurement_record_targets_for(keys) to get the measurement record 
        # targets corresponding to the provided keys. These targets are references to prior measurement outcomes in 
        # the circuit.
        ms = self.tracker.current_measurement_record_targets_for(keys)

        # If there are measurement targets (ms is not empty), it appends an OBSERVABLE_INCLUDE instruction to the circuit.
        if ms:
            # The OBSERVABLE_INCLUDE instruction specifies that the parity (XOR) of the specified measurement outcomes 
            # contributes to the logical observable identified by obs_index.
            self.circuit.append(
                'OBSERVABLE_INCLUDE',
                ms,
                obs_index,
            )

    def tick(self) -> None:
        """
        The tick method is used to advance the circuit's time step, separating different layers of operations. 
        Each tick represents a new time step where subsequent operations are considered to happen after the 
        previous operations. This is crucial in quantum circuits to maintain the correct temporal sequence of 
        operations.
        """
        self.circuit.append('TICK')

    def cz(self, pairs: List[Tuple[complex, complex]]) -> None:
        """The cz function adds Controlled-Z (CZ) gates between specified pairs of qubits in the circuit.
        
            args:
                pairs: A list of tuples, where each tuple contains two qubits (represented as complex numbers) between which a CZ gate should be applied.
        """
        sorted_pairs = []
        # For each pair (a, b), if a has a higher sort order than b (based on complex_key), they are swapped to ensure a consistent order.
        for a, b in pairs:
            if complex_key(a) > complex_key(b):
                a, b = b, a
            sorted_pairs.append((a, b))
        # The pairs are sorted based on the complex_key of the qubits to maintain a consistent and deterministic order of gate application.
        sorted_pairs = sorted(sorted_pairs, key=lambda e: (complex_key(e[0]), complex_key(e[1])))
        for a, b in sorted_pairs:
            # For each sorted pair, the function appends a 'CZ' gate to the circuit, mapping the qubits to their integer indices using self.q2i
            self.circuit.append('CZ', [self.q2i[a], self.q2i[b]])

    def swap(self, pairs: List[Tuple[complex, complex]]) -> None:
        """The swap function adds SWAP gates between specified pairs of qubits, effectively exchanging their quantum states.
        
        args:
            pairs: A list of tuples containing pairs of qubits to be swapped.
        """
        sorted_pairs = []
        for a, b in pairs:
            if complex_key(a) > complex_key(b):
                a, b = b, a
            sorted_pairs.append((a, b))
        sorted_pairs = sorted(sorted_pairs, key=lambda e: (complex_key(e[0]), complex_key(e[1])))
        for a, b in sorted_pairs:
            self.circuit.append('SWAP', [self.q2i[a], self.q2i[b]])

    def classical_paulis(self,
                         *,
                         control_keys: Iterable[Any],
                         targets: Iterable[complex],
                         basis: str) -> None:
        
        """
        The classical_paulis function applies classically controlled Pauli gates (CX, CY, or CZ) to 
        target qubits based on prior measurement results.

        args:
            control_keys: An iterable of measurement keys that serve as classical control bits.
		    targets: An iterable of qubits (complex numbers) on which the gates will be applied.
		    basis: A string specifying the Pauli basis ('X', 'Y', or 'Z') for the gate.
          
        """
        # Constructs the gate name by prefixing 'C' to the specified basis, resulting in 'CX', 'CY', or 'CZ'.
        gate = f'C{basis}'
        # Converts target qubits to their integer indices and sorts them for consistency.
        indices = [self.q2i[q] for q in sorted_complex(targets)]
        # Gets the measurement record targets corresponding to the control keys.
        for rec in self.tracker.current_measurement_record_targets_for(control_keys):
            for i in indices:
                # For each control measurement result (rec) and each target qubit index (i),
                # appends the controlled Pauli gate to the circuit.
                self.circuit.append(gate, [rec, i])

    def plan_interactions(
            self,
            *,
            layer_count: int,
            start_orientations: Optional[Dict[complex, str]] = None,
            end_orientations: Optional[Dict[complex, str]] = None,
    ) -> 'InteractionPlanner':
        """
        The plan_interactions function initializes an InteractionPlanner to schedule and 
        plan multi-qubit interactions over multiple layers in the circuit.

        args:
            layer_count: The number of layers (time steps) to plan interactions over.
		    start_orientations: An optional dictionary mapping qubits to their initial orientations (states).
		    end_orientations: An optional dictionary mapping qubits to their desired orientations at the end of the interactions.
        """

        # Imports the InteractionPlanner class and instantiates it with the provided parameters and the current builder instance.
        from gen._interaction_planner import InteractionPlanner
        # The planner can then be used to define and schedule specific interactions.
        return InteractionPlanner(
            layer_count=layer_count,
            builder=self,
            start_orientations=start_orientations,
            end_orientations=end_orientations,
        )

    def plan_rotations(self) -> 'SingleQubitGatesPlanner':
        """The plan_rotations function initializes a SingleQubitGatesPlanner to schedule and 
        plan single-qubit rotations (gates) in the circuit."""
        from gen._interaction_planner import SingleQubitGatesPlanner
        # The planner allows you to specify and schedule single-qubit gates.
        return SingleQubitGatesPlanner(self)
