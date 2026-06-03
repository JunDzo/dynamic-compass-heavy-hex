"""
Heavy-hex circuit generator (less-gates variant).

This module is a drop-in replacement for `_heavy_hex.py` with the same public
API (`make_heavy_hex_circuit`, `Tile`, etc.), but uses a modified X-check
schedule that removes two CX gates per 4-body tile.

Key differences vs `_heavy_hex.py`:
- The 4-body X-check routine `_add_x_circuit` uses a shorter CX schedule:
  some “cleanup” CX gates are omitted.
- Correctness is intended under the assumption of mid-circuit resets on
  measurement and flag qubits; without resets the same pattern can produce
  non-deterministic detectors and `stim.detector_error_model(...)` issues.
- Geometry (tiles, coordinates, qubit layout), 2-body Z-check construction,
  and detector construction logic are kept identical to `_heavy_hex.py` to
  isolate the effect of removing CX gates.

Use this file when you want to study the trade-off between fewer entangling
gates and increased reliance on noisy reset operations, while keeping the rest
of the heavy-hex construction unchanged.
"""

import stim  # type: ignore
from dataclasses import dataclass
from typing import List, Union, Dict, Set, Tuple, Callable
import gen
from collections import defaultdict
import warnings
from functools import lru_cache


@dataclass(frozen=True)
class Tile:
    measure: complex                 # ID of measurement qubit
    flags: Tuple[complex, ...]       # Two flag qubits (variable length)
    data: Tuple[complex, ...]        # 2 or 4 data qubits (variable length)
    phase: int = 0                   # Phase of the tile (0 or 1)


@lru_cache(maxsize=13)
def get_tiles(code_distance: int):
    x_tiles = FourBodyTileGenerator(code_distance).generate_4body_check()
    z_tiles = TwoBodyTileGenerator([], code_distance).generate_2body_checks()
    return x_tiles, z_tiles

def collect_qubit_sets(code_distance: int):
    x_tiles, z_tiles = get_tiles(code_distance)

    data = {q for tile in x_tiles for q in tile.data}
    measure = {tile.measure for tile in x_tiles}
    flags_z = {q for tile in z_tiles for q in tile.flags}
    flags_x = {q for tile in x_tiles for q in tile.flags}
    flags = flags_x | flags_z

    return data, measure, flags


def get_dataset(code_distance: int) -> List[complex]:
    x_tiles, z_tiles = get_tiles(code_distance)
    all_coords = {
        q
        for tile in x_tiles
        for q in tile.data + tile.flags + ((tile.measure,) if tile.measure is not None else ())
    }.union(
        q for tile in z_tiles for q in tile.flags
    )
    return sorted(all_coords, key=gen.complex_key)


class TwoBodyTileGenerator:
    """
    Prepare tiles for 2-body checks in the heavy hex code.

    This class generates tiles for 2-body checks based on a provided table of coordinates.
    The table is a list of complex numbers representing coordinates (where the real part is the row and the imaginary part is the column)
    of the 2-body measurement stabilizers (squares) to be retained from the previous step (after performing 4-body check measurements).

    A tile is a data structure containing:
        - flags: the qubits used to check the parity of the data qubits (these sit between data qubits),
        - data: the data qubits involved in the check,
        - phase: the phase for the 2-body check.

    The 'generate_2body_checks' function generates a list of Tile objects representing the required 2-body checks.
    """
    def __init__(self, table: List[complex], code_distance: int):
        # here table is a list of complex numbers representing coordinates (real is the row, imag is the column)
        # of the 2-body measurement formed stabilizers (squares) we want to keep from the previous step (when do 4-body check measurement)
        self.table = table
        self.code_distance = code_distance
        self.tile_group: List[Tile] = self._initialize_2body_check_tiles()
        
        

    def _invalid_coordinates(self) -> set:
        """
        Returns a set of flags that do not match the expected pattern.
        These coordinates are not selected for 2-body checks.
        """
        bad_flags = []
        for square in self.table:
            col = square.imag
            row = square.real
            if (col + row) % 2 == 0:
                raise ValueError(f"Invalid 2-body square coordinate in table: {square}. Expected odd sum of real and imaginary parts.")
            col_offset = col + 1
            row_offset = row + 0.5
            bad_flags.extend([complex(row_offset - 0.5, col_offset), complex(row_offset + 0.5, col_offset)])
        bad_flags_set = set(bad_flags)
        return bad_flags_set

    def _remove_invalid_coordinates(self) -> List[Tile]:
        """
        Removes invalid flag coordinates from the table that do not match the expected pattern.
        Removes the tiles that we do not want to use for 2-body checks.
        """
        bad_flags_set = self._invalid_coordinates()
        if bad_flags_set:
            return [
                tile for tile in self.tile_group
                if bad_flags_set.isdisjoint(tile.flags)
            ]
        return list(self.tile_group)

    def _initialize_2body_check_tiles(self) -> List[Tile]:
        """
        Generates 2-body check tiles based on the provided table of coordinates.

        Returns:
            List[Tile]: A list of Tile objects representing the selected 2-body checks.
        """
        tile_group_2body = []
        cols = self.code_distance - 1
        for i in range(self.code_distance):
            for j in range(cols):
                col_offset = j + 1
                data = [complex(i, col_offset-0.5), complex(i, col_offset + 0.5)]
                flags = [complex(i, col_offset)]
                phase = 0 if col_offset % 2 == 0 else 1
                tile_group_2body.append(Tile(measure=None, flags=flags, data=data, phase=phase)) # type: ignore
        if tile_group_2body:
            return tile_group_2body
        else:
            raise ValueError("Bad function generate_initial_2body_checks: no tiles generated.")

    def generate_2body_checks(self, with_table:bool=True) -> List[Tile]:
        if self.table and with_table:
            self.tile_group = self._remove_invalid_coordinates()
        else:
            self.tile_group = self._initialize_2body_check_tiles()

        return self.tile_group

class FourBodyTileGenerator:
    def __init__(self, code_distance: int):
        self.code_distance = code_distance
        self.tile_group: Set[Tile] = set()

    @staticmethod
    def _form_4body_tile(row_offset: float, col_offset: float, tile_type: int) -> Tile:
        """
        Creates a Tile object representing a tile in a heavy hex lattice, with configurable position, type, and phase.

        Args:
            row_offset (float): The row coordinate (center) of the tile.
            col_offset (float): The column coordinate (center) of the tile.
            tile_type (int): The type of tile to form.
                0: Full tile (covers both left and right halves).
                1: Half-right tile (covers only the right half).
                -1: Half-left tile (covers only the left half).

        Returns:
            Tile: A Tile object with the specified geometry and phase.

        Raises:
            ValueError: If tile_type is not one of 0, 1, or -1.
        """
        # tile_type: 0 = full, 1 = half-right, -1 = half-left
        if tile_type == 0:
            col_index = [col_offset - 0.5, col_offset + 0.5]
        elif tile_type == 1:
            col_index = [col_offset - 0.5]
        elif tile_type == -1:
            col_index = [col_offset + 0.5]
        else:
            raise ValueError(f"Invalid tile type: {tile_type}. Expected 0 (full), 1 (half-right), or -1 (half-left).")
        row_index = [row_offset - 0.5, row_offset + 0.5]
        data = tuple(complex(r, c) for r in row_index for c in col_index)
        flags = tuple(complex(r, col_offset) for r in row_index)
        measure = complex(row_offset, col_offset)
        phase = 0 if col_offset % 2 == 0 else 1 
        return Tile(measure, flags, data, phase)

    def generate_4body_check(self) -> Set[Tile]:
        """Generates tiles for 4-body checks in the heavy hex code."""
        m = 0
        rows = cols = self.code_distance - 1

        for i in range(rows):
            for j in range(cols):
                row_offset = i + 0.5
                if (i + j) % 2 == 0:
                    col_offset = j + 1
                    tile_type = 0  # full tile
                elif j == cols - 1:
                    col_offset = j + 2
                    tile_type = 1  # half-right tile
                elif j == 0:
                    col_offset = j
                    tile_type = -1  # half-left tile
                else:
                    continue  # skip non-tile positions

                self.tile_group.add(self._form_4body_tile(row_offset, col_offset, tile_type))
                m += 1

        return self.tile_group



class HeavyHexCircuitBuilder:
    def __init__(self, 
                 builder: gen.Builder, 
                 code_distance: int, 
                 step: int, 
                 tile_group_z: List[Tile], 
                 tile_group_x: Set[Tile],                  
                 recorded_measurement_keys: set):
        self.code_distance = code_distance
        self.tile_group_x = tile_group_x
        self.tile_group_z = tile_group_z
        self.step = step  # roundstep for the circuit
        self.recorded_measurement_keys = recorded_measurement_keys  # Set of measurement qubits to record
        self.coord_to_index: Dict[complex, int] = builder.q2i  # Mapping from complex coordinates to indices
        self.builder = builder  # Builder for the circuit
        self.data, self.measure, self.flags = collect_qubit_sets(code_distance)
        self.flags_z = {f for group in tile_group_z for f in group.flags}
        self.true_flags = self.flags - self.flags_z

    def _add_x_circuit(self):
        """
        Adds an optimized 4-Body X-check circuit for tiles.

        Operations across tiles in the same spatial region and phase are grouped together
        and applied in parallel to minimize circuit depth.
        """

        def operations(op: int, f: complex, d: complex, m: complex):
            if op == 0:
                self.builder.gate2("CX", [(f, d)])
            elif op == 1:
                self.builder.gate2("CX", [(m, f)])
            elif op == 2:
                self.builder.gate2("CX", [(f, d)])
            elif op == 3:
                self.builder.gate2("CX", [(m, f)])

        # Define selector functions for the 4 tile quadrants
        def selectors(data, flags, pos):
            if pos == "TL":
                return len(data) == 4, data[0].imag < flags[0].imag and data[0].real == flags[0].real, flags[0], data[0]
            if pos == "BL":
                return len(data) == 4, data[1].imag < flags[1].imag and data[1].real == flags[1].real, flags[1], data[1]
            if pos == "TR":
                if len(data) == 4:
                    valid = data[2].imag > flags[0].imag and data[2].real == flags[0].real
                    return True,valid, flags[0], data[2]
                else:
                    valid = data[0].imag > flags[0].imag and data[0].real == flags[0].real
                    return False, valid, flags[0], data[0]
            if pos == "BR":
                return len(data) == 4, data[-1].imag > flags[1].imag and data[-1].real == flags[1].real, flags[1], data[-1]

        schedule = [
        [("TL", 1)],
        [("TL", 2), ("BL", 1)],
        [("TL", 3), ("BL", 2), ("TR", 3)],
        [("BL", 3), ("TR", 0), ("BR", 3)],
        [("TR", 1), ("BR", 0)],
        [("BR", 1)],
        ]
        # schedule = [
        # [("TL", 1)],
        # [("TL", 2), ("BL", 1)],
        # [("TL", 3), ("BL", 2)],
        # [("BL", 3), ("TR", 0)],
        # [("TR", 1), ("BR", 0)],
        # [("TR", 2), ("BR", 1)],
        # [("TR", 3), ("BR", 2)],
        # [("BR", 3)],
        # ]
        
        for moment in schedule:
            for region, op in moment:
                for tile in self.tile_group_x:
                    if len(tile.data) == 2 or (len(tile.data) == 4 and op != 3):
                        m = tile.measure
                        data = sorted(tile.data, key=lambda a: (a.imag, a.real))
                        flags = sorted(tile.flags, key=lambda a: a.real)
                        result = selectors(data, flags, region)
                        if result is not None:
                            full, valid, f, d = result
                            if valid:
                                operations(op, f, d, m)

            self.builder.shift_coords(dt=1)
            self.builder.tick()


    def _add_z_circuit(self):
        """
        Adds a 2-body Z-check circuit for the given tile.

        This method appends controlled-X (CX) gates to the circuit for a Z-check operation.
        It expects the tile to have exactly two data qubits. For each data qubit, a CX gate is
        added with the data qubit as control and the flag qubit as target.

        Args:
            tile (Tile): The tile object containing data and flag qubit coordinates.

        Raises:
            ValueError: If the number of data qubits in the tile is not exactly two.
        """
        # Define selector functions for the 4 tile quadrants
        def selectors(data, flags, pos):
            if pos == "L":
                return data[0].imag < flags[0].imag and data[0].real == flags[0].real, flags[0], data[0]
            if pos == "R":
                return data[1].imag > flags[0].imag and data[1].real == flags[0].real, flags[0], data[1]

        regions = ["L", "R"]
        for region in regions:
            for tile in self.tile_group_z:
                data_length = len(tile.data)
                if data_length != 2:
                    raise ValueError(f"Invalid number of data qubits for Z-check: {data_length}")
                data = sorted(tile.data, key=lambda a: a.imag)
                flags = sorted(tile.flags, key=lambda a: a.real)
                result = selectors(data, flags, region)
                if result is not None:
                    valid, f, d = result
                    if valid:
                        self.builder.gate2("CX", [(d, f)])
            self.builder.shift_coords(dt=1)
            self.builder.tick()




    def apply_circuit_stage(self, stage: str, basis: str) -> None:
        """
        Applies a specific initialization or measurement stage to the circuit.

        Args:
            stage (str): One of {"init", "R", "M"}
                - "init": Initializes the circuit RX on data qubits, R on measure and flag
                - "R": Resets the measurement and flag qubits to the |0> state. R on measure and flag only
                - "M": Measures all measurement and flag qubits. M on measure and flag only
        Raises:
            ValueError: If the stage is not one of "init", "R", or "M".
        """

        # Apply the operator to all measurement and flag qubits
        if stage == "init":
            self.builder.gate("R", self.flags)
            self.builder.gate("R" if basis == 'Z' else "RX", self.data)
            self.builder.gate("RX", self.measure)
            self.builder.tick()


        elif stage == "R":
            self.builder.gate("RX", self.measure)
            self.builder.gate("R", self.flags)
            self.builder.shift_coords(dt=1)
            self.builder.tick()

        elif stage == "M":
            # M on measure qubits with tracking            
            self.builder.measure(
                self.measure,
                basis='X',
                tracker_key=lambda q: q,
                save_layer=f'x_round_{self.step}',
            )
            
            self.builder.measure(
                self.flags_z,
                basis='Z',
                tracker_key=lambda q: q,
                save_layer=f'z_round_{self.step}',
            )

            self.builder.measure(
                self.true_flags,
                basis='Z',
                tracker_key=lambda q: q,
                save_layer=f'f_round_{self.step}',
            )

            # Save keys
            for m in self.measure:
                key = gen.AtLayer(m, f'x_round_{self.step}')
                self.recorded_measurement_keys.add(key)
            for f in self.flags_z:
                key = gen.AtLayer(f, f'z_round_{self.step}')
                self.recorded_measurement_keys.add(key)
            for f in self.true_flags:
                key = gen.AtLayer(f, f'f_round_{self.step}')
                self.recorded_measurement_keys.add(key)
        else:
            raise ValueError(f"Invalid stage: {stage}. Expected 'init', 'R', or 'M'.")
    
    def get_recorded_measurement_keys(self) -> set:
        """
        Returns the set of recorded measurement keys for the current circuit step.
        
        This set contains keys for all measurement and flag qubits that are recorded
        during the circuit generation process.
        
        Returns:
            set: A set of recorded measurement keys.
        """
        return self.recorded_measurement_keys


    def generate_round_circuit(self, basis: str, is_end: bool = False):
        """
        Generates a round of the heavy hex circuit. Performing one X-check and one Z-check is one round.

        This method adds X-checks and Z-checks for each tile in the respective tile groups,
        grouped by phase, and adds measurement operations for all measurement and flag qubits.
        """
        # Reset the circuit stage
        if self.step == 0:
            self.apply_circuit_stage("init", basis=basis)

        else:
            self.apply_circuit_stage("R", basis=basis)         

        if basis == 'Z':
            self._add_z_circuit()
            self._add_x_circuit()
        else:
            self._add_x_circuit()
            self._add_z_circuit()

        # Apply the circuit stage
        self.apply_circuit_stage("M", basis=basis)
        if is_end:
            self.do_end_measures(basis=basis)

        return self.builder.circuit

    def do_end_measures(self, basis: str):
        layer = 'x_round_end' if basis == 'X' else 'z_round_end'
        self.builder.measure(self.data, basis=basis, save_layer=layer)

        self.recorded_measurement_keys.update(gen.AtLayer(q, layer) for q in self.data)

        if basis == 'X':
            qubits = [gen.AtLayer(q, layer) for q in self.data if q.real == 0]
        elif basis == 'Z':
            qubits = [gen.AtLayer(q, layer) for q in self.data if q.imag == 0.5]
        else:
            warnings.warn(f"Unknown basis {basis} for end measures. Using X basis by default.")
            qubits = [gen.AtLayer(q, 'x_round_end') for q in self.data if q.real == 0]

        self.builder.obs_include(qubits, obs_index=0)
    


class ConstructDetectors:
    def __init__( self, 
                  builder: gen.Builder,
                  code_distance: int,
                  step: int,
                  tile_group_z: List[Tile],
                  tile_group_x: Set[Tile],
                  recorded_measurement_keys: set 
                  ):
        self.builder = builder
        self.code_distance = code_distance
        self.tile_group_z = tile_group_z
        self.tile_group_x = tile_group_x
        self.recorded_measurement_keys = recorded_measurement_keys
        self.step = step
        _, _, self.flags = collect_qubit_sets(code_distance)
        self.true_flags = self.flags - {f for group in tile_group_z for f in group.flags}

    def init_detectors(self, basis: str):
        group = self.tile_group_x if basis == 'X' else self.tile_group_z
        for tile in group:
            keys = {gen.AtLayer(tile.measure if basis == 'X' else tile.flags[0], f'{basis.lower()}_round_0')}
            if keys:
                min_pos = tile.measure if basis == 'X' else min(tile.flags, key=gen.complex_key)
                self.builder.detector(
                    keys,
                    pos=min_pos,
                    extra_coords=[1] if basis == 'X' else [2]
                )
        self.flag_detectors()  # Generate flag detectors
        # Shift coordinates and tick the builder to finalize the initialization
        self.builder.shift_coords(dt=1)
        self.builder.tick()
    
    def flag_detectors(self):
        """Generates detectors for the flag qubits (not participating in Z checks)."""

        for flag in self.true_flags:
            key = gen.AtLayer(flag, f'f_round_{self.step}')
            # self.builder.detector(keys=[key], pos=flag, extra_coords=[2])
            if key in self.recorded_measurement_keys:
                self.builder.detector(keys=[key], pos=flag, extra_coords=[2])
            else:
                raise ValueError(f"flag detector not build")

    def end_detectors(self, table: List[complex], basis: str, org_tile_group_z: List[Tile] = []):  # should be the table in current time step
        """Generates detectors for the end of the circuit."""
        if self.step == 0:
            raise ValueError("Step must be greater than 0 to generate end detectors.")
        if basis == 'X':
            self.general_x_detectors(table, basis, is_end=True) 
        elif basis == 'Z':
            self.general_z_detectors(basis, is_end=True, org_tile_group_z=org_tile_group_z)
        else:
            raise ValueError(f"Unknown basis {basis} for end detectors. Use 'X' or 'Z'.")
        self.flag_detectors()  # Generate flag detectors
        
    def general_x_detectors(self, table: List[complex], basis: str, is_end: bool = False):  # should be the table in previous time step
        """Performs the X detectors in the circuit. The table is a list of complex numbers representing coordinates (where the real part is the row and the imaginary part is the column)
        of the 2-body measurement formed stabilizers (squares) we want to keep from the previous step (when do 4-body check measurement).
        """
        if self.step == 0:
            raise ValueError("Step must be greater than 0 to generate X detectors otherwise use function init_x_detectors.")
        # Group tiles and squares by row (real part)
        tiles_by_row = defaultdict(list)
        squares_by_row = defaultdict(list)
        for tile in self.tile_group_x:
            row = tile.measure.real - 0.5  # Adjust row to match square row
            tiles_by_row[row].append(tile)
        
        for square in table:
            row = square.real
            squares_by_row[row].append(square.imag)
            
        # For each row, process tiles and separate them by square-imag dividers
        for row in sorted(tiles_by_row.keys()):
            tiles = sorted(tiles_by_row[row], key=lambda t: t.measure.imag)
            dividers = sorted(squares_by_row.get(row, []))  # squares in the same row as tiles behave as dividers

            # Build the divider bins (with endpoints)
            groups: list[list[Tile]] = []
            current_group: list[Tile] = []
            divider_idx = 0

            for tile in tiles:
                tile_col = tile.measure.imag # get the column of the measurement qubit as the col of the tile
                # divider works
                while divider_idx < len(dividers) and tile_col > dividers[divider_idx]:
                    # Finish current group before divider
                    if current_group:
                        groups.append(current_group)
                        current_group = []
                    divider_idx += 1
                # Add tile to the current group before the next divider works
                current_group.append(tile)

            # Add the final group
            if current_group:
                groups.append(current_group)

            # Now do something with each group
            for group in groups:
                keys = set()
                keys_end = set()
                for tile in group:
                    if not is_end:
                        keys.add(gen.AtLayer(tile.measure, f'x_round_{self.step}')) #current round step
                        if gen.AtLayer(tile.measure, f'x_round_{self.step - 1}') in self.recorded_measurement_keys:
                            keys.add(gen.AtLayer(tile.measure, f'x_round_{self.step - 1}')) # last round step
                        else:
                            raise ValueError(f"Measurement key {gen.AtLayer(tile.measure, f'x_round_{self.step - 1}')} not found in recorded keys.")
                    elif is_end and basis == 'X':  # set up the detectors for the end step
                        for data in tile.data:
                            keys_end.add(gen.AtLayer(data, f'x_round_end')) #current round step
                        if gen.AtLayer(tile.measure, f'x_round_{self.step}') in self.recorded_measurement_keys:
                            keys_end.add(gen.AtLayer(tile.measure, f'x_round_{self.step}'))

                if keys:
                    min_pos = min((tile.measure for tile in group), key=gen.complex_key)
                    self.builder.detector(
                        keys,
                        pos=min_pos,
                        extra_coords=[1]
                    )
                if keys_end:
                    min_pos = min((tile.measure for tile in group), key=gen.complex_key)
                    self.builder.detector(
                        keys_end,
                        pos=min_pos,
                        extra_coords=[1]
                    )

    def seperate_tiles_z(self,tile_group_z: List[Tile]) -> List[List[Tile]]:
        tiles_by_rc = defaultdict(list)
        for tile in tile_group_z:
            flag = tile.flags[0]
            row = int(flag.real)
            col = int(flag.imag - 1)  # match square column
            keys =[(row, col) , (row - 1, col)]
            if row == 0:
                tiles_by_rc[keys[0]].append(tile)
            elif row == self.code_distance - 1:
                tiles_by_rc[keys[1]].append(tile)
            else:
                tiles_by_rc[keys[0]].append(tile)
                tiles_by_rc[keys[1]].append(tile)

        groups = []

        for (row, col), tiles in tiles_by_rc.items():
            if (row + col) % 2 != 0:
                groups.append(tiles)  # accept all tiles in this group
            else:
                edge_tiles = [tile for tile in tiles if tile.flags[0].real in {0, self.code_distance - 1}]
                if edge_tiles:
                    groups.append(edge_tiles)
        return groups


    def general_z_detectors(self,basis: str, is_end: bool = False, org_tile_group_z: List[Tile]=[]):
        if self.step == 0:
            raise ValueError("Step must be greater than 0 to generate Z detectors.")
        if basis == 'Z' and is_end:
            tile_group = org_tile_group_z
        else:
            tile_group = self.tile_group_z
        groups = self.seperate_tiles_z(tile_group)
        
        if not is_end:
            check_steps = [self.step - 1, self.step - 2] if self.step > 1 else [self.step - 1]
        elif is_end and basis == 'Z':
            check_steps = [self.step, self.step - 1]

        for group in groups:
            flags = {tile.flags[0] for tile in group}
            prev_step0 = {gen.AtLayer(f, f"z_round_{check_steps[0]}") for f in flags}
            prev_step1 = (
                            {gen.AtLayer(f, f"z_round_{check_steps[1]}") for f in flags}
                            if len(check_steps) > 1
                            else set()
                        )
            recorded = self.recorded_measurement_keys
            have_0   = prev_step0 <= recorded           # all flags present at step - 1
            have_1   = prev_step1 and prev_step1 <= recorded  # all flags present at step - 2 
            if not is_end:
                keys  = {gen.AtLayer(flag, f"z_round_{self.step}") for flag in flags}
                if  have_0:
                    keys.update(prev_step0)
                elif have_1 and not have_0:
                    keys.update(prev_step1)
                else:
                    keys = keys if basis == "Z" else set()
                    warnings.warn(f"No previous step keys found for z-detectors at tail {group} at timestep {self.step} ")
                if keys:
                    min_pos = min(flags, key=gen.complex_key)
                    self.builder.detector(keys, pos=min_pos, extra_coords=[2])
            elif is_end and basis == 'Z':
                keys_end = {gen.AtLayer(d, 'z_round_end') for tile in group for d in tile.data}
                data = {d for tile in group for d in tile.data}
                if have_0:
                    keys_end.update(prev_step0)
                elif have_1 and not have_0:
                    keys_end.update(prev_step1)
                else:
                    keys_end = set()
                    warnings.warn(f"No previous step keys found for z-detectors at tail {group} at timestep {self.step} ")
                if keys_end:
                    min_pos = min(data, key=gen.complex_key)
                    self.builder.detector(keys_end, pos=min_pos, extra_coords=[2])

    def detector_generator(self , table: List,basis: str, is_end: bool = False):
        self.general_x_detectors(table, basis, is_end)
        self.general_z_detectors(basis, is_end)
        self.flag_detectors()  # Generate flag detectors
        self.builder.shift_coords(dt=1)
        self.builder.tick()

def make_heavy_hex_circuit(
    *,
    table: List,
    diameter: int,
    rounds: int,
    basis: str,
) -> stim.Circuit:
    all_coords = get_dataset(diameter)
    builder = gen.Builder.for_qubits(all_coords)
    table = filter_table_by_diameter(table, diameter)
    recorded_measurements: set[str] = set()

    fourbody = FourBodyTileGenerator(diameter)
    fourbody_tiles = fourbody.generate_4body_check()

    for step in range(rounds):
        is_end = step == rounds - 1
        twobody = TwoBodyTileGenerator(table=table[step % 2], code_distance=diameter)
        twobody_tiles = twobody.generate_2body_checks()
        
        circuitbuilder = HeavyHexCircuitBuilder(
            builder,
            diameter,
            step,
            twobody_tiles,
            fourbody_tiles,
            recorded_measurements,
        )
        circuitbuilder.generate_round_circuit(basis, is_end=is_end)
        detector = ConstructDetectors(
            builder,
            diameter,
            step,
            twobody_tiles,
            fourbody_tiles,
            recorded_measurements,
        )
        if step == 0:
            detector.init_detectors(basis)
        else:
            if basis == 'X':
                detector.detector_generator(table[(step - 1) % 2], basis)
            elif basis == 'Z':
                detector.detector_generator(table[(step) % 2], basis)
        if is_end:
            detector.end_detectors(
                table=table[step % 2],
                basis=basis,
                org_tile_group_z=twobody.generate_2body_checks(with_table=False),
            )
    return builder.circuit

def filter_table_by_diameter(table: List[List[complex]], diameter: int) -> List[List[complex]]:
    """
    Filters the input table so that only coordinates with both real (row) and imaginary (col)
    parts less than (diameter - 2) are kept.

    Args:
        table (List[List[complex]]): The input table, a list of lists of complex numbers.
        diameter (int): The code distance.

    Returns:
        List[List[complex]]: The filtered table.
    """
    limit = diameter - 1
    filtered = []
    for step_list in table:
        step_filtered = [
            sq for sq in step_list
            if sq.real < limit and sq.imag < limit
        ]
        filtered.append(step_filtered)
    return filtered
