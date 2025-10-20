import stim
import sinter
import re
import gen
import utils
import ast
import heavyhex._heavy_hex as hh
from beliefmatching import BeliefMatching
import pymatching
import beliefmatching.belief_matching as bm

from typing import Iterable, Tuple, Set, List, Union


# table_file = []
table_file = "out/plot/t1-th/table.txt"
diameter = 3
rounds = 4

if table_file:
    with open(table_file) as f:
        table = ast.literal_eval(f.read())
        table = hh.filter_table_by_diameter(table, diameter)
else:
    table = [[], []]
    
circuit = hh.make_heavy_hex_circuit(
        table=table,
        diameter=diameter,
        rounds=rounds,
        basis="Z"
    )
circuit = gen.NoiseModel.uniform_depolarizing(1e-3).noisy_circuit(circuit)

flags = []
for step in range(2):
    twobody = hh.TwoBodyTileGenerator(table=table[step], code_distance=diameter)
    tile_group_z = twobody.generate_2body_checks()
    data, measure, flags_all = hh.collect_qubit_sets(diameter)
    flags_z = {f for group in tile_group_z for f in group.flags}
    true_flags = flags_all - flags_z
    flags.append(true_flags)

with open("test/result/flags.txt", "w") as f:
    f.write(repr(flags))

def dem(circuit):
    dem = circuit.detector_error_model(decompose_errors=True)
    return dem

def _to_dem_str(circuit_or_dem: Union[stim.Circuit, stim.DetectorErrorModel, str], *, decompose=True) -> str:
    if isinstance(circuit_or_dem, stim.Circuit):
        return str(circuit_or_dem.detector_error_model(decompose_errors=decompose))
    if isinstance(circuit_or_dem, stim.DetectorErrorModel):
        return str(circuit_or_dem)
    if isinstance(circuit_or_dem, str):
        return circuit_or_dem
    raise TypeError("Expected stim.Circuit, stim.DetectorErrorModel, or str")

def _clean_tokens(tokens, remove_ids_set):
    # 1) remove targeted Dk tokens
    kept = []
    for tok in tokens:
        m = re.fullmatch(r"D(\d+)", tok)
        if m and int(m.group(1)) in remove_ids_set:
            continue
        kept.append(tok)

    # 2) normalize separators: remove duplicate ^, and trim leading/trailing ^
    # collapse consecutive ^
    norm = []
    for tok in kept:
        if tok == "^" and (not norm or norm[-1] == "^"):
            continue
        norm.append(tok)
    # trim leading/trailing ^
    if norm and norm[0] == "^":
        norm = norm[1:]
    if norm and norm[-1] == "^":
        norm = norm[:-1]

    # 3) if only ^ or empty => empty
    non_sep = [t for t in norm if t != "^"]
    if not non_sep:
        return None  # means: drop the line

    return norm

def remove_detectors_from_dem(
    circuit_or_dem: Union[stim.Circuit, stim.DetectorErrorModel, str],
    remove_ids: Iterable[int],
    *,
    decompose_errors: bool = True,
) -> stim.DetectorErrorModel:
    dem_str = _to_dem_str(circuit_or_dem, decompose=decompose_errors)
    remove_ids_set = set(int(x) for x in remove_ids)

    new_lines = []
    for line in dem_str.splitlines():
        if line.startswith("error("):
            # Split "error(prob) ..." into prefix and tokens
            # find the first ')' that closes error(...)
            i = line.find(")")
            if i == -1:
                # malformed; just keep as-is
                new_lines.append(line)
                continue
            prefix = line[: i + 1]  # includes "error(prob)"
            rest = line[i + 1 :].strip()

            if not rest:
                # already no targets; drop
                continue

            # Tokenize by spaces; DEM uses spaces between targets and '^'
            tokens = rest.split()

            cleaned = _clean_tokens(tokens, remove_ids_set)
            if cleaned is None:
                # no targets remain => drop line
                continue

            new_lines.append(prefix + " " + " ".join(cleaned))
        else:
            # keep non-error lines verbatim (detectors, shifts, etc.)
            new_lines.append(line)

    return stim.DetectorErrorModel("\n".join(new_lines))

DetectorCoord = Tuple[float, float]

_DETECTOR_RE = re.compile(
    r"""^detector\(\s*
        (?P<x>-?\d+(?:\.\d+)?)\s*,\s*
        (?P<y>-?\d+(?:\.\d+)?)
        (?:\s*,\s*-?\d+(?:\.\d+)?\s*,\s*-?\d+(?:\.\d+)?)?
        \)\s+D(?P<idx>\d+)\s*$""",
    re.X
)

def _float_pair_equal(a: DetectorCoord, b: DetectorCoord, tol: float = 1e-9) -> bool:
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol

def _contains_coord(coord_set: Set[DetectorCoord], xy: DetectorCoord, tol: float = 1e-9) -> bool:
    # Coordinate sets are small; tolerant membership check
    return any(_float_pair_equal(xy, c, tol) for c in coord_set)

def select_detectors_by_coords(
    dem_text: Union[str, "stim.DetectorErrorModel"],
    coords_even: Iterable[complex],
    coords_odd: Iterable[complex],
    *,
    tol: float = 1e-9,
) -> List[int]:
    """
    Parse a Stim DEM (string or stim.DetectorErrorModel), infer timestep t for each
    detector block via shift_detectors(..., 1), and collect detector indices Dk whose
    (x,y) match coords_even when t%2==0 or coords_odd when t%2==1.

    - coords_even / coords_odd: iterables of complex numbers; we use (c.real, c.imag) as (x,y).
    - Matching is tolerant up to `tol` to handle float formatting.
    """
    if hasattr(dem_text, "__str__"):  # stim.DetectorErrorModel object or similar
        dem_text = str(dem_text)

    # Build coordinate sets as (x, y) float pairs
    even_set: Set[DetectorCoord] = {(complex(c).real, complex(c).imag) for c in coords_even}
    odd_set: Set[DetectorCoord]  = {(complex(c).real, complex(c).imag) for c in coords_odd}

    selected: List[int] = []

    # We detect timesteps by noticing when we encounter detectors after at least one shift_detectors.
    t = -1
    shift_index = 0
    pending_shift = False

    for line in dem_text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith("shift_detectors"):
            # Any number of these may occur between detector groups; they advance time for future detectors.
            pending_shift = True
            shift_index += 1
            continue

        if line.startswith("error("):
            continue

        m = _DETECTOR_RE.match(line)
        if m:
            # If we just passed through a group of shifts before detectors, we’re at a new time step.
            if pending_shift and shift_index > 1:
                t += 1
                pending_shift = False
                shift_index = 0
            x = float(m.group("x"))
            y = float(m.group("y"))
            idx = int(m.group("idx"))

            target_set = even_set if (t % 2 == 0) else odd_set
            if _contains_coord(target_set, (x, y), tol=tol):
                selected.append(idx)

            continue

    return selected

# indices = select_detectors_by_coords(dem(circuit), coords_even=flags[0], coords_odd=flags[1])
# print(indices)
# flags = [4, 5, 6, 7, 8, 9, 16, 17, 18, 19, 20, 21, 28, 29, 30, 31, 32, 33, 40, 41, 42, 43, 44, 45, 50, 51, 52, 53, 54, 55]

check_matrix = bm.detector_error_model_to_check_matrices(dem(circuit))

with open("test/result/check_matrix.txt", "w") as f:
    f.write(repr(check_matrix))

# decoder = BeliefMatching(dem(circuit), max_bp_iters=20) 

# circuit = stim.Circuit.from_file("out_t1-flag/circuits/r=12,d=3,p=0.001,noise=IBM1907,c=heavy_hex,q=23,b=X,g=all.stim")
utils.export_pdf(circuit, diagram_type='timeline-svg', pdf_path='test/result/timeline.pdf')
with open("test/result/detector_error_model.txt", "w") as f:
    f.write(repr(dem(circuit)))
# matching_corr = pymatching.Matching.from_detector_error_model(remove_detectors_from_dem(circuit, flags), enable_correlations=True)
