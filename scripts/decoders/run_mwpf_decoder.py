"""
MWPF decoding runner for .stim circuits using sinter.

This script wires up the `mwpf` decoder (Minimum-Weight Parity Factor, HyperBlossom/Hyperion)
to sinter's collection pipeline. It selects a per-task `cluster_node_limit` based on the code
distance `d` embedded in each circuit's filename (parsed via `sinter.comma_separated_key_values`).

Usage examples:
    python mwpfdecoder.py \
    --circuits "$circuits_dir"/*.stim \
    --save-resume "stats.csv" \
    --processes 4 \
    --max-shots 10_000_000 \
    --max-errors 1000

Notes:
- The heuristic in `cluster_limit_from_d` is conservative for d>=4. Tune per noise model if needed.
- If debugging multiprocessing issues, try `--processes 1` to run single-process.
"""
import argparse
import glob
import os
import stim
import sinter
from mwpf import SinterMWPFDecoder

def cluster_limit_from_d(d: int) -> int:
    """
    Map code distance d -> MWPF `cluster_node_limit`.

    Rationale:
    - Clusters in circuit-level noise tend to grow super-linearly with distance.
    - These presets aim to avoid over-pruning for d>=4 while keeping runtime sane.
    - For d>7 we fall back to a gentle ~d^2 growth (adjust if accuracy suffers).

    Returns:
        int: suggested per-cluster node limit for the MWPF decoder.
    """
    preset = {2: 40, 3: 60, 4: 90, 5: 120, 6: 160, 7: 200}
    if d in preset:
        return preset[d]
    # fallback heuristic for d > 7
    return max(200, int(18 * (d ** 2)))  # gentle ~d^2 growth


def main():
    parser = argparse.ArgumentParser(description="Run SinterMWPFDecoder on .stim circuits using sinter.")
    parser.add_argument('--circuits',
                        nargs='+',
                        required=True,
                        help='Circuit files to sample from and decode.\n'
                             'This parameter can be given multiple arguments.')
    parser.add_argument('--max-shots', type=int, required=True, help='Maximum number of shots to run.')
    parser.add_argument('--max-errors', type=int, required=True, help='Maximum number of errors to stop at.')
    parser.add_argument('--processes', type=int, default=1, help='Number of parallel processes.')
    parser.add_argument('--save-resume', type=str, required=True, help='Path to save or resume sinter results.')
    parser.add_argument('--out', type=str, default=None, help='Optional CSV output file.')
    

    args = parser.parse_args()

        # Expand multiple patterns/files passed to --circuits
    circuit_files_set = set()
    for pat in args.circuits:
        matched = glob.glob(pat)
        if matched:
            circuit_files_set.update(matched)
        elif os.path.isfile(pat):
            circuit_files_set.add(pat)
    circuit_files = sorted(circuit_files_set)
    if not circuit_files:
        raise SystemExit("No .stim files found for the given --circuits arguments.")

    tasks = []
    needed_limits = set()

    for filename in circuit_files:
        circuit = stim.Circuit.from_file(filename)
        json_metadata = sinter.comma_separated_key_values(filename)
        d = int(json_metadata.get("d", 0))
        limit = cluster_limit_from_d(d)
        needed_limits.add(limit)

        tasks.append(
            sinter.Task(
                circuit=circuit,
                json_metadata=json_metadata,
            )
        )
    
    # one MWPF instance per unique limit (top-level class; picklable)
    custom_decoders = {
        f"mwpf": SinterMWPFDecoder(cluster_node_limit=L)
        for L in needed_limits
    }

    samples = sinter.collect(
        num_workers=args.processes,
        tasks=tasks,
        save_resume_filepath=args.save_resume,
        max_shots=args.max_shots,
        max_errors=args.max_errors,
        decoders=["mwpf"],
        custom_decoders=custom_decoders,
    )

if __name__ == '__main__':
    main()
