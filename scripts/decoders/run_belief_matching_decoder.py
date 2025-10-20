import argparse
import glob
import os
import csv
import stim
import sinter
from beliefmatching import BeliefMatchingSinterDecoder

def main():
    parser = argparse.ArgumentParser(description="Run BeliefMatching decoder on .stim circuits using sinter.")
    parser.add_argument('--circuits',
                        nargs='+',
                        required=True,
                        help='Circuit files to sample from and decode.\n'
                             'This parameter can be given multiple arguments.')
    parser.add_argument('--max-shots', type=int, required=True, help='Maximum number of shots to run.')
    parser.add_argument('--max-errors', type=int, required=True, help='Maximum number of errors to stop at.')
    parser.add_argument('--processes', type=int, default=1, help='Number of parallel processes.')
    parser.add_argument('--max_bp_iters', type=int, default=20, help='Maximum number of belief propagation iterations.')
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

    def generate_tasks(circuit_files):
        for filename in circuit_files:
            circuit = stim.Circuit.from_file(filename)
            # Use sinter.comma_separated_key_values to parse metadata from filename
            json_metadata = sinter.comma_separated_key_values(filename)
            yield sinter.Task(
                circuit=circuit,
                json_metadata=json_metadata,
            )

    samples = sinter.collect(
        num_workers=args.processes,
        tasks=generate_tasks(circuit_files),
        save_resume_filepath=args.save_resume,
        max_shots=args.max_shots,
        max_errors=args.max_errors,
        decoders=["beliefmatching"],
        custom_decoders={"beliefmatching": BeliefMatchingSinterDecoder(max_bp_iters=args.max_bp_iters)},
    )

if __name__ == '__main__':
    main()
