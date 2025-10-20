#!/usr/bin/env python3
"""
Heavy-Hex Circuit Generator

This tool generates Stim circuit files for Heavy-Hex quantum error correction codes
with configurable noise models. It supports both uniform depolarizing noise and
correlated mixing noise models for multi-dimensional threshold studies.

The generated circuits can be used with decoders to estimate logical error rates
and compute error correction thresholds.

Usage Examples
--------------
Generate circuits with uniform noise:
    python tools/generate_circuit.py \\
        --diameter 5 \\
        --rounds auto \\
        --noise_model uniform \\
        --noise_strength 0.001 \\
        --b Z

Generate circuits with correlated mixing noise:
    python tools/generate_circuit.py \\
        --diameter 7 \\
        --rounds 28 \\
        --noise_model corr_mixing \\
        --noise_strength 0.001 \\
        --base_prob 0.0002 \\
        --sweep_slot idle \\
        --bias_slot_1 D2 \\
        --bias_slot_1_prob 0.0007 \\
        --bias_slot_2 M \\
        --bias_slot_2_prob 0.005 \\
        --table_file data/breaker_table.txt

Author: Heavy-Hex Research Team
"""
import heavyhex._heavy_hex as hh
import gen
import argparse
import os
import ast


def main():
    """
    Main entry point for circuit generation.

    Parses command-line arguments, generates a Heavy-Hex circuit with the specified
    parameters and noise model, and saves it to a .stim file with metadata-rich naming.
    """
    parser = argparse.ArgumentParser(
        description="Generate Heavy-Hex quantum error correction circuits with noise",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Noise Models:
  uniform       - Uniform depolarizing noise with single parameter p
  corr_mixing   - Correlated mixing noise with slot-specific error rates
                  (requires sweep_slot, bias_slot_1, bias_slot_2 to be distinct)

Output:
  Circuits are saved to --out_dir with filenames encoding all parameters.
  Example: r=20,d=5,p=0.001,noise=uniform,c=heavy_hex,b=Z,g=all.stim
        """
    )
    parser.add_argument("--out_dir", type=str, default="data/circuits",
                        help="Output directory for generated circuits (default: data/circuits)")
    parser.add_argument("--diameter", type=int, required=True,
                        help="Code distance (grid diameter) for the Heavy-Hex lattice")
    parser.add_argument("--rounds", type=str, required=True,
                        help="Number of QEC rounds (int or 'auto' for 4*diameter)")
    parser.add_argument("--noise_model", type=str, default="uniform",
                        choices=["uniform", "corr_mixing"],
                        help="Noise model: 'uniform' or 'corr_mixing'")
    parser.add_argument("--noise_strength", type=float, default=0.001,
                        help="Primary noise parameter (swept in threshold studies)")
    parser.add_argument("--style", type=str, default="heavy_hex",
                        help="Code style identifier (for filename metadata)")
    parser.add_argument("--b", "--basis", type=str, default="Z", dest="b",
                        choices=["X", "Z"],
                        help="Logical basis for encoding (X or Z)")

    # Heavy-Hex breaker table configuration
    parser.add_argument("--table_file", type=str,
                        help="Path to breaker table file (Python literal with complex coordinates)")

    # Correlated mixing noise model parameters
    parser.add_argument("--base_prob", type=float, default=0.001,
                        help="Base error probability for unspecified slots (corr_mixing only)")
    parser.add_argument("--sweep_slot", type=str, default="D2",
                        choices=["D1", "D2", "M", "idle"],
                        help="Slot name for swept channel (corr_mixing only)")
    parser.add_argument("--bias_slot_1", type=str, default="D1",
                        choices=["D1", "D2", "M", "idle"],
                        help="First biased slot for 3D threshold study (corr_mixing only)")
    parser.add_argument("--bias_slot_1_prob", type=float, default=0.0,
                        help="Error probability for bias_slot_1 (corr_mixing only)")
    parser.add_argument("--bias_slot_2", type=str, default="M",
                        choices=["D1", "D2", "M", "idle"],
                        help="Second biased slot for 3D threshold study (corr_mixing only)")
    parser.add_argument("--bias_slot_2_prob", type=float, default=0.0,
                        help="Error probability for bias_slot_2 (corr_mixing only)")

    args = parser.parse_args()

    # Determine number of QEC rounds
    if args.rounds == "auto":
        rounds = 4 * args.diameter  # Standard convention: rounds = 4d
    else:
        rounds = int(args.rounds)

    # Load breaker table if specified (for custom Heavy-Hex configurations)
    if args.table_file:
        with open(args.table_file) as f:
            table = ast.literal_eval(f.read())
    else:
        table = [[], []]  # Empty table = standard Heavy-Hex layout

    # Ensure output directory exists
    os.makedirs(args.out_dir, exist_ok=True)

    # Generate the noiseless Heavy-Hex circuit
    circuit = hh.make_heavy_hex_circuit(
        table=table,
        diameter=args.diameter,
        rounds=rounds,
        basis=args.b
    )

    # Apply noise model to the circuit
    if args.noise_model == "uniform":
        # Uniform depolarizing noise: single parameter applies to all gates
        circuit = gen.NoiseModel.uniform_depolarizing(args.noise_strength).noisy_circuit(circuit)
        filename = (
            f"{args.out_dir}/r={rounds},d={args.diameter},"
            f"p={args.noise_strength},noise={args.noise_model},"
            f"c={args.style},b={args.b},g=all.stim"
        )

    elif args.noise_model == "corr_mixing":
        # Correlated mixing noise: slot-specific error rates for multi-dimensional studies
        # Validate that the three slots are distinct (required for 3D threshold studies)
        if len({args.sweep_slot, args.bias_slot_1, args.bias_slot_2}) != 3:
            raise ValueError(
                f"corr_mixing requires three distinct slots: "
                f"sweep_slot='{args.sweep_slot}', "
                f"bias_slot_1='{args.bias_slot_1}', "
                f"bias_slot_2='{args.bias_slot_2}'."
            )

        circuit = gen.NoiseModel.corr_mixing(
            base_prob=args.base_prob,
            sweep_slot=args.sweep_slot,
            sweep_prob=args.noise_strength,
            bias_slot_1=args.bias_slot_1,
            bias_slot_1_prob=args.bias_slot_1_prob,
            bias_slot_2=args.bias_slot_2,
            bias_slot_2_prob=args.bias_slot_2_prob
        ).noisy_circuit(circuit)

        # Filename encodes all noise parameters for tracking
        filename = (
            f"{args.out_dir}/r={rounds},d={args.diameter},"
            f"p={args.noise_strength},p_base={args.base_prob},"
            f"noise={args.noise_model},c={args.style},"
            f"sweep={args.sweep_slot},"
            f"b1={args.bias_slot_1},b1p={args.bias_slot_1_prob},"
            f"b2={args.bias_slot_2},b2p={args.bias_slot_2_prob},"
            f"b={args.b},g=all.stim"
        )
    else:
        raise ValueError(f"Unsupported noise model: {args.noise_model}")

    # Write the circuit to file
    with open(filename, "w") as f:
        f.write(str(circuit))

    print(f"✓ Generated circuit: {filename}")

if __name__ == "__main__":
    main()