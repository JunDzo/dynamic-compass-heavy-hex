import _heavy_hex as hh
import gen
import argparse
import os
import ast

parser = argparse.ArgumentParser()
parser.add_argument("--out_dir", type=str, default="out/circuits")
parser.add_argument("--diameter", type=int, required=True)
parser.add_argument("--rounds", type=int, required=True)
parser.add_argument("--noise_model", type=str, default="uniform")
parser.add_argument("--noise_strength", type=float, default=0.001)
parser.add_argument("--style", type=str, default="heavy_hex")
parser.add_argument("--b", type=str, default="Z")
parser.add_argument("--table_file", type=str, help="Path to a file containing the table variable as a Python literal")
args = parser.parse_args()


if args.table_file:
    with open(args.table_file) as f:
        table = ast.literal_eval(f.read())
else:
    table = [[], []]

os.makedirs(args.out_dir, exist_ok=True)

d = args.diameter
all_coords = hh.get_dataset(d)
r = args.rounds
p = args.noise_strength
noise = args.noise_model
c_name = args.style
q = len(all_coords)
b = args.b
g = "all.stim"

circuit = hh.make_heavy_hex_circuit(table=table, diameter=d, rounds=r, basis=b)
circuit = gen.NoiseModel.uniform_depolarizing(p).noisy_circuit(circuit)

filename = f"{args.out_dir}/r={r},d={d},p={p},noise={noise},c={c_name},q={q},b={b},g={g}"

with open(filename, "w") as f:
    f.write(str(circuit))



