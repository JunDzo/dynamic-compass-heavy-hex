import _heavy_hex as hh
import gen


import os

os.makedirs("out/circuits", exist_ok=True)


for d in range(3, 11):
    table = [[], []]
    all_coords = hh.get_dataset(d)
    r = d * 4
    p = 0.001
    noise = "uniform"
    c_name = "heavy_hex"
    q = len(all_coords)
    b = "Z"
    g = "all.stim"
    
    # print(sorted(all_coords, key=lambda a: (a.real, a.imag)))

    # builder = gen.Builder.for_qubits(all_coords)

    # recorded_measurements = set()
    # fourbody = hh.FourBodyTileGenerator(d)
    # for step in range(r):
    #     is_end = True if step == r - 1 else False
    #     twobody = hh.TwoBodyTileGenerator(table=table[step % 2], code_distance=d)
    #     circuitbuilder = hh.HeavyHexCircuitBuilder(
    #         builder,
    #         d,
    #         step,
    #         fourbody.generate_4body_check(),
    #         twobody.generate_2body_checks(),
    #         recorded_measurements
    #     )
    #     circuitbuilder.generate_round_circuit(b,is_end)
    #     detector = hh.ConstructDetectors(
    #         builder,
    #         d,
    #         step,
    #         twobody.generate_2body_checks(),
    #         fourbody.generate_4body_check(),
    #         recorded_measurements,
    #     )
    #     if step == 0:
    #         a = detector.init_detectors(b)
    #         # print(a)
    #     else:
    #         detector.detector_generator(table[(step - 1) % 2],b, is_end)

    # circuit = builder.circuit
    circuit = hh.make_heavy_hex_circuit(table=table, diameter=d, rounds=r, basis=b)
    circuit = gen.NoiseModel.uniform_depolarizing(p).noisy_circuit(circuit)

    filename = f"out/circuits/r={r},d={d},p={p},noise={noise},c={c_name},q={q},b={b},g={g}"

    with open(filename, "w") as f:
        f.write(str(circuit))

