import _heavy_hex as hh
import gen
from collections import defaultdict


d = 2
rounds = 2
# Global dictionary to store all square tables by step


# table = [[0+1j, 1+2j, 2+1j, 3+2j],[0 + 3j , 1 + 0j, 2 + 3j , 3 + 0j]]
table=[[],[]]
all_coords=hh.get_dataset(d)
print(sorted(all_coords, key=lambda a: (a.real, a.imag)))

builder = gen.Builder.for_qubits(all_coords)
builder.tick()
recorded_measurements = set()
for step in range(rounds):
    twobody = hh.TwoBodyTileGenerator(table = table[step % 2], code_distance=d)
    fourbody = hh.FourBodyTileGenerator(d)
    circuitbuilder = hh.HeavyHexCircuitBuilder(
        builder,
        d,
        step,
        fourbody.generate_4body_check(),
        twobody.generate_2body_checks(),
        recorded_measurements
    )
    if step == 0:
        print(circuitbuilder._collect_qubit_sets())
        circuitbuilder.apply_circuit_stage("init")
    circuitbuilder.generate_round_circuit()
    detector = hh.ConstructDetectors(
        builder,
        d,
        step,
        twobody.generate_2body_checks(),
        fourbody.generate_4body_check(),
        recorded_measurements,
        table[step % 2]
    )
    if step == 0:
        detector.init_x_detectors()
    else:
        detector.detector_generator()
    
circuit = builder.circuit
circuit_test = circuit.diagram('timeline-svg')
if circuit_test is None:
    raise ValueError("diagram('timeline-svg') returned None — circuit may be empty or malformed.")


with open("circuit-timeline.svg", "w") as svg_file:
    svg_file.write(str(circuit_test))

# from gen._util import sorted_complex
# print(sorted_complex(all_coords))
