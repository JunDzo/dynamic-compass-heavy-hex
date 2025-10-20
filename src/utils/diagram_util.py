import stim
import os
import cairosvg

def load_circuit(path: str) -> stim.Circuit:
    """Load a circuit from a .stim file."""
    return stim.Circuit.from_file(path)

def export_svg(
    circuit: stim.Circuit,
    diagram_type: str,
    svg_path: str,
    noise: bool = False,
    tick: int | range | None = None,
):
    """Generate a diagram SVG and save it."""
    os.makedirs(os.path.dirname(svg_path), exist_ok=True)
    if noise:
        diagram = circuit.diagram(diagram_type, tick=tick)
    else:
        diagram = circuit.without_noise().diagram(diagram_type, tick=tick)
    with open(svg_path, "w") as svg_file:
        svg_file.write(str(diagram))


def export_pdf(
    circuit: stim.Circuit,
    diagram_type: str,
    pdf_path: str,
    noise: bool = False,
    tick: int | range | None = None,
):
    """Generate a diagram PDF (via temporary SVG).
    eg: export_pdf(
    circuit = load_circuit("path/to/circuit.stim"),
    diagram_type='timeline-svg',
    pdf_path='viewer/timeline.pdf',
    tick=range(0,27)
)
    """
    os.makedirs(os.path.dirname(pdf_path), exist_ok=True)

    # Temporary SVG (same basename as PDF, just with .svg)
    svg_temp = os.path.splitext(pdf_path)[0] + ".svg"

    export_svg(circuit, diagram_type, svg_temp, noise=noise, tick=tick)
    cairosvg.svg2pdf(url=svg_temp, write_to=pdf_path)

def filter_detectors_and_measurements(circuit: stim.Circuit, path: str) -> stim.Circuit:
    """Return the circuit text with only DETECTOR and M instructions."""
    lines = str(circuit).splitlines()
    filtered = [line for line in lines if line.startswith("DETECTOR") or line.startswith("M") or line.startswith("MX") or line.startswith("QUBIT_COORDS") ]  #or line.startswith("TICK")
    with open(path, "w") as f:
        f.write("\n".join(filtered))
    return load_circuit(path)

def truncate_circuit_by_ticks(circuit: stim.Circuit, max_ticks: int) -> stim.Circuit:
    """Return a truncated circuit containing only up to `max_ticks` TICKs."""
    result = []
    tick_count = 0
    for line in str(circuit).splitlines():
        result.append(line)
        if line.strip() == "TICK":
            tick_count += 1
            if tick_count >= max_ticks:
                break
    return stim.Circuit("\n".join(result))

path = "out-test/circuits/r=16,d=4,p=0.001,noise=uniform,c=heavy_hex,b=X,g=all.stim"
circuit = load_circuit(path)
path = "viewer/filtered_circuit.stim"
# circuit = filter_detectors_and_measurements(circuit, path)

export_pdf(
    circuit,
    diagram_type='timeline-svg',
    pdf_path='viewer/timeline1.pdf',
    noise=True,
    tick=range(0,27)
)


