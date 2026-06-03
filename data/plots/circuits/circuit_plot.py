import stim
import os
import cairosvg

def load_circuit(path: str) -> stim.Circuit:
    """Load a circuit from a .stim file."""
    return stim.Circuit.from_file(path)

def dem(
    circuit: stim.Circuit,
    output_path: str | None = None,
):
    dem = circuit.detector_error_model()
    if output_path is not None:
        with open(output_path, "w") as f:
            f.write(str(dem))
    else:
        return dem

def export_svg(
    circuit: stim.Circuit,
    diagram_type: str,
    svg_path: str,
    noise: bool = False,
    tick: int | range | None = None,
    filter_coords: list | None = None,
):
    """Generate a diagram SVG and save it."""
    os.makedirs(os.path.dirname(svg_path), exist_ok=True)
    if noise:
        diagram = circuit.diagram(diagram_type, tick=tick, filter_coords=filter_coords)
    else:
        diagram = circuit.without_noise().diagram(diagram_type, tick=tick, filter_coords=filter_coords)
    with open(svg_path, "w") as svg_file:
        svg_file.write(str(diagram))

def export_pdf(
    circuit: stim.Circuit,
    diagram_type: str,
    pdf_path: str,
    noise: bool = False,
    tick: int | range | None = None,
    filter_coords: list | None = None,
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

    export_svg(circuit, diagram_type, svg_temp, noise=noise, tick=tick, filter_coords=filter_coords)
    cairosvg.svg2pdf(url=svg_temp, write_to=pdf_path)
    # Remove temporary SVG
    if os.path.exists(svg_temp):
        os.remove(svg_temp)

def make_output_path_from_circuit(
    circuit_path: str,
    diagram_type: str,
    out_dir: str,
    ext: str = ".pdf",
    tick: int | range | None = None,
    noise: bool = False,
) -> str:
    """
    Construct an output file path based on the circuit filename and diagram type.

    Example output:
        out_dir / "<circuit_stem>_<diagram_type>[_tick-0-25][_noise].pdf"
    """
    os.makedirs(out_dir, exist_ok=True)

    circuit_stem = os.path.splitext(os.path.basename(circuit_path))[0]
    diagram_tag = diagram_type.replace("-svg", "").replace("_svg", "")

    tick_tag = ""
    if isinstance(tick, range):
        tick_tag = f"_tick-{tick.start}-{tick.stop}"
    elif isinstance(tick, int):
        tick_tag = f"_tick-{tick}"

    noise_tag = "_noise" if noise else ""

    filename = f"{circuit_stem}_{diagram_tag}{tick_tag}{noise_tag}{ext}"
    return os.path.join(out_dir, filename)

def filter_detectors_and_measurements(circuit: stim.Circuit, path: str) -> stim.Circuit:
    """Return the circuit text with DETECTOR and QUBIT_COORDS lines removed."""
    lines = str(circuit).splitlines()
    # filtered = [line for line in lines if line.startswith("DETECTOR") or line.startswith("QUBIT_COORDS") or line.startswith("TICK") or line.startswith("M") or line.startswith("MX") ] #
    filtered = [
        line for line in lines
        if not line.lstrip().startswith("QUBIT_COORDS")
        and not line.lstrip().startswith("DETECTOR")
    ]
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

circuit_path = "figures/manuscript/no-reset-0/stats.csv"
circuit = load_circuit(circuit_path)
# filtered_circuit_path ="out-test/circuits/[filtered]r=2,d=2,p=0.001,noise=uniform,c=heavy_hex,b=X,g=all.stim"

path = "data/plots/circuits/plot"
# circuit = filter_detectors_and_measurements(circuit, filtered_circuit_path)

diagram_type = 'timeline-svg'
# tick_range = range(0, 5)
# tick_range=range(29, 39)
tick_range= None
# filter_coords=['D11', 'D17', ]
# range(25, 49)

auto_pdf_path = make_output_path_from_circuit(
    circuit_path=circuit_path,
    diagram_type=diagram_type,
    out_dir=path,
    ext='.pdf',
    tick=tick_range,
    noise=False,
)

export_pdf(
    circuit,
    diagram_type=diagram_type,
    pdf_path=auto_pdf_path,
    noise=False,
    tick=tick_range,
    # filter_coords=filter_coords
)

# dem(
#     circuit,
#     output_path='data/plots/circuits/r=8,d=5,noise=uniform,c=heavy_hex,b=Z,dem.txt'
# ) detslice-svg. detslice-with-ops-svg
# diagram_type='timeslice-svg',
# diagram_type='timeline-svg',
