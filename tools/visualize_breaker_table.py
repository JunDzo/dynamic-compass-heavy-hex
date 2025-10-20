#!/usr/bin/env python3
"""
Visualize Breaker Table

Renders saved breaker configurations as PDF plots for verification and documentation.

Usage:
    python tools/visualize_breaker_table.py

Configuration:
    Edit the 'd' and 'table_file' variables in the main section.
"""
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

# === Visualization Function ===
d = int(input("Enter grid size d: "))

def visualize_table(file_path, output_dir="data/plots/breakers"):
    """
    Visualize breaker table as grid plots.

    Args:
        file_path: Path to breaker table file (Python literal)
        output_dir: Directory to save visualization PDFs

    Output:
        visualized_step_<N>.pdf for each step in the table
    """
    import ast

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    with open(file_path, "r") as f:
        data = ast.literal_eval(f.read())

    for step, step_data in enumerate(data):
        fig, ax = plt.subplots()
        ax.set_xlim(0, d-1)
        ax.set_ylim(0, d-1)
        ax.set_aspect('equal')
        ax.set_xticks(range(d))
        ax.set_yticks(range(d))
        ax.grid(True)
        ax.invert_yaxis()
        plt.title(f"Visualize Step {step}")

        for row in range(d-1):
            for col in range(d-1):
                facecolor = 'white'
                if (row + col) % 2 == 0:
                    facecolor = '#cccccc'
                rect = patches.Rectangle((col, row), 1, 1, facecolor=facecolor, edgecolor='black')
                ax.add_patch(rect)

        for pos in step_data:
            row, col = int(pos.real), int(pos.imag)
            rect = patches.Rectangle((col, row), 1, 1, facecolor='red', edgecolor='black')
            ax.add_patch(rect)

        output_file = os.path.join(output_dir, f"visualized_step_{step}.pdf")
        plt.savefig(output_file, bbox_inches='tight')
        plt.close()
        print(f"✓ Saved: {output_file}")

# === Main ===
if __name__ == "__main__":
    table_file = input("Enter breaker table file path (default: data/analysis/breaker_table.txt): ")
    if not table_file.strip():
        table_file = "data/analysis/breaker_table.txt"

    visualize_table(table_file)