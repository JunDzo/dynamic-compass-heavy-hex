import matplotlib.pyplot as plt
import matplotlib.patches as patches

# === USER INPUT ===
d = int(input("Enter grid size d: "))
steps = int(input("Enter number of steps: "))
output_file = input("Output file name (e.g., table.txt, press Enter for default): ")
if not output_file.strip():
    output_file = "table.txt"

# === Function to run one selection ===
def run_selection(step):
    selected = set()
    fig, ax = plt.subplots()
    ax.set_xlim(0, d-1)
    ax.set_ylim(0, d-1)
    ax.set_aspect('equal')
    ax.set_xticks(range(d))
    ax.set_yticks(range(d))
    ax.grid(True)
    ax.invert_yaxis()  # (0,0) at top-left
    plt.title(f"Step {step+1}: Click squares to toggle, only (row+col) odd allowed.")

    rects = {}
    for row in range(d-1):
        for col in range(d-1):
            rect = patches.Rectangle((col, row), 1, 1, facecolor='white', edgecolor='black')
            ax.add_patch(rect)
            rects[(row, col)] = rect
            # Optionally shade invalid squares:
            if (row + col) % 2 == 0:
                rect.set_facecolor('#cccccc')  # gray out even squares

    def onclick(event):
        if event.inaxes != ax:
            return
        col = int(event.xdata)
        row = int(event.ydata)
        if 0 <= col < d and 0 <= row < d:
            if (row + col) % 2 == 0:
                print(f"Square ({row},{col}) not allowed (row+col even).")
                return
            key = (row, col)
            if key in selected:
                selected.remove(key)
                rects[key].set_facecolor('white')
            else:
                selected.add(key)
                rects[key].set_facecolor('red')
        fig.canvas.draw()

    fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()

    return [complex(row, col) for (row, col) in sorted(selected)]

# === Run all steps ===
all_steps = []
for step in range(steps):
    result = run_selection(step)
    all_steps.append(result)

# === Save ===
with open(output_file, "w") as f:
    f.write(str(all_steps))

print(f"Saved {steps} selection(s) to {output_file}")