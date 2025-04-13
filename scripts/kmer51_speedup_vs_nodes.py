import matplotlib.pyplot as plt
import numpy as np

# Serial baseline time
serial_time = 42.201469

# Parallel assembled times for limited dataset (60 and 64 tasks per node)
parallel_assembly_times = {
    60: {
        1: 13.883080,
        2: 9.946318,
        3: 6.978182,
        4: 4.890504,
    },
    64: {
        1: 6.791644,
        2: 5.598765,
        3: 4.315666,
        4: 3.783803,
    },
}

node_counts = np.array([1, 2, 3, 4])
color = "blue"  # Common color for both lines

# Compute speedup for each dataset
speedup_data = {}
for tasks, data in parallel_assembly_times.items():
    speedup_data[tasks] = {N: serial_time / data[N] for N in data}

# **Plot Speedup for 60 and 64 Tasks**
plt.figure(figsize=(8, 6))
for tasks, speedups in speedup_data.items():
    linestyle = ':' if tasks == 64 else '-'  # Dotted line for 64 tasks
    plt.plot(node_counts, list(speedups.values()), marker='o', linestyle=linestyle, color=color)

# Annotate second-to-last data point
for tasks, speedups in speedup_data.items():
    plt.annotate(
        f"human-chr14-synthetic.txt (Tasks={tasks})",
        (node_counts[1], speedups[node_counts[1]]),
        textcoords="offset points",
        xytext=(5,-10),
        ha='left',
        fontsize=10,
        color=color
    )

# Labels, grid, and title
plt.xlabel("Node Count (N)")
plt.ylabel("Speedup (Tref / Talg) [log scale]")
plt.xticks(node_counts)  # Ensure only whole numbers are labeled
plt.yscale('log')  # Log scale for Y-axis
plt.title("Speedup vs. Node Count (parallel kmer51)")
plt.grid(True, which="both", linestyle="--")

# Save plot as PNG
plt.savefig("kmer51_speedup_vs_nodes.png", dpi=300)
plt.show()
