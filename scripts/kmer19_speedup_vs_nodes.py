import matplotlib.pyplot as plt
import numpy as np

# Serial baseline assembly times
serial_assembly_times = {
    "test.txt": 4.729069,
    "small.txt": 0.828955,
    "little.txt": 0.080258,
    "verysmall.txt": 0.011761,
    "tiny.txt": 0.001908,
}

# Parallel assembly times for limited dataset (60 and 64 tasks per node)
parallel_assembly_times = {
    60: {
        "test.txt": [0.528498, 0.484606, 0.432433, 0.329213],
        "small.txt": [0.127048, 0.133171, 0.130418, 0.119585],
        "little.txt": [0.031970, 0.033045, 0.041487, 0.044262],
        "verysmall.txt": [0.018429, 0.032551, 0.038666, 0.035282],
        "tiny.txt": [0.012586, 0.022592, 0.027050, 0.028342],
    },
    64: {
        "test.txt": [0.475044, 0.365957, 0.300998, 0.249688],
        "small.txt": [0.110340, 0.103425, 0.108469, 0.086336],
        "little.txt": [0.033846, 0.033826, 0.033654, 0.033911],
        "verysmall.txt": [0.017927, 0.026297, 0.028532, 0.027775],
        "tiny.txt": [0.013288, 0.017380, 0.019777, 0.020475],
    },
}

node_counts = np.array([1, 2, 3, 4])

# Compute speedup for each dataset
speedup_data = {}
for tasks, data in parallel_assembly_times.items():
    speedup_data[tasks] = {}
    for file_name, parallel_times in data.items():
        speedup_data[tasks][file_name] = [serial_assembly_times[file_name] / t for t in parallel_times]

# Sort entries by final speedup value (highest speedup first)
sorted_entries = []
for tasks, file_speeds in speedup_data.items():
    for file_name, speedups in file_speeds.items():
        sorted_entries.append((file_name, speedups, tasks))

sorted_entries = sorted(sorted_entries, key=lambda x: x[1][-1], reverse=True)

# Define consistent colors per file name
file_colors = {
    "test.txt": "blue",
    "small.txt": "green",
    "little.txt": "red",
    "verysmall.txt": "purple",
    "tiny.txt": "orange",
}

# **Plot Speedup for 60 and 64 Tasks**
plt.figure(figsize=(8, 6))
for file_name, speedups, tasks in sorted_entries:
    linestyle = ':' if tasks == 64 else '-'  # Dotted line for 64 tasks
    plt.plot(node_counts, speedups, marker='o', linestyle=linestyle, color=file_colors[file_name])
    
    # Annotate point at node = 3
    plt.annotate(f"{file_name} (Tasks={tasks})", (node_counts[-2], speedups[-2]), textcoords="offset points", xytext=(5,-10), ha='left', fontsize=10, color=file_colors[file_name])

# Labels, grid, and title
plt.xlabel("Node Count (N)")
plt.ylabel("Speedup (Tref / Talg) [log scale]")
plt.xticks(node_counts)  # Ensure only whole numbers are labeled on the x-axis
plt.yscale('log')  # Log scale for Y-axis
plt.title("Speedup vs. Node Count (parallel kmer19)")
plt.grid(True, which="both", linestyle="--")

# Save plot as PNG
plt.savefig("kmer19_speedup_vs_nodes.png", dpi=300)
plt.show()
