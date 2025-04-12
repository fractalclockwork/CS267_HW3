import matplotlib.pyplot as plt
import numpy as np

# Serial baseline time
serial_time = 42.201469

# Parallel assembled times from results_parallel_kmer51.txt
data = {
    1: [136.002856, 138.320711, 108.482213, 60.895227, 35.947635, 19.719521, 10.485733, 7.727779],
    2: [216.270044, 141.978644, 78.036363, 38.999279, 22.620205, 12.206918, 6.791644, 3.783803],
    3: [261.829604, 138.418709, 70.848890, 35.702402, 18.854000, 10.323634, 5.598765, 3.053621],
    4: [215.369635, 112.417470, 55.271559, 28.489452, 15.041629, 8.094177, 4.315666, 2.461360]
}

ntasks_per_node = np.array([1, 2, 4, 8, 16, 32, 64, 128])

# Compute speedup for each N
speedup_data = {N: [serial_time / t for t in times] for N, times in data.items()}

# Find minimum initial speedup value (first entry of each N)
min_initial_speedup = min(speedup_data[N][0] for N in speedup_data)

# Adjust O(n) reference line to start at min_initial_speedup
ideal_speedup = (ntasks_per_node / ntasks_per_node[0]) * min_initial_speedup

# Create log-log speedup plot
plt.figure(figsize=(8, 6))
for N, speedup in speedup_data.items():
    plt.loglog(ntasks_per_node, speedup, marker='o', linestyle='-', label=f'N={N}')

# Add normalized reference line for O(n) starting at min initial speedup
plt.loglog(ntasks_per_node, ideal_speedup, '--', color='gray', label='O(n) Ideal Speedup')

# Labels, legend, and grid
plt.xlabel('Number of tasks per node')
plt.ylabel('Speedup')
plt.title('Speedup Plot (parallel kmer51 human-chromosome)')
plt.legend()
plt.grid(True, which="both", linestyle="--")

# Save plot as PNG
plt.savefig("speedup_plot.png", dpi=300)
plt.show()
