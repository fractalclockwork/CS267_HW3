import matplotlib.pyplot as plt
import numpy as np

# Data from results_parallel_kmer51.txt
data = {
    1: [136.002856, 138.320711, 108.482213, 60.895227, 35.947635, 19.719521, 10.485733, 7.727779],
    2: [216.270044, 141.978644, 78.036363, 38.999279, 22.620205, 12.206918, 6.791644, 3.783803],
    3: [261.829604, 138.418709, 70.848890, 35.702402, 18.854000, 10.323634, 5.598765, 3.053621],
    4: [215.369635, 112.417470, 55.271559, 28.489452, 15.041629, 8.094177, 4.315666,  2.461360] 
}

ntasks_per_node = [1, 2, 4, 8, 16, 32, 64, 128]

# Create log-log plot
plt.figure(figsize=(8, 6))
for N, times in data.items():
    plt.loglog(ntasks_per_node, times, marker='o', linestyle='-', label=f'N={N}')

# Labels, legend, and grid
plt.xlabel('Number of tasks per node')
plt.ylabel('Assembled Time (s)')
plt.title('Raw Performance (parallel kmer51)')
plt.legend()
plt.grid(True, which="both", linestyle="--")

# Save plot as PNG
plt.savefig("kmer51_performance_plot.png", dpi=300)
#plt.show()

