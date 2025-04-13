import matplotlib.pyplot as plt
import numpy as np

# Parallel assembly times
parallel_assembly_times = {
    1: {
        "test.txt": [5.042901, 6.043239, 4.510081, 2.543830, 1.537401, 0.897974, 0.475044, 0.292781],
        "small.txt": [0.867483, 1.007858, 0.839256, 0.521688, 0.359272, 0.204488, 0.110340, 0.086642],
        "little.txt": [0.084174, 0.141291, 0.133076, 0.100398, 0.074768, 0.043707, 0.033846, 0.029188],
        "verysmall.txt": [0.012007, 0.022241, 0.027699, 0.032243, 0.026502, 0.016592, 0.017927, 0.022166],
        "tiny.txt": [0.002076, 0.005519, 0.007461, 0.008355, 0.009506, 0.011215, 0.013288, 0.016019],
    },
}

ntasks_per_node = np.array([1, 2, 4, 8, 16, 32, 64, 128])

# **Plot Raw Performance**
plt.figure(figsize=(8, 6))
for file_name, times in parallel_assembly_times[1].items():
    plt.loglog(ntasks_per_node, times, marker='o', linestyle='-', label=file_name)

plt.xlabel("Number of Tasks per Node")
plt.ylabel("Assembly Time (s)")
plt.title("Raw Performance vs. Number of Tasks (parallel kmer19)")
plt.legend()
plt.grid(True, which="both", linestyle="--")
plt.savefig("kmer19_performance_plot.png", dpi=300)
plt.show()
