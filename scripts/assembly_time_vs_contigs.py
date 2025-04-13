import matplotlib.pyplot as plt

# Data mapping filenames to contigs and assembly times
file_data = {
    "human-chr14-synthetic.txt": (860329, 42.117088),
    "test.txt": (5736, 1.648523),
    "small.txt": (1000, 0.326189),
    "little.txt": (100, 0.042616),
    "verysmall.txt": (10, 0.005782),
    "tiny.txt": (1, 0.001006),
}

# Extract contigs and times
contigs = [data[0] for data in file_data.values()]
assembly_times = [data[1] for data in file_data.values()]
file_names = list(file_data.keys())

# Create log-log plot
plt.figure(figsize=(8, 6))
plt.loglog(contigs, assembly_times, marker='o', linestyle='-', label="Assembly Time")

# Annotate each point with the file name
for i, filename in enumerate(file_names):
    plt.annotate(filename, (contigs[i], assembly_times[i]), textcoords="offset points", xytext=(5,5), ha='right')

# Labels, legend, and grid
plt.xlabel('Number of Contigs')
plt.ylabel('Assembly Time (s)')
plt.title('Assembly Time vs. Contigs')
plt.grid(True, which="both", linestyle="--")

# Save plot as PNG
plt.savefig("assembly_time_vs_contigs.png", dpi=300)
plt.show()
