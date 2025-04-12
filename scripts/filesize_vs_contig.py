import matplotlib.pyplot as plt

# Data mapping filenames to file sizes and contigs
file_data = {
    "human-chr14-synthetic.txt": (4934090810, 860329),
    "test.txt": (103826531, 5736),
    "small.txt": (22473576, 1000),
    "little.txt": (3354849, 100),
    "verysmall.txt": (470672, 10),
    "tiny.txt": (82340, 1),
}

# Extract file sizes and contigs
file_sizes = [data[0] for data in file_data.values()]
contigs = [data[1] for data in file_data.values()]
file_names = list(file_data.keys())

# Create log-log plot
plt.figure(figsize=(8, 6))
plt.loglog(file_sizes, contigs, marker='o', linestyle='-', label="File Data")

# Annotate each point with the file name
for i, filename in enumerate(file_names):
    plt.annotate(filename, (file_sizes[i], contigs[i]), textcoords="offset points", xytext=(5,5), ha='right')

# Labels, legend, and grid
plt.xlabel('File Size (bytes)')
plt.ylabel('Number of Contigs')
plt.title('File Size vs. Contigs')
plt.grid(True, which="both", linestyle="--")

# Save plot as PNG
plt.savefig("file_size_vs_contigs.png", dpi=300)
plt.show()
