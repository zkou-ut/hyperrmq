import sys

from collections import defaultdict

folder = "./result/"

N = 200 * 2**20 + 1

labels = [
    "dblp.xml.200MB",
    "dna.200MB",
    "english.200MB",
    "sources.200MB",
]

method_time = defaultdict(list)
method_space = defaultdict(list)

for label in labels:
    result_file_name = folder + f"lcp_{label}.txt"

    with open(result_file_name) as f:
        lines = f.readlines()

    print("\n".join(lines[:6]))
    lines = lines[6:]
    for block_start in range(0, len(lines), 4):
        method = lines[block_start].rstrip()
        total_time = float(lines[block_start + 1].split()[-1])
        total_space = float(lines[block_start + 2].split()[-1])
        total_time = float(lines[block_start + 1].split()[-1])
        num_queries = float(lines[block_start + 3].split("=")[-1])
        method_time[method].append(total_time / num_queries)
        method_space[method].append(total_space / N)

print("Space table (for LaTeX)")
print("&", " & ".join(labels), r"\\")
for method, space in method_space.items():
    print(method, "&", " & ".join(map(lambda s: f"{s:.3f}", space)), r"\\")

print()
print("Time table (for LaTeX)")
print("&", " & ".join(labels), r"\\")
for method, time in method_time.items():
    print(method, "&", " & ".join(map(lambda t: str(round(t * 1e9)), time)), r"\\")
