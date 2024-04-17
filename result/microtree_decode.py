import sys

import numpy as np
import matplotlib.pyplot as plt

for_paper = True


py_name = sys.argv[0]
txt_name = py_name.replace(".py", ".txt")
fig_name = py_name.replace(".py", ".pdf" if for_paper else ".png")

plt.rcParams["font.size"] = 22
if for_paper:
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "Times New Roman"

with open(txt_name) as f:
    lines = f.readlines()

for _ in range(2):
    print(lines.pop(0).rstrip())

ns = []
times = []
for line in lines:
    n, _, time = line.split()
    ns.append(int(n))
    times.append(float(time))

fig, ax = plt.subplots(figsize=(6, 5))

ax.plot(ns, times, marker="o")
ax.grid()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1 / 2, 2 * 1e6)
ax.set_xticks(2 ** np.arange(0, 21, 5))
ax.set_xticklabels("$2^{" + str(n) + "}$" for n in range(0, 21, 5))
ax.set_ylim(0.1, 1e5)
ax.set_yticks(10.0 ** np.arange(-1, 6))
ax.set_xlabel("Micro-tree size $n$")
ax.set_ylabel("Decoding time [$\\mu$s]")

plt.tight_layout()
plt.savefig(fig_name)
