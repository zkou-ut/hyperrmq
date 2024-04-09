import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

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

for _ in range(3):
    print(lines.pop(0).rstrip())

labels = ["F\&N", "NewRMQ", "H"]
s = 64
while s <= 1024:
    labels.append("B" + str(s))
    s *= 2

Ns = []
label_to_const = defaultdict(list)
label_to_query = defaultdict(list)
label_to_space = defaultdict(list)
line_index = 0
while line_index < len(lines):
    N = int(lines[line_index].split()[-1])
    line_index += 1
    Ns.append(N)
    for label in labels:
        line_index += 1
        label_to_const[label].append(float(lines[line_index].split()[-1]))
        line_index += 1
        label_to_query[label].append(float(lines[line_index].split()[-1]))
        line_index += 1
        label_to_space[label].append(float(lines[line_index].split()[-1]))
        line_index += 1


markers = [
    "o",
    "D",
    "p",
    "x",
    "v",
    "^",
    "<",
    ">",
]


fig, axes = plt.subplots(ncols=2, figsize=(13, 5))
ax = axes[0]
Ns = np.array(Ns)
midx = 0
for label in labels:
    space = np.array(label_to_space[label])
    ax.plot(Ns, space, label=label, marker=markers[midx])
    midx += 1
ax.grid()
ax.set_xscale("log")
ax.set_xticks(10 ** np.arange(4, 10))
ax.set_ylim(0, 5)
ax.set_xlabel("Array length $n$")
ax.set_ylabel("Space consumption (bpe)")


ax = axes[1]
Ns = np.array(Ns)
midx = 0
for label in labels:
    time = np.array(label_to_query[label])
    ax.plot(Ns, time / 1000, label=label, marker=markers[midx])
    midx += 1
ax.grid()
ax.set_xscale("log")
ax.set_xticks(10 ** np.arange(4, 10))
ax.set_xlabel("Array length $n$")
ax.set_ylabel("Average query time [$\\mu$s]")
ax.set_ylim(0, 20)
ax.legend(bbox_to_anchor=(1.1, 1), loc="upper left", borderaxespad=0)

plt.tight_layout()
plt.savefig(fig_name)
