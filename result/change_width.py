import sys
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

for_paper = True


name_no_ext = sys.argv[0][:-3]
txt_name = name_no_ext + ".txt"
fig_ext = ".pdf" if for_paper else ".png"

plt.rcParams["font.size"] = 22
if for_paper:
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "Times New Roman"

with open(txt_name) as f:
    lines = f.readlines()


for _ in range(2):
    print(lines.pop(0).rstrip())

N = int(lines.pop(0).split()[-1])
Q = int(lines.pop(0).split()[-1])

print(f"{N = }")
print(f"{Q = }")

ws = [2]
while ws[-1] * 2 <= N:
    ws.append(ws[-1] * 2)

step = len(ws) + 2

labels = ["RMM", "REC", "H"]
s = 64
while s <= 1024:
    labels.append("B" + str(s))
    s *= 2
labels.append("D512")
labels.append("BNP512")
labels.append("DNP512")

space_list = []
times_list = []
for index, label in enumerate(labels):
    start = index * step
    space_list.append(int(lines[start + 1]))
    times_list.append(
        list(map(lambda s: float(s.split()[-1]), lines[start + 2 : start + step]))
    )

# comparison

used_indices = list(range(8))
markers = ["o", "D", "p", "x", "v", "^", "<", ">"]

fig, ax = plt.subplots(figsize=(8, 5))

for index, marker in zip(used_indices, markers):
    label = labels[index]
    times = times_list[index]
    ax.plot(
        np.arange(1, len(ws) + 1),
        times,
        marker=marker,
        label=label,
    )

ax.set_xticks(np.arange(0, 31, 5))
ax.set_xticklabels(["$2^{" + str(x) + "}$" for x in range(0, 31, 5)])
ax.grid()
ax.set_xlabel("Query width")
ax.set_ylabel("Average query time [$\\mu$s]")
ax.set_ylim((0, 23))
ax.legend(bbox_to_anchor=(1.1, 1), loc="upper left", borderaxespad=0)

plt.tight_layout()
plt.savefig(name_no_ext + "_comparison" + fig_ext)


# depth breadth

fig, ax = plt.subplots(figsize=(8, 5))

depth_time = times_list[8]
breadth_time = times_list[6]
depth_no_pruning_time = times_list[10]
breadth_no_pruning_time = times_list[9]

ax.plot(
    np.arange(1, len(ws) + 1),
    depth_time,
    marker="o",
    color="C0",
    label="Depth-first with pruning",
)
ax.plot(
    np.arange(1, len(ws) + 1),
    breadth_time,
    marker="o",
    color="C1",
    label="Breadth-first with pruning",
)
ax.plot(
    np.arange(1, len(ws) + 1),
    depth_no_pruning_time,
    marker="o",
    color="C0",
    label="Depth-first without pruning",
    linestyle=":",
)
ax.plot(
    np.arange(1, len(ws) + 1),
    breadth_no_pruning_time,
    marker="o",
    color="C1",
    label="Breadth-first without pruning",
    linestyle=":",
)

ax.set_xticks(np.arange(0, 31, 5))
ax.set_xticklabels(["$2^{" + str(x) + "}$" for x in range(0, 31, 5)])
ax.grid()
ax.set_xlabel("Query width")
ax.set_ylabel("Average query time [$\\mu$s]")
ax.set_ylim(0, 40)
ax.legend(fontsize=18)

plt.tight_layout()
plt.savefig(name_no_ext + "_depth_breadth" + fig_ext)
