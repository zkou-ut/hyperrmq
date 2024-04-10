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


markers = ["o", "D", "p", "x", "v", "^", "<", ">", "+"]
labels = ["F\&N", "SDSLNEWRMQ", "H"]
s = 64
while s <= 1024:
    labels.append("B" + str(s))
    s *= 2
labels.append("D512")

space_list = []
times_list = []
for index, label in enumerate(labels):
    start = index * step
    space_list.append(int(lines[start + 1]))
    times_list.append(
        list(map(lambda s: float(s.split()[-1]), lines[start + 2 : start + step]))
    )


xs = np.arange(len(labels))

# fig, ax = plt.subplots()  # figsize=(6, 5))

# ax.bar(xs, space_list, 0.8)
# ax.set_xticks(xs)
# ax.set_xticklabels(labels, rotation=90)
# ax.grid()
# ax.set_xlabel("RMQs")
# ax.set_ylabel("Space consumption (bpe)")

# plt.tight_layout()
# plt.savefig("hist_space" + ext)

fig, ax = plt.subplots(figsize=(8, 5))

for label, times, marker in zip(labels, times_list, markers):
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
# ax.set_yscale("log")
ax.legend(bbox_to_anchor=(1.1, 1), loc="upper left", borderaxespad=0)

plt.tight_layout()
plt.savefig(fig_name)
