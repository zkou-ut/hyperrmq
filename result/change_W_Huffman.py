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


print(lines.pop(0).rstrip())  # date
lines.pop(0)
n = int(lines.pop(0).split()[-1])  # n
q = int(lines.pop(0).split()[-1])  # q


step = 19


def lines_to_np(lines, start, step, index=-1):
    return np.array(
        list(map(lambda line: int(line.split()[index]), lines[start::step]))
    )


Ws = lines_to_np(lines, 0, step)
total = lines_to_np(lines, 4, step)
query = lines_to_np(lines, 18, step, -2)

total = total / n
query = query / q * 1e3

necessary_indices = [6, 10, 12]
necessary_arr = 0
for idx in necessary_indices:
    necessary_arr += lines_to_np(lines, idx, step)

necessary = necessary_arr[0]
assert np.all(necessary == necessary_arr)
print(necessary / n)

fig, ax = plt.subplots()

ax.set_xlabel("Space consumption (bpe)")
ax.set_xlim(2, 5)
ax.set_yscale("log")
ax.set_ylabel("Average query time [$\mu$s]")
ax.set_ylim(10 / 12, 12 * 100)
ax.grid()
ax.plot(total, query, marker="o")

for i, (t, q, w) in enumerate(zip(total, query, Ws)):
    shiftx = 0.05
    prody = 1.1 ** (1 - i / 7)
    if w == 256:
        ax.text(t + shiftx, q * prody, str(w), ha="left", va="top")
    elif w >= 32:
        ax.text(t + shiftx, q * prody, str(w), ha="left", va="center")
    else:
        ax.text(t + shiftx, q * prody, str(w), ha="left", va="bottom")


ax.axvline(necessary / n, linestyle="dotted")
# ax.set_xlabel("Tree-covering parameter $B$")
# ax.set_ylabel("Space consumption (bpe)")

plt.tight_layout()
plt.savefig(fig_name)
