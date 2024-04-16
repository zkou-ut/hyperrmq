import sys

import numpy as np
import matplotlib.pyplot as plt

for_paper = False


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


huffman_lines = lines[: 10 * 18]

huffman_step = 18


def lines_to_np(lines, start, step, index=-1):
    return np.array(
        list(map(lambda line: int(line.split()[index]), lines[start::step]))
    )


# huffman
Bs = lines_to_np(huffman_lines, 1, huffman_step)
total = lines_to_np(huffman_lines, 3, huffman_step)
rmm = lines_to_np(huffman_lines, 4, huffman_step)
main = lines_to_np(huffman_lines, 9, huffman_step) + lines_to_np(
    huffman_lines, 10, huffman_step
)
table = lines_to_np(huffman_lines, 11, huffman_step)
sample = lines_to_np(huffman_lines, 12, huffman_step)

querytime = lines_to_np(huffman_lines, 17, huffman_step, -2)

assert np.all(total == rmm + main + table + sample)

stack_plot = {
    "Sequence of Huffman codes": main / n,
    "BP of top-tier tree": rmm / n,
    "Decoding table": table / n,
    "Sparse sampling": sample / n,
}

fig, axes = plt.subplots(ncols=2, figsize=(12, 7))
bottom = np.zeros(Bs.shape)

colors = ["C0", "C1", "C2", "C3"]
assert len(colors) == len(stack_plot)

ax = axes[0]
ax.grid(axis="y")
for color, (label, weights) in zip(colors, stack_plot.items()):
    p = ax.bar(Bs, weights, 0.8, label=label, bottom=bottom, color=color)
    bottom += weights

ax.legend(bbox_to_anchor=(0, -0.25), loc="upper left", borderaxespad=0, fontsize=18)
ax.set_xlim((min(Bs) - 0.6, max(Bs) + 0.6))
ax.set_ylim(0, 4)
ax.set_xticks(Bs)
ax.set_xlabel("Tree-covering parameter $B$")
ax.set_ylabel("Space consumption (bpe)")

ax = axes[1]
ax.plot(Bs, querytime / q * 1e3, marker="o")
ax.grid()
ax.set_xticks(Bs)
ax.set_ylim(0, 20)
ax.set_xlabel("Tree-covering parameter $B$")
ax.set_ylabel("Average query time [$\\mu$s]")

plt.tight_layout()
plt.savefig(fig_name)
