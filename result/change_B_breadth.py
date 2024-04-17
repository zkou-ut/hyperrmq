import sys

import matplotlib.pyplot as plt
import numpy as np

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
print(f"{n = }, {q = }")


step = 19


def lines_to_np(lines, start, step, index=-1):
    return np.array(
        list(map(lambda line: int(line.split()[index]), lines[start::step]))
    )


lines = lines[:]

# huffman
Bs = lines_to_np(lines, 1, step)
total = lines_to_np(lines, 3, step)
rmm = lines_to_np(lines, 4, step)
main = lines_to_np(lines, 9, step) + lines_to_np(lines, 10, step)
split_rank = lines_to_np(lines, 11, step)
node_count = lines_to_np(lines, 12, step)
sample = lines_to_np(lines, 13, step)
querytime = lines_to_np(lines, 18, step, -2)

assert np.all(total == rmm + main + split_rank + node_count + sample)

stack_plot = {
    "Sequence of arithmetic codes": main / n,
    "BP of top-tier tree": rmm / n,
    "Micro-tree sizes": node_count / n,
    "Split ranks": split_rank / n,
    "Sparse sampling": sample / n,
}

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 5.5))
bottom = np.zeros(Bs.shape)

xs = np.arange(len(Bs)) + 1

colors = ["C0", "C1", "C4", "C6", "C3"]
assert len(colors) == len(stack_plot)

ax = axes[0]
ax.grid(axis="y")
for color, (label, weights) in zip(colors, stack_plot.items()):
    p = ax.bar(xs, weights, 0.8, label=label, bottom=bottom, color=color)
    bottom += weights

xticks = list(range(5, 21, 5))

ax.legend(fontsize=18)
# ax.legend(bbox_to_anchor=(0, -0.2), loc="upper left", borderaxespad=0)
ax.set_xlim((min(xs) - 0.6, max(xs) + 0.6))
ax.set_yticks(np.arange(0, 7))
ax.set_xticks(xticks)
ax.set_xticklabels(["$2^{" + str(x) + "}$" for x in xticks])
# ax.set_xticklabels(Bs, rotation=60)
ax.set_xlabel("Tree-covering parameter $B$")
ax.set_ylabel("Space consumption (bpe)")

ax = axes[1]
ax.plot(xs, querytime / q * 1e3, marker="o")
ax.grid()
ax.set_ylim(0, 10)
ax.set_xticks(xticks)
ax.set_xticklabels(["$2^{" + str(x) + "}$" for x in xticks])
ax.set_xlabel("Tree-covering parameter $B$")
ax.set_ylabel("Average query time [$\\mu$s]")

plt.tight_layout()
plt.savefig("change_B_arith_plots.pdf")

plt.tight_layout()
plt.savefig(fig_name)
