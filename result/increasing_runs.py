import sys
from functools import lru_cache
from math import log2

import matplotlib.pyplot as plt
import numpy as np

N = 10**8
expected_rs = [N // 1000, N // 100, N // 10]

for_paper = True


plt.rcParams["font.size"] = 16

if for_paper:
    plt.rcParams["text.usetex"] = True
    plt.rcParams["font.family"] = "Times New Roman"

fig_ext = ".pdf" if for_paper else ".png"

name_no_ext = sys.argv[0][:-3]
file_names = [f"{name_no_ext}_{N}_{r}.txt" for r in expected_rs]

item_line_count = 16
label_to_line_indices = {
    "Sequence of Huffman codes": [7, 8],
    "BP of top-tier tree": [2],
    "Decoding table": [9],
    "Sparse sampling": [10],
}

Bs = list(range(5, 101, 5))


def log_combination(n, k):
    ret = 0
    for i in range(1, k + 1):
        ret += log2(n - i + 1)
        ret -= log2(i)
    return ret


@lru_cache(None)
def log_narayana(n, k):
    return 2 * log_combination(n, k)  # + log2(k) - log2(n - k + 1) - log2(n)


rs = []
data_map = dict()
time_map = dict()
for expected_r, file_name in zip(expected_rs, file_names):

    with open(file_name) as f:
        lines = f.readlines()

    idx = 2

    def check(name, value):
        global idx
        assert (
            lines[idx].rstrip() == f"{name} = {value}"
        ), f"{lines[idx].rstrip()} / {name} = {value}"
        idx += 1

    check("N", N)
    check("Q", 10**6)

    check("r(expected)", expected_r)
    r = int(lines[idx].split()[-1])
    rs.append(r)
    check("r", r)

    for B in Bs:
        idx += 1
        check("B", B)
        block = lines[idx : idx + item_line_count]
        idx += item_line_count
        total_actual = 0
        for label, line_indices in label_to_line_indices.items():
            val = 0
            for line_index in line_indices:
                val += int(block[line_index].split()[-1])
            data_map[(label, r, B)] = val
            total_actual += val
        total_expected = int(block[1].split()[-1])
        assert (
            total_actual == total_expected
        ), f"{total_actual} {total_expected} \n {block}"
        # print(block[-1])
        time_map[(r, B)] = int(block[-1].split()[-2])

# print(time_map)

fig, axes = plt.subplots(2, 3, figsize=(10, 5), gridspec_kw={"height_ratios": [4, 1]})
fig.subplots_adjust(wspace=0.4)
xs = np.arange(len(Bs))
titles = [f"$r = {r}" for r in rs]
titles[0] += r"= 10^5$"
titles[1] += r"\sim 10^6$"
titles[2] += r"\sim 10^7$"
for r, ax, title in zip(rs, axes[0], titles):
    stack_plot = {
        key: [data_map[(key, r, B)] / N for B in Bs] for key in label_to_line_indices
    }

    # ax.stackplot(Bs, stack_plot.values(), labels=stack_plot.keys(), alpha=0.7)
    ax.grid(axis="x")
    bottom = np.zeros(len(Bs))
    for label, weights in stack_plot.items():
        p = ax.bar(xs, weights, 0.8, label=label, bottom=bottom)
        bottom += weights

    ax.axhline(
        log_narayana(N, r) / N,
        0,
        1,
        color="black",
        linestyle="dotted",
        linewidth=2,
    )

    ax.grid()
    ax.set_xlabel("$B$")
    ax.set_ylim(0, log_narayana(N, r) / N * 8)
    ax.set_xlim(xs[0] - 0.6, xs[-1] + 0.6)
    ax.set_xticks(xs[3::4])
    ax.set_xticklabels(Bs[3::4])
    ax.set_title(title, fontsize=16)

for ax in axes[1]:
    fig.delaxes(ax)

axes[0, 0].set_ylabel("Space consumption (bpe)")
axes[0, 0].legend(bbox_to_anchor=(0, -0.25), loc="upper left", borderaxespad=0, ncol=2)
# plt.tight_layout()
plt.savefig(name_no_ext + fig_ext)
