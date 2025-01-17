import subprocess
import tempfile
from multiprocessing import Pool

N = 10**8
rs = [N // 1000, N // 100, N // 10]
parallel_execution = True


def execute(r):
    print(r, "start")
    tmp = tempfile.TemporaryFile()
    p = subprocess.Popen(
        ["./build/experiment/increasing_runs", str(N), str(r)], stdout=tmp
    )
    print(r, "finished", p.wait())
    out_name = f"result/increasing_runs_{N}_{r}.txt"
    tmp.seek(0)
    with open(out_name, "wb") as f:
        f.write(tmp.read())
    print(r, "written")


if parallel_execution:
    with Pool() as p:
        p.map(execute, rs)
else:
    for r in rs:
        execute(r)
