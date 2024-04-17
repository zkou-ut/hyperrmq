import subprocess
import tempfile
import glob
from multiprocessing import Pool

parallel_execution = False


files = glob.glob("pizza&chilli/*.200MB")

print(files)


def execute(file):
    print(file, "start")
    tmp = tempfile.TemporaryFile()
    p = subprocess.Popen(["./build/experiment/lcp", file, "pizza&chilli"], stdout=tmp)
    print(file, "finished", p.wait())
    out_name = f"result/lcp_{file.split('/')[-1]}.txt"
    tmp.seek(0)
    with open(out_name, "wb") as f:
        f.write(tmp.read())
    print(file, "written")


if parallel_execution:
    with Pool() as p:
        p.map(execute, files)
else:
    for file in files:
        execute(file)
