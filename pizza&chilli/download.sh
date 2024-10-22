# This file is retrieved from the following repository:
# https://github.com/kittobi1992/rmq-experiments/blob/master/pizza%26chilli/download.sh
# The baseurl is modified in this file.

baseurl='https://pizzachili.dcc.uchile.cl/texts'

files='
./code/sources.200MB.gz
./xml/dblp.xml.200MB.gz
./dna/dna.200MB.gz
./nlang/english.200MB.gz
'

mkdir build

for f in ${files}; do
    echo "${f}"
    dir=$(dirname "${f}")
    name=$(basename "${f}")

    echo "curl -O ${baseurl}/${f}"
    curl -O ${baseurl}/${f}
done

gunzip *


