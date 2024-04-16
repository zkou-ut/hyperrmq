set -e

# build
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make

# test
# ./simple-test
# ./test/hyperrmq-test
