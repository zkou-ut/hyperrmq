set -e

# build
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j

# run test
./test/hyperrmq-test