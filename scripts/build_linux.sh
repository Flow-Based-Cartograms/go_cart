
# This echoes the commands to STDOUT, so the user can see what's going on.
set -x
sudo apt-get update
sudo apt-get install libfftw3-3 libfftw3-dev && \
sudo apt-get install build-essential pkg-config autoconf automake

set +x
./scripts/build_cjson.sh || { echo "Installing dependency cJSON failed."; exit 1; }
set -x

./autogen.sh && ./configure && make clean && make && cp cartogram_generator/cartogram .
# Turn off echoing commands
set +x
