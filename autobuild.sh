# Checking for Operating System
if [ "$(uname)" == "Darwin" ]; then
    # Install dependencies and Generate Makefile on Mac
    bash scripts/build_mac.sh
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Install dependencies and Generate Makefile on Linux
    bash scripts/build_linux.sh
else
    echo "Warning: Your environment is currently not supported by autobuild."
    ./configure && make
fi
