#!/bin/sh

# Checking for Operating System
if [ "$(uname)" = "Darwin" ]; then
    # Install dependencies and Generate Makefile on Mac
    ./scripts/build_mac.sh
elif [ "$(expr substr $(uname -s) 1 5)" = "Linux" ]; then
    # Install dependencies and Generate Makefile on Linux
    ./scripts/build_linux.sh
else
    echo "Warning: Your environment is currently not supported by autobuild."
    ./configure && make
fi
