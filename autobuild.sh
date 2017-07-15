# Checking for Operating System
if [ "$(uname)" == "Darwin" ]; then
    # Install dependencies and Generate Makefile on Mac
    bash scripts/prepare_env_mac.sh
elif [ "$(expr substr $(uname -s) 1 5)" == "Linux" ]; then
    # Install dependencies and Generate Makefile on Linux
    bash scripts/prepare_env_linux.sh
#elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW32_NT" ]; then
    # Do something under 32 bits Windows NT platform
#elif [ "$(expr substr $(uname -s) 1 10)" == "MINGW64_NT" ]; then
    # Do something under 64 bits Windows NT platform
fi

# Going into source code directory
cd cartogram_generator

# Building cartogram generator and placing in root directory
make clean &>/dev/null
echo "Cleaned up directory." 
make && \
cp cartogram ../ && \
echo "Successfully built cartogram generator." && \
echo "Run using: ./cartogram [.gen file] [.dat file]." || \
echo "An error occured. Please refer to troubleshooting instructions."

# Moving back to root directory
cd ..
