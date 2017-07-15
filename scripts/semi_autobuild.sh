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
