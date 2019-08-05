# This echoes the commands to STDOUT, so the user can see what's going on.
set -x
brew install cjson fftw pkg-config autoconf automake && \
gcc_version=$(brew install gcc 2>&1 | tee /dev/tty | egrep -o "(gcc)(\D*)(\d+)" | head -1 | sed -E 's|^[^0-9]*||') && \
./autogen.sh && CC=gcc-$gcc_version ./configure && make clean && make
set +x
# Turn off echoing commands