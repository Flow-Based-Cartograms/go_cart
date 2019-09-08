#!/bin/sh

# Script to download, build, and install cJSON library.
# This is only necessary on Linux, because there is a Homebrew package


PKG_CONFIG=`which pkg-config`

echo -n "Checking whether cJSON is installed... "
if !  $PKG_CONFIG --exists libcjson; then

echo "no"
echo "Installing cJSON now..."

# This echoes the commands to STDOUT, so the user can see what's going on.
set -x

sudo apt-get install cmake git || { exit 1; }

TMP_DIR=`mktemp -d`

cd $TMP_DIR

git clone https://github.com/DaveGamble/cJSON || { exit 1; }

cd cJSON
mkdir build
cd build

cmake ..

make || { exit 1; }
sudo make install || { exit 1; }

sudo ldconfig || { exit 1; }

rm -rf $TMP_DIR

# Turn off echoing commands
set +x

echo -n "Checking whether cJSON is installed... "

if $PKG_CONFIG --exists libcjson; then
    echo "yes"
else
    echo "no"
    exit 1
fi

else

echo "yes"

fi
