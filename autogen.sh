#!/bin/sh
mkdir -p m4
mkdir -p config
autoreconf --force --install -I config -I m4
