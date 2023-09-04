#!/usr/bin/env bash
# NOTE: in mingw64 it is not LD_LIBRARY_PATH, but path that
# must be appended to find dynamic libraries!
export FAKEROOT="${HOME}/fakeroot"
export PATH="${FAKEROOT}/lib:${PATH}"
export PATH="${FAKEROOT}/casadi:${PATH}"

make

./hs107.exe
./gri30.exe

make clean
