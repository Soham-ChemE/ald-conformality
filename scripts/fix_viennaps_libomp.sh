#!/bin/sh
# macOS: the ViennaPS and ViennaLS PyPI wheels each bundle their own libomp.dylib. Two OpenMP runtimes
# in one process segfault inside the first parallel region (e.g. MakeTrench). Upstream issue:
# https://github.com/ViennaTools/ViennaPS/issues/167. Fix: make ViennaPS load ViennaLS's copy.
set -e
VENV="${1:-.venv}"
SP="$("$VENV/bin/python" -c 'import sysconfig; print(sysconfig.get_paths()["purelib"])')"
cd "$SP/viennaps/.dylibs"
if [ ! -L libomp.dylib ]; then
  mv libomp.dylib libomp.dylib.bundled
  ln -s ../../viennals/.dylibs/libomp.dylib libomp.dylib
  echo "linked viennaps libomp -> viennals libomp"
else
  echo "already linked"
fi
