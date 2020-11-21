#!/bin/sh
eval "$(/usr/bin/garden-exec)" 
set -e

garden env-keep-only HOME

PYTHON3=desres-python/3.6.6-04c7
SCONS=scons/3.1.2-01c7
BOOST=boost/1.57.0-02c7
PYBIND11=pybind11/2.3.0-01c7

MSYS_PREFIX=$HOME/git/gerrit/msys/build

loadmodules() {
    garden load \
        $PYTHON3/bin \
        $SCONS/bin \
        $BOOST/lib \
        $PYBIND11/lib \
        enscons/0.23.0-01c7/lib-python37 \

    garden prepend-path DESRES_MODULE_CXXFLAGS "-fpermissive -I$MSYS_PREFIX/include"
    garden prepend-path DESRES_MODULE_LDFLAGS "-L$MSYS_PREFIX/lib -Wl,-rpath,$MSYS_PREFIX/lib"
    garden prepend-path PYTHONPATH $(readlink -f $(dirname $0))    # for sconsutils
    export PYTHONVER=36
}

loadmodules
BUILD_WHEEL=1 BUILD_WHEEL_VERSION=$PYTHONVER DESRES_LOCATION= scons "$@"

version=$(src/version.py)
garden load auditwheel/3.1.1-01c7/bin
garden load patchelf/0.9-01c7/bin
auditwheel repair \
    --plat manylinux2014_x86_64 \
    -w build/wheel \
    build/wheel/viparr-${version}-cp${PYTHONVER}-cp${PYTHONVER}m-linux_x86_64.whl

