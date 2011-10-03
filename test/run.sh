#!/usr/bin/env bash

eval "$(/usr/bin/garden-exec)"
set -e

garden env-keep-only PREFIX
garden load $(cat $(dirname $0)/../modules.txt)
PREFIX=${PREFIX:-build}
garden prepend-path PATH $PREFIX/bin
garden prepend-path PYTHONPATH $PREFIX/lib/python
exec pytest "$@"

