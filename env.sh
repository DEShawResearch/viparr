garden env-keep-only VIPARR_FFPATH --user
garden load `cat $(dirname $0)/../share/modules.txt`
garden prepend-path PYTHONPATH `readlink -f $(dirname $0)/../lib/python`

