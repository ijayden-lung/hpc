#!/bin/bash


# Use OpenBLAS with 1 thread only as it seems to be using too many
# on the CIs apparently.
export OPENBLAS_NUM_THREADS=1

# # Depending on our platform, shared libraries end with either .so or .dylib
if [[ $(uname) == 'Darwin' ]]; then
    export LDFLAGS="$LDFLAGS -undefined dynamic_lookup"
else
    export LDFLAGS="$LDFLAGS -shared"
fi

if [ "$blas_impl" = "mkl" ]; then

cat > site.cfg <<EOF
[mkl]
library_dirs = $PREFIX/lib
include_dirs = $PREFIX/include
mkl_libs = mkl_rt
lapack_libs =

EOF

else # Non-MKL

cat > site.cfg <<EOF
[DEFAULT]
library_dirs = $PREFIX/lib
include_dirs = $PREFIX/include

[atlas]
atlas_libs = openblas
libraries = openblas

[openblas]
libraries = openblas
library_dirs = $PREFIX/lib
include_dirs = $PREFIX/include

EOF

fi

$PYTHON setup.py install --single-version-externally-managed --record=record.txt
