git submodule update --recursive --init --remote
#$PYTHON setup.py install --single-version-externally-managed --record=record.txt
mkdir build
cd build

cmake \
 -DCMAKE_INSTALL_PREFIX="$PREFIX" \
 -DCMAKE_PREFIX_PATH="$PREFIX" \
 -DBUILD_PYTHON_INTERFACE=ON \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_LIBRARY_OUTPUT_DIRECTORY="$SP_DIR" \
 -DCMAKE_SWIG_OUTDIR="$SP_DIR" \
 -DPYTHON_VERSION=$(python -c 'import platform; print(platform.python_version())')\
 ..

make
make install
