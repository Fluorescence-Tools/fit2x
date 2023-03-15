git submodule update --recursive --init --remote

md build
cd build

cmake .. -G "NMake Makefiles" ^
 -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
 -DCMAKE_PREFIX_PATH="%PREFIX%" ^
 -DBUILD_PYTHON_INTERFACE=ON ^
 -DCMAKE_BUILD_TYPE=Release ^
 -DCMAKE_LIBRARY_OUTPUT_DIRECTORY="%SP_DIR%" ^
 -DCMAKE_SWIG_OUTDIR="%SP_DIR%"

nmake
nmake install

cd ..
rmdir build /s /q
