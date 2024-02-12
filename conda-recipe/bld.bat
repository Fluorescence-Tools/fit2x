git submodule update --recursive --init --remote

md build
cd build

REM Call Python with the --version flag to get the version information
for /f "tokens=2 delims= " %%v in ('%PYTHON% --version 2^>^&1') do set PYTHON_VERSION=%%v

REM Extract only the numeric part of the version
for /f "tokens=1-3 delims=." %%a in ("%PYTHON_VERSION%") do set PYTHON_VERSION_NUMERIC=%%a.%%b.%%c

cmake .. -G "NMake Makefiles" ^
 -DCMAKE_INSTALL_PREFIX="%LIBRARY_PREFIX%" ^
 -DCMAKE_PREFIX_PATH="%PREFIX%" ^
 -DBUILD_PYTHON_INTERFACE=ON ^
 -DCMAKE_BUILD_TYPE=Release ^
 -DCMAKE_LIBRARY_OUTPUT_DIRECTORY="%SP_DIR%" ^
 -DPYTHON_VERSION="%PYTHON_VERSION_NUMERIC%" ^
 -DCMAKE_SWIG_OUTDIR="%SP_DIR%"

nmake
nmake install

cd ..
rmdir build /s /q
