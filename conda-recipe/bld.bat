rmdir build /s /q
"%PYTHON%" setup.py install --single-version-externally-managed --record=record.txt

cd build
cmake .. -G "NMake Makefiles" ^
  -DCMAKE_INSTALL_PREFIX=%PREFIX%/Library ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DBUILD_PYTHON_INTERFACE=OFF ^
  -DCMAKE_PREFIX_PATH=%PREFIX%

nmake
if errorlevel 1 exit 1
nmake install
if errorlevel 1 exit 1
