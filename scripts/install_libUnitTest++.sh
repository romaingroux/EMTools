# install the boost library
library_dir='lib/UnitTest++'
lib_dir="$library_dir/lib"
include_dir="$library_dir/include"

# download src
git clone https://github.com/unittest-cpp/unittest-cpp.git

mkdir -p $library_dir
mkdir -p $lib_dir
mkdir -p $include_dir

cd unittest-cpp/

# install
cmake3 . && make
find UnitTest++/ -name "*.cpp" -type f -delete
mv ./libUnitTest++.a ../$lib_dir/
mv UnitTest++/*      ../$include_dir/

# clean
cd ..
rm -rf unittest-cpp
