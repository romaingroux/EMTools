# install the boost library

dir=$(pwd)
library_dir="$dir/lib/boost"
mkdir -p $library_dir

# download src
wget https://dl.bintray.com/boostorg/release/1.70.0/source/boost_1_70_0.tar.gz
tar -xzvf boost_1_70_0.tar.gz

cd boost_1_70_0/

# build and install
./bootstrap.sh --prefix=$library_dir
./b2 install link=static # program_options

# clean
cd ..
rm -r boost_1_70_0
rm    boost_1_70_0.tar.gz
