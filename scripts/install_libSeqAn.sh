# install the header only SeqAn library

library_dir='lib/seqan'

# clone git
git clone https://github.com/seqan/seqan.git
cd seqan
#install
mkdir -p           ../$library_dir
## header files
mv * ../$library_dir
cd ..
# clean
rm -rf seqan
