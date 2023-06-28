#!/bin/bash

cwd=$(pwd)

# Resolve dependencies

echo 'Pip installables (scipy, numpy, mappy, Cython)'
python3 -m pip install --user --upgrade scipy numpy mappy Cython editdistance

echo 'conk'
python3 -m pip install --user --upgrade wheel setuptools
git clone https://github.com/rvolden/conk
cd conk && make
cd $cwd

echo 'Racon'
git clone --recursive https://github.com/isovic/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd $cwd

echo 'abpoa'
wget https://github.com/yangao07/abPOA/releases/download/v1.4.1/abPOA-v1.4.1.tar.gz
tar -zxvf abPOA-v1.4.1.tar.gz
cd abPOA-v1.4.1; make
cd $cwd
rm abPOA-v1.4.1.tar.gz

echo 'blat'
mkdir blat
cd blat
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat
chmod +x blat
cd $cwd



echo 'Done'
