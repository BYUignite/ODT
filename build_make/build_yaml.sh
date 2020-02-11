#!/bin/bash

rm -rf ../source/yaml/include
rm -rf ../source/yaml/lib
rm -rf ../source/yaml/yaml-cpp
cd ../source/yaml
git clone https://github.com/jbeder/yaml-cpp.git
cd yaml-cpp
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../..
make
make install
