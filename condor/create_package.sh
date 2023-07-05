#!/bin/bash

SCRAMARCH=slc7_amd64_gcc10
CMSSWVERSION=CMSSW_12_6_0_patch1

rm tmp_create_package/ -rf
mkdir -p tmp_create_package
cd tmp_create_package

mkdir -p ScoutingRun3
cp ../*.txt ../*.py ../*.json ScoutingRun3/. -r # Copy relevant folders

tar -chJf package.tar.gz ScoutingRun3

mv package.tar.gz ../.
#rm -rf tmp_create_package # Leave the package be to able to easily check what is done
