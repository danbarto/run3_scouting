#!/usr/bin/env bash

#tar cvz -f package.tar.gz CMSSW_12_4_12/ pythia8/ psets/
wget https://pythia.org/download/pythia83/pythia8309.tgz
tar cvz -f package.tar.gz pythia8309.tgz psets/ signal/ data/
