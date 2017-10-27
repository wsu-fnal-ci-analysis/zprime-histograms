#!/bin/bash

g++ -I $ROOTSYS/include RunZprimeEleEleAnalysis.C ZprimeEleElePatMiniAodNewData.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o RunZprimeEleEleAnalysis
