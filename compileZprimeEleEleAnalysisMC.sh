#!/bin/bash

g++ -g3 -O0 -fno-inline -I $ROOTSYS/include RunZprimeEleEleAnalysisMC.C ZprimeEleElePatMiniAodNewMC.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o RunZprimeEleEleAnalysisMC
