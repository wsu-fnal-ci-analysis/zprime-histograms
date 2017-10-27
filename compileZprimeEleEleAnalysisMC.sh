#!/bin/bash

g++ -I $ROOTSYS/include RunZprimeEleEleAnalysisMC.C ZprimeEleElePatMiniAodNewMC.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o RunZprimeEleEleAnalysisMC
