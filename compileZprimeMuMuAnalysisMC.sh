#!/bin/bash

g++ -I $ROOTSYS/include RunZprimeMuMuAnalysisMC.C ZprimeMuMuPatMiniAodNewMC.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o RunZprimeMuMuAnalysisMC
