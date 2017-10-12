#!/bin/bash

g++ -I $ROOTSYS/include RunZprimeMuMuAnalysis.C ZprimeMuMuPatMiniAodNewData.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o RunZprimeMuMuAnalysis
