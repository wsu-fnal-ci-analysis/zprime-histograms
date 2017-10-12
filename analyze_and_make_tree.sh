#!/usr/bin/env zsh

# Run this file to analyze a given tree

if [ -z "$1" ]
    then
        OUTFILENAME="CMSSW803-Analyze_ZtoMuMu_13TeV_M5000_miniaod_filledhistos.root"
        echo $0: optional usage: $0 outputfilename inputfilename
        #exit 3
    else
        OUTFILENAME=$1
        if [ -z "$2" ]
            then
                INFILENAME="/nfs/dust/cms/user/sonnevej/CMSSW803_MC_DYtoMuMu_13TeV_pattuple200.root"
                echo $0: optional usage: $0 outputfilename inputfilename
                #echo $0: usage: $0 outputfilename
                #exit 3
            else
                INFILENAME=$2
        fi
fi



if [[ -a $OUTFILENAME ]]
    then
        echo will not analyze your tree $INFILENAME
        echo outputfile $OUTFILENAME already exists
        echo please remove or rename $OUTFILENAME
        exit 1
    else
        echo "will write output to" $OUTFILENAME
        echo will now analyze your tree $INFILENAME
fi
export QUOTES='"'
export WRITE_TO_FILE=$QUOTES$OUTFILENAME$QUOTES
export ANALYZE_FILE=$QUOTES$INFILENAME$QUOTES



root -b -l << EOF
Char_t name[300];
sprintf(name,${ANALYZE_FILE});
TFile *file0 = TFile::Open(name)
TTree *tree3 = (TTree*)file0->Get("tree");
.L ZprimeMuMuPatMiniAod.C+
ZprimeMuMuPatMiniAod b(tree3);
b.Loop(${WRITE_TO_FILE});
.q
EOF


if [[ -a $OUTFILENAME ]]
    then
        echo file $OUTFILENAME has been created
    else
        echo did not succeed to write to $OUTFILENAME "-- Sorry!"
        echo "Sorry!"
fi
