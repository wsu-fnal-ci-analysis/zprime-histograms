#!/bin/bash


mkdir -p /lustre/cms/store/user/defilip/ZprimeAnalysis/80X/jobdir
mkdir -p /lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histodir

echo "Running ZprimeMuMu Analysis with executables RunZprimeMuMuAnalysis"
source /cvmfs/cms.cern.ch/cmsset_default.sh

export LD_LIBRARY_PATH=lib:$LD_LIBRARY_PATH
export PATH=path:$PATH


if [ -d "$_CONDOR_SCRATCH_DIR" ]; then
    workdir=`echo $_CONDOR_SCRATCH_DIR`;
    cd ${workdir};
else 
    workdir=`echo $PWD`;
    cd ${workdir};
fi


savedir=`echo /lustre/cms/store/user/defilip/ZprimeAnalysis/80X/histodir`

echo "Working dir is $workdir"
#echo "Executable dir is $exedir"
echo "Saving dir is $savedir"

echo "Compiling the macros"
bash compileZprimeMuMuAnalysis.sh


./RunZprimeMuMuAnalysis which sig_input.txt 1 bkg_input.txt 1 data_input.txt 1 site year mc >& ${workdir}/RunZprimeMuMuAnalysis.log

mv -f ${workdir}/RunZprimeMuMuAnalysis.log /lustre/cms/store/user/defilip/ZprimeAnalysis/80X/jobdir/output.log
mv -f ${workdir}/ZprimeToMuMu_13TeV.root     ${savedir}/output.root
mv -f ${workdir}/ZprimeToMuMu_13TeV_cand.txt ${savedir}/output_cand.txt

