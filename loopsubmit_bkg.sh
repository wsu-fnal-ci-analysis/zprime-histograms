#!/bin/bash

if [ "$1" == "" ] || [ "$2" == "" ] || [ "$4" == "" ] || [ "$4" == "" ] || [ "$5" == "" ]; then
    echo "Please provide arguments to the script: site configuration, data type and MC type, Ele/Mu for di electron vs dimuon"
    echo "Usage bash loopsubmit_bkg_<finalstate>.sh <arg1> <arg2> <arg3> <arg4> <arg5>"
    exit
fi

echo "$1 configuration";
echo "$2 data"
echo "$3 simulation"
echo "$4 site"
echo "$5 lflav"

conf=${1}
data=${2}
simu=${3}
site=${4}
lflav=${5}

SCERN="CERN";
SFNAL="FNAL";
SDESY="DESY";
SBARI="BARI";

mkdir -p jobs;

###### Background
n=0;
m=0;

subflav=
if [ "${lflav}" = "Ele" ]
then
    subflav="ee"
elif [ "${lflav}" = "Mu" ]
then
    subflav="mumu"
else
    echo "Invalid lepton flavour specified: ${lflav} is not one of 'Mu' or 'Ele'"
    exit 1
fi

echo "Reading bkg_input_${simu}_${subflav}_AN.txt file"

cp -f bkg_input_${simu}_${subflav}_AN.txt bkg_input.txt
nlines=`wc -l bkg_input_${simu}_${subflav}_AN.txt | awk '{print $1}'`;

while [ $n -lt ${nlines} ]; do
  (( n = n + 1 ))
  (( m = ${nlines} - n ))
  echo $n $m
  mkdir -p BkgCards${simu}_${subflav}
  rm -f BkgCards${simu}_${subflav}/bkg_input_${n}.txt
  cat bkg_input.txt | head -1 > BkgCards${simu}_${subflav}/bkg_input_${n}.txt
  samplename=`cat BkgCards${simu}_${subflav}/bkg_input_${n}.txt | awk '{print $1}'`
  echo $samplename

  skip=0
  res=$(echo ${samplename} | fgrep "#")
  if [ $res ]
  then
      echo "Ignoring commented out ${samplename}"
      skip=1
  elif [ ${lflav} = "Ele" ]
  then
      res=$(echo ${samplename} | fgrep "MuMu")
      if [ $res ]
      then
	  echo "skipping ${samplename} for di-electron analysis"
	  skip=1
      fi
  elif [ ${lflav} = "Mu" ]
  then
      res=$(echo ${samplename} | fgrep "EE")
      if [ $res ]
      then
	  echo "skipping ${samplename} for di-muon analysis"
	  skip=1
      fi
  fi

  rm -f jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
  rm -f jobs/condor_template.cfg
  cat bkg_input.txt | tail -n $m >  bkg_input_tmp.txt
  mv  bkg_input_tmp.txt bkg_input.txt

  if [ "${skip}" = "1" ]
  then
      rm -f BkgCards${simu}_${subflav}/bkg_input_${n}.txt
      continue
  fi
  
  if [ ${site} = ${SCERN} ]; then
      cat submit_Zprime${lflav}${lflav}Analysis_CERN.sh | \
	  sed "s?which?bkg?g" | \
	  sed "s?site?${site}?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?Zprime${lflav}${lflav}Analysis?RunZprime${lflav}${lflav}Analysis?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?bkg_input.txt?BkgCards${simu}_${subflav}/bkg_input_${n}.txt?g" | \
	  sed "s?s.log?s_${samplename}.log?g" > jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
  elif  [ ${site} = ${SFNAL} ]; then
      cat submit_Zprime${lflav}${lflav}Analysis_FNAL.sh | \
	  sed "s?Zprime${lflav}${lflav}Analysis?Zprime${lflav}${lflav}AnalysisMC?g" | \
	  sed "s?path?$PATH?g"  | \
	  sed "s?lib?$LD_LIBRARY_PATH?g" | \
	  sed "s?which?bkg?g" | \
	  sed "s?site?${site}?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?bkg_input.txt?bkg_input_${n}.txt?g" | \
	  sed "s?sMC.log?s_${samplename}.log?g" > jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
      cat condor_templateMC.cfg  | \
	  sed "s?submit_ZprimeAnalysis_FNAL?submit_Zprime${lflav}${lflav}Analysis_${samplename}?g" | \
	  sed "s?ZprimeAnalysis?Zprime${lflav}${lflav}Analysis?g" | \
	  sed "s?ZprimePat?Zprime${lflav}${lflav}Pat?g" | \
	  sed "s?sig_input_h150.txt?BkgCards${simu}_${subflav}/bkg_input_${n}.txt?g" | \
	  sed "s?mail?`whoami`?g" > jobs/condor_Zprime${lflav}${lflav}Analysis_${samplename}.cfg
  elif  [ ${site} = ${SDESY} ]; then
      cat submit_Zprime${lflav}${lflav}Analysis_DESY.sh | \
	  sed "s?which?bkg?g" | \
	  sed "s?site?${site}?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?bkg_input.txt?BkgCards${simu}_${subflav}/bkg_input_${n}.txt?g" | \
	  sed "s?s.log?s_${samplename}.log?g" > jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
  elif  [ ${site} = ${SBARI} ]; then
      cat submit_Zprime${lflav}${lflav}Analysis_BARI.sh  | \
	  sed "s?path?$PATH?g"  | \
	  sed "s?lib?$LD_LIBRARY_PATH?g" | \
	  sed "s?which?bkg?g" | \
	  sed "s?site?${site}?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?bkg_input.txt?bkg_input_${n}.txt?g" | \
	  sed "s?s.log?s_${samplename}.log?g" > jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
      cat condor_templateMC.cfg  | \
	  sed "s?submit_ZprimeAnalysis_BARI?submit_Zprime${lflav}${lflav}Analysis_${samplename}?g" | \
	  sed "s?ZprimeAnalysis?Zprime${lflav}${lflav}Analysis?g" | \
	  sed "s?ZprimePat?Zprime${lflav}${lflav}Pat?g" | \
	  sed "s?sig_input_h150.txt?BkgCards${simu}_${subflav}/bkg_input_${n}.txt?g" | \
	  sed "s?mail?`whoami`?g" > jobs/condor_Zprime${lflav}${lflav}Analysis_${samplename}.cfg
  else
      cat submit_Zprime${lflav}${lflav}Analysis.sh | \
	  sed "s?which?bkg?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?bkg_input.txt?BkgCards${simu}_${subflav}/bkg_input_${n}.txt?g" | \
	  sed "s?s.log?s_${samplename}.log?g" > jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
  fi

  chmod u+xr jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh

  cd jobs

  if [ ${site} = ${SCERN} ]; then
      echo "Submitting jobs via LSF at CERN"
      bsub -q 8nh  submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
  elif  [ ${site} = ${SFNAL} ]; then
      echo "Submitting jobs via CONDOR at FNAL"
      condor_submit  condor_Zprime${lflav}${lflav}Analysis_${samplename}.cfg
  elif  [ ${site} = ${SDESY} ]; then
      echo "Submitting jobs via SGE"
      qsub  -l h_rt="10:00:00" submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
   elif  [ ${site} = ${SBARI} ]; then
      echo "Submitting jobs via CONDOR at BARI"
      condor_submit  -name ettore condor_Zprime${lflav}${lflav}Analysis_${samplename}.cfg
  else
      echo "Submitting jobs via PBS"
      qsub -q local submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
  fi
  cd ..
done

