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

echo "Reading data_input_${data}_AN.txt file"

cp -f data_input_${data}_AN.txt data_input.txt
nlines=`wc -l data_input_${data}_AN.txt | awk '{print $1}'`;

while [ $n -lt ${nlines} ]; do
  (( n = n + 1 ))
  (( m = ${nlines} - n ))
  echo $n $m
  mkdir -p DataCards${data}
  rm -f DataCards${data}/data_input_${n}.txt
  cat data_input.txt | head -1 > DataCards${data}/data_input_${n}.txt
  samplename=`cat DataCards${data}/data_input_${n}.txt | awk '{print $1}'`
  echo $samplename
  cat data_input.txt | tail -n $m >  data_input_tmp.txt
  mv  data_input_tmp.txt data_input.txt
  rm -f jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
  rm -f jobs/condor_template.cfg

  if [ ${site} = ${SCERN} ]; then
      cat submit_Zprime${lflav}${lflav}Analysis_CERN.sh | \
	  sed "s?which?data?g" | \
	  sed "s?site?${site}?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?Zprime${lflav}${lflav}Analysis?RunZprime${lflav}${lflav}Analysis?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?data_input.txt?DataCards${data}/data_input_${n}.txt?g" | \
	  sed "s?s.log?s_${samplename}.log?g" > jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
  elif  [ ${site} = ${SFNAL} ]; then
      cat submit_Zprime${lflav}${lflav}Analysis_FNAL.sh  | \
	  sed "s?path?$PATH?g"  | \
	  sed "s?lib?$LD_LIBRARY_PATH?g" | \
	  sed "s?which?data?g" | \
	  sed "s?site?${site}?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?data_input.txt?data_input_${n}.txt?g" | \
	  sed "s?s.log?s_${samplename}.log?g" > jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
      cat condor_template.cfg  | \
	  sed "s?submit_ZprimeAnalysis_FNAL?submit_Zprime${lflav}${lflav}Analysis_${samplename}?g" | \
	  sed "s?ZprimeAnalysis?Zprime${lflav}${lflav}Analysis?g" | \
	  sed "s?ZprimePat?Zprime${lflav}${lflav}Pat?g" | \
	  sed "s?sig_input_h150.txt?DataCards${data}/data_input_${n}.txt?g" | \
	  sed "s?mail?`whoami`?g" > jobs/condor_Zprime${lflav}${lflav}Analysis_${samplename}.cfg
  elif  [ ${site} = ${SDESY} ]; then
      cat submit_Zprime${lflav}${lflav}Analysis_DESY.sh | \
	  sed "s?which?data?g" | \
	  sed "s?site?${site}?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?data_input.txt?DataCards${data}/data_input_${n}.txt?g" | \
	  sed "s?s.log?s_${samplename}.log?g" > jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
  elif  [ ${site} = ${SBARI} ]; then
      cat submit_Zprime${lflav}${lflav}Analysis_BARI.sh  | \
	  sed "s?path?$PATH?g"  | \
	  sed "s?lib?$LD_LIBRARY_PATH?g" | \
	  sed "s?which?data?g" | \
	  sed "s?site?${site}?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?data_input.txt?data_input_${n}.txt?g"  | \
	  sed "s?s.log?s_${samplename}.log?g" > jobs/submit_Zprime${lflav}${lflav}Analysis_${samplename}.sh
      cat condor_template.cfg  | \
	  sed "s?submit_ZprimeAnalysis_BARI?submit_Zprime${lflav}${lflav}Analysis_${samplename}?g" | \
	  sed "s?ZprimeAnalysis?Zprime${lflav}${lflav}Analysis?g" | \
	  sed "s?ZprimePat?Zprime${lflav}${lflav}Pat?g" | \
	  sed "s?sig_input_h150.txt?DataCards${data}/data_input_${n}.txt?g" | \
	  sed "s?mail?`whoami`?g" > jobs/condor_Zprime${lflav}${lflav}Analysis_${samplename}.cfg
  else
      cat submit_Zprime${lflav}${lflav}Analysis.sh | \
	  sed "s?which?data?g" | \
	  sed "s?mc?${simu}?g" | \
	  sed "s?year?${data}?g" | \
	  sed "s?jobdir?jobs/jobsZprime${lflav}${lflav}?g" | \
	  sed "s?histodir?histos/histosZprime${lflav}${lflav}?g" | \
	  sed "s?output?output_${samplename}?g" | \
	  sed "s?data_input.txt?DataCards${data}/data_input_${n}.txt?g" | \
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

