#!/bin/bash

dir=${1:-plots}
outdir=${dir}/cmsswval-pu70-extrectracks
base=SNB_CMSSW_PU70

echo "Moving plots and text files locally to ${outdir}"

mkdir -p ${outdir}
mv ${base}_*.png ${outdir}
for build in BH STD CE
do
    vbase=validation_${base}_${build}
    mv ${vbase}/totals_${vbase}_cmssw.txt ${outdir}
done

host=kmcdermo@lxplus.cern.ch
whost=${host}":~/www"
echo "Moving plots and text files remotely to ${whost}"
scp -r ${dir} ${whost}

echo "Executing remotely ./makereadable.sh ${outdir}"
ssh ${host} bash -c "'
cd www
./makereadable.sh ${outdir}
exit
'"

echo "Removing local files"
for build in BH STD CE
do
    testbase=${base}_${build}
    rm -rf validation_${testbase}
    rm -rf log_${testbase}_NVU8int_NTH24_cmsswval.txt 
done

rm -rf ${dir}
