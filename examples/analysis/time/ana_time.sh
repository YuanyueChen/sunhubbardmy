#!/bin/bash

source cal_para.sh
WORKDIR="$PWD"
echo $WORKDIR
cd  $WORKDIR
for code in ${codearray}; do
for beta in ${betaarray}; do
for u in ${uarray}; do
for rv in ${varray}; do
for plqu in ${plquarray}; do
for nu in ${nuarray}; do
for nflr in ${nflrarray}; do
# Create a new data file for each code
datafile="time_${code}.b${beta}.U${u}.plqU${plqu}.V${rv}.Nf${nflr}.dat"
echo -n > $datafile
rm $datafile
  for L  in ${Larray}; do
       cd $WORKDIR
        maindir=${code}.b${beta}.U${u}.plqU${plqu}.V${rv}.Nf${nflr}.L${L}
        echo ${maindir}
         if [ -f $maindir/dqmc.out ]; then
           cd $maindir
           time=$(grep -E "The_time_of_.*_update"  dqmc.out|awk -F: '{printf "%.2f\n", $2}')
           time1=$(grep 'The_time_of_sweep_in_total'  dqmc.out|awk -F: '{printf "%.2f\n", $2}')
           echo ${L} ${time} ${time1} >> $WORKDIR/$datafile
         fi
done
done
done
done
done
done
done
done
