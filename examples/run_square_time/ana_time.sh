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
for L  in ${Larray}; do
# Create a new data file for each code
datafile="time_${code}.b${beta}.U${u}.plqU${plqu}.V${rv}.nu${nu}.Nf${nflr}.L${L}.dat"
rm $WORKDIR/$datafile
  for nublock  in ${nublockarray}; do
       cd $WORKDIR
        maindir=${code}.b${beta}.U${u}.plqU${plqu}.V${rv}.nu${nu}.Nf${nflr}.L${L}.nublock${nublock}
         if [ -f $maindir/dqmc.out ]; then
           echo ${maindir}
           cd $maindir
           time=$(grep -E "The_time_of_update_ratio.*_update"  dqmc.out|awk -F: '{printf "%.2f\n", $2}')
           time1=$(grep -E "The_time_of_update_gfunc.*_update"  dqmc.out|awk -F: '{printf "%.2f\n", $2}')
           time2=$(echo "$time + $time1" | bc)
           time3=$(grep 'The_time_of_sweep_in_total'  dqmc.out|awk -F: '{printf "%.2f\n", $2}')
          #  echo ${nublock} ${time} ${time1} ${time2} ${time3} #>> $WORKDIR/$datafile
          #  echo $WORKDIR/$datafile
          #  echo -n > $datafile
           echo ${nublock} ${time} ${time1} ${time2} ${time3} >> $WORKDIR/$datafile
         fi
done
done
done
done
done
done
done
done
done
