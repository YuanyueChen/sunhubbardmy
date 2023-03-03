#!/bin/bash

source cal_para.sh
WORKDIR="$PWD"
datadir=$WORKDIR/../../run_cubic/
#datadir=/dssg/home/acct-phyxxy/phyxxy/mywork/sun-hubbard/run_cubic
#datadir="$PWD"
echo $WORKDIR
cd  $WORKDIR
rm time.dat
for u in ${uarray}; do
for plqu in ${plquarray}; do
for nu in ${nuarray}; do
for nflr in ${nflrarray}; do
for beta in ${betaarray}; do
for L  in ${Larray}; do
       cd $datadir
        maindir=${code}.b${beta}.dtau${dtau}.t${rt}.U${u}.plqU${plqu}.alpha${alpha}.theta${theta}.Nf${nflr}.nu${nu}.L${L}.proj${lprojplqu}${lproju}
        #echo ${maindir}
         if [ -f $maindir/dqmc.out ]; then
           cd $maindir
           time=$(grep 'The time of update'  dqmc.out|awk -F: '{printf "%.2f\n", $2}')
           time1=$(grep 'Total time spent'  dqmc.out|awk -F: '{printf "%.2f\n", $2}')
           echo ${L} ${time} ${time1} >> $WORKDIR/time.dat
         fi
done
done
done
done
done
done
