#!/bin/bash

NP=4
ND=1

source cal_para.sh

WORKDIR="$PWD"
echo $WORKDIR
cd $WORKDIR
for beta in ${betaarray}; do
for u in ${uarray}; do
for rv in ${varray}; do
for plqu in ${plquarray}; do
for nu in ${nuarray}; do
for nflr in ${nflrarray}; do
  for L  in ${Larray}; do
  for nublock  in ${nublockarray}; do 
for code in ${codearray}; do
        EXE=$HOME/Projects/sun-hubbard/bin/${code}
        cd $WORKDIR
        #beta=$( echo "$L" |awk '{printf( "%4i", $1*2)}' | tr -d '[:space:]' )
        # maindir=${code}.b${beta}.U${u}.plqU${plqu}.V${rv}.Nf${nflr}.L${L}.nublock${nublock}
        maindir=${code}.b${beta}.U${u}.plqU${plqu}.V${rv}.nu${nu}.Nf${nflr}.L${L}.nublock${nublock}
        if [ ! -d $maindir ]; then
            mkdir $maindir
        fi
        cd $maindir

        if [ ! -f energy.bin ]; then
        #if [ ! -f dqmc.in ]; then
        echo "submit job $maindir ... "

        if [ -f conf_plqu.out ]; then
            cp conf_plqu.out conf_plqu.in
        fi
        if [ -f conf_u.out ]; then
            cp conf_u.out conf_u.in
        fi

        if [ -f conf_v.out ]; then
            cp conf_v.out conf_v.in
        fi

        #if [ -f energy.bin ]; then
        #    nbin=15
        #else
        #    nbin=20
        #fi

cat>dqmc.in<<endin
&model_para
rt = $rt,
la = $L,
lb = $L,
beta = $beta,
dtau = $dtau,
nu = $nu,
rhub = $u,
rhub_plq = $plqu,
rv = $rv
alpha = $alpha,
theta = $theta,
nflr = $nflr,
lprojplqu = $lprojplqu,
lproju = $lproju,
xmag = $xmag,
flux_x = $flux_x,
flux_y = $flux_y,
rndness = $rndness,
/
&ctrl_para
ltau = $ltau,
dyntau = $dyntau,
obs_eqt_mid_len = $obs_eqt_mid_len,
nwrap = $nwrap,
nsweep = $nsweep,
nbin = $nbin,
nublock = $nublock
/
endin

if [ -f dqmc.out ]; then
    cat dqmc.out >> dqmc.out.bak
fi
mpirun -np $NP $EXE >logs 2>&1

       #fi
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