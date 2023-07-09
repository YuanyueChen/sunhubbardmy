#!/bin/bash

NP=16
ND=1
NT=4 # number of threads

source cal_para.sh

WORKDIR="$PWD"
echo $WORKDIR
EXE=$HOME/mywork/sun-hubbard/bin/${code}
cd $WORKDIR
for beta in ${betaarray}; do
for u in ${uarray}; do
for rv in ${varray}; do
for plqu in ${plquarray}; do
for nu in ${nuarray}; do
for nflr in ${nflrarray}; do
  for L  in ${Larray}; do
        cd $WORKDIR
        #beta=$( echo "$L" |awk '{printf( "%4i", $1*2)}' | tr -d '[:space:]' )
        maindir=${code}.b${beta}.U${u}.plqU${plqu}.V${rv}.Nf${nflr}.L${L}
        echo "submit job $maindir ... "
        if [ ! -d $maindir ]; then
            mkdir $maindir
        fi
        cd $maindir

        if [ ! -f energy.bin ]; then
        #if [ ! -f dqmc.in ]; then

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
rv = $rv,
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

cat>${maindir}.sub<<endsub
#!/bin/bash
#SBATCH -J $maindir
#SBATCH -p 64c512g
#SBATCH --exclusive
#SBATCH -n $NP
#SBATCH --ntasks-per-node=$NP
#SBATCH -t 168:00:00
#SBATCH -o job.%j.%N.out
#SBATCH -e job.%j.%N.err

module load intel-oneapi-compilers
module load intel-oneapi-mkl
module load intel-oneapi-mpi
export OMP_NUM_THREADS=$NT

mpirun -np $NP $EXE
endsub

       #nb=$( wc -l energy.bin |awk '{print $1}' )
       #if [ "$nb" = "30" ]; then
       #if [ ! -f energy.bin ]; then
          #rm -f job*.out
          #rm -f job*.err
          cat dqmc.out >> dqmc.out.bak
          sbatch ${maindir}.sub
       #fi
       fi
  done
done
done
done
done
done
