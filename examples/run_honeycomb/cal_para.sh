#!/bin/bash
plquarray=$( echo "0.00")
uarray=$(awk 'BEGIN{for(i=0.00;i<=0.001;i+=1.00) printf("%6.2f",i)}')
#uarray=$( echo "2.00 3.00 4.00 5.00")
varray=$(awk 'BEGIN{for(i=1.00;i<=1.001;i+=1.00) printf("%6.2f",i)}')
#varray=$( echo "2.00 3.00 4.00 5.00")
betaarray=$( echo "40.000")
Larray=$( echo "3")
alpha=0.00
theta=0.00
nflrarray=$( echo "2")
nuarray=$( echo "0.00")
code="proj_FFT_delay"
xmag=0.0
flux_x=0.0
flux_y=0.0
rndness=0.00001
dtau=0.1
rt=1.0
lprojplqu=F
lproju=F
ltau=F
dyntau=0.0
nwrap=10
nsweep=20
nbin=20
nublock=0
obs_eqt_mid_len=21
echo " rt = " $rt
echo " L = " $Larray
echo " beta = " $betaarray
echo " dtau = " $dtau
echo " U = " $uarray
echo " plqU = " $plquarray
echo " alpha = $alpha" 
echo " theta = $theta" 
echo " nflr = $nflrarray" 
echo " nu = $nuarray" 
echo " lprojplqu = " $lprojplqu
echo " lproju = " $lproju
echo " ltau = " $ltau
echo " nwrap = " $nwrap
echo " nsweep = " $nsweep
echo " nbin = " $nbin
echo " code = " $code

for plqutmp in $plquarray; do
    let num_plqu=num_plqu+1
done
echo " num_plqu = " $num_plqu

for betatmp in $betaarray; do
    let num_beta=num_beta+1
done
echo " num_beta = " $num_beta

for Ltmp in $Larray; do
    let num_L=num_L+1
done
echo " num_L = " $num_L

for nutmp in $nuarray; do
    let num_nu=num_nu+1
done
echo " num_nu = " $num_nu

for vtmp in $varray; do
    let num_v=num_v+1
done
echo " num_V = " $num_v

for utmp in $uarray; do
    let num_u=num_u+1
done
echo " num_U = " $num_u
