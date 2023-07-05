#!/bin/bash
plquarray=$( echo "0.00")
uarray=$(awk 'BEGIN{for(i=2.00;i<=5.001;i+=1.00) printf("%6.2f",i)}')
#uarray=$( echo "2.00 3.00 4.00 5.00")
betaarray=$( echo "40.000")
Larray=$( echo "3 6")
alpha=0.00
theta=0.00
nflrarray=$( echo "2")
nuarray=$( echo "0.00")
code="honeycomb_suN_proj_FFT_delay"
xmag=0.0
flux_x=0.0001
flux_y=0.0
dtau=0.1
rt=1.0
lprojplqu=F
lproju=F
ltau=F
nwrap=10
nsweep=20
nbin=20
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
