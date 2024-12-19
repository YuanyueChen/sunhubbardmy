#!/bin/bash
plquarray=$( echo "0.00")
# uarray=$(awk 'BEGIN{for(i=0.00;i<=0.001;i+=1.00) printf("%6.2f",i)}')
uarray=$( echo "-1.00")
# varray=$(awk 'BEGIN{for(i=1.00;i<=1.001;i+=1.00) printf("%6.2f",i)}')
varray=$( echo "0.00")
betaarray=$( echo "1.000")
# Larray=$(awk 'BEGIN{for(i=60.00;i<=60.001;i+=4.00) printf("%6.2f",i)}')
#Larray=$(awk 'BEGIN{for(i=20.00;i<=80.001;i+=6.00) printf("%6.2f",i)}')
# Larray=$( echo "38")
Larray=$( echo "44")
# Larray=$( echo "30")
alpha=0.00
theta=0.00
nflrarray=$( echo "2")
# nuarray=$( echo "0.00 -1.00 -2.00 -3.00")
nuarray=$( echo "-1.00")
# nuarray=$( echo "0.00")
# nuarray=$( echo "-0.50")
codearray=$( echo "proj_square_submatrixG_timing proj_square_submatrixLR_timing")
# codearray=$( echo "proj_square_submatrixG_timing proj_square_submatrixLR_timing proj_square_fast_timing")
# codearray=$( echo "proj_square_fast_timing")
# codearray=$( echo "proj_square_submatrixLR1-1_timing proj_square_submatrixLR1-3_timing")
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
nsweep=5
nbin=2
#nublock=0
#nublockarray=$(awk 'BEGIN{for(i=2.00;i<=2.001;i+=64.00) printf("%8.2f",i)}')
# nublockarray=$(awk 'BEGIN{for(i=32.00;i<=544.001;i+=64.00) printf("%8.2f",i)}')
# nublockarray=$( echo "8 16 32 96 160 224 288 352 416 480 544")
nublockarray=$( echo "4 8 16 32 96 160 224 288")
# nublockarray=$( echo "4")
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
echo " nublock = $nublockarray"  
echo " lprojplqu = " $lprojplqu
echo " lproju = " $lproju
echo " ltau = " $ltau
echo " nwrap = " $nwrap
echo " nsweep = " $nsweep
echo " nbin = " $nbin
echo " code = " $codearray

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

for nublocktmp in $nublockarray; do
    let num_nublock=num_nublock+1
done
echo " num_nublock = " $num_nublock


for vtmp in $varray; do
    let num_v=num_v+1
done
echo " num_V = " $num_v

for utmp in $uarray; do
    let num_u=num_u+1
done
echo " num_U = " $num_u
