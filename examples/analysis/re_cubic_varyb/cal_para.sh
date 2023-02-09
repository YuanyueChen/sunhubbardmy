#!/bin/bash
#plquarray=$(awk 'BEGIN{for(i=0.00;i<=4.001;i+=4.00) printf("%6.2f",i)}')
#plquarray=$( echo "0.01 0.10 1.00 6.00")
#plquarray=$( echo "0.10 1.00 6.00")
plquarray=$( echo "0.00")
#plquarray=$( echo "99.00")
#uarray=$(awk 'BEGIN{for(i=28.00;i<=64.001;i+=8.00) printf("%6.2f",i)}')
#uarray=$( echo "0.10 1.00 6.00")
#uarray=$( echo "0.00")
#uarray=$( echo "4.00")
uarray=$( echo "4.00")
#betaarray=$( echo "0.2")
#betaarray=$( echo "2.000 2.200 2.500 2.800 3.300 4.000 5.000")
#betaarray=$( echo "1.800 2.400 3.000 3.600 4.200 4.800")
#betaarray=$( echo "2.000 2.100 2.200 2.350 2.500 2.650 2.850 3.100 3.350 3.650 4.000 4.450 5.000") # U=8
#betaarray=$( echo "2.500 2.650 2.850 3.100 3.350 3.650 4.000 4.450 5.000")  # U=12
#betaarray=$( echo "3.350 3.650 4.000 4.450 5.000 5.700 6.650 8.000 10.000") # U=4
#betaarray=$( echo "2.800")
betaarray=$(awk 'BEGIN{for(i=1.00;i<=10.001;i+=0.005) printf("%7.3f",i)}')
#Larray=$( echo "18 21")
#Larray=$( echo "8 12 16")
Larray=$( echo "16")
#Larray=$( echo "4 6 8 10 12 14 16")
#Larray=$( echo "4 6")
#Larray=$( echo "8 12 16 20")
#Larray=$( echo "4 8 12 16 20")
#Larray=$( echo "6 9 12 15 18")
alpha=0.00
theta=0.00
#nflrarray=$( echo "2 4 6 8")
#nflrarray=$( echo "2 4 6")
nflrarray=$( echo "2")
#nuarray=$(awk 'BEGIN{for(i=0.00;i<=4.001;i+=4.00) printf("%6.2f",i)}')
nuarray=$( echo "0.00")
#nuarray=$( echo "0.00 -1.00 -2.00 -3.00")
#nuarray=$( echo "0.00 -2.00")
#code="square_piflux_suN_dqmc_ft_hist"
#code="square_suN_dqmc_ft_hist"
#code="honeycomb_suN_dqmc_ft_hist"
#code="honeycomb_suN_dqmc_proj_hist"
#code="square_suN_dqmc_proj_hist"
#code="square_piflux_suN_dqmc_proj_hist"
#code="cubic_suN_ft"
#code="cubic_suN_ft_delay"
#code="cubic_suN_ft_brT"
code="cubic_suN_ft_brT_delay"
xmag=0.0
flux_x=0.0001
flux_y=0.0
dtau=0.05
rt=1.0
lprojplqu=F
lproju=F
ltau=T
nwrap=10
nsweep=10
nbin=25
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
