#!/bin/bash
#plquarray=$(awk 'BEGIN{for(i=0.00;i<=4.001;i+=4.00) printf("%6.2f",i)}')
#plquarray=$( echo "0.01 0.10 1.00 6.00")
#plquarray=$( echo "0.10 1.00 6.00")
plquarray=$( echo "0.00")
#plquarray=$( echo "99.00")
#uarray=$(awk 'BEGIN{for(i=28.00;i<=64.001;i+=8.00) printf("%6.2f",i)}')
#uarray=$( echo "0.10 1.00 6.00")
#uarray=$( echo "0.00")
uarray=$( echo "4.00")
#uarray=$( echo "10.00")
#betaarray=$( echo "0.2")
#betaarray=$( echo "2.000 2.100 2.200 2.350 2.500 2.650 2.850 3.100 3.350 3.650 4.000 4.450 5.000") # U=8
#betaarray=$( echo "2.000 2.100 2.200 2.250 2.300 2.350 2.400 2.450 2.500 2.550 2.600 2.650 2.700 2.800 2.850 2.950 3.100 3.350 3.650 4.000 4.450 5.000") # U=8
#betaarray=$( echo "2.250 2.300 2.400 2.450 2.550 2.600 2.700 2.800 2.950") # U=8
#betaarray=$( echo "3.350 3.650 4.000 4.450 5.000 5.700 6.650 8.000 10.000") # U=4
#betaarray=$( echo "3.450 3.550 3.750 3.850 4.150 4.300 4.650 4.800 5.300") # U=4
#betaarray=$( echo "2.500 2.650 2.850 3.100 3.350 3.650 4.000 4.450 5.000")  # U=12
#betaarray=$( echo "2.550 2.600 2.700 2.800 2.900 3.000 3.200 3.300 3.450 3.550 3.800")  # U=12
#betaarray=$( echo "2.500 2.550 2.600 2.650 2.700 2.800 2.850 2.900 3.000 3.100 3.200 3.300 3.350 3.450 3.550 3.650 3.700 3.800 3.900 4.000 4.450 5.000")  # U=12
#betaarray=$( echo "2.850 2.900 3.000 3.100 3.200 3.300 3.350 3.450 3.550 3.650 3.700 3.800 3.900")  # U=12
#betaarray=$( echo "3.300 3.350 3.450 3.550 3.650 3.700 3.800 3.900")  # U=12
#betaarray=$( echo "3.350 3.450 3.550 3.650 3.750 3.850 4.000 4.150 4.300 4.450 4.650 4.800 4.900 5.000 5.150 5.300 5.500 5.700 6.650 8.000 10.000") # U=4
#b
#betaarray=$( echo "4.300 4.450 4.650 4.800 4.900 5.000 5.150 5.300 5.500 5.700") # U=4
betaarray=$( echo "5.000") # U=4
#betaarray=$( echo "4.800 4.900 5.150 5.300") # U=4
#betaarray=$( echo "2.500 2.550 2.600 2.650 2.700 2.800 2.850 2.950 3.100 3.350") # U=8
#betaarray=$( echo "2.600 2.650 2.700 2.800 2.850 2.950 3.100 3.350 3.650 4.000 4.450 5.000") # U=6
#betaarray=$( echo "2.600 2.650 2.700 2.800 2.850 2.950 3.100 3.350 3.650 4.000 4.450 5.000") # U=10
#betaarray=$( echo "2.650 2.700 2.800 2.850 2.950 3.100 3.350") # U=8 L=16
#betaarray=$( echo "0.2 0.25 0.5 1.000 2.000 4.000 5.000 10.000")
#Larray=$( echo "18 21")
Larray=$( echo "4")
#Larray=$( echo "16")
#Larray=$( echo "4 6 8 10 12")
#Larray=$( echo "4 6 8 10")
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
nwrap=20
nsweep=10
nbin=30
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
