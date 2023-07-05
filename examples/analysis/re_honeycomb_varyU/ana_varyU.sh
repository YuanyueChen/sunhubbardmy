#!/bin/bash

source cal_para.sh
nskip=1
nskps=0
nskpl=0

WORKDIR="$PWD"
datadir=$WORKDIR/../../run_honeycomb/
echo $WORKDIR
cd $WORKDIR
for beta in ${betaarray}; do
for plqu in ${plquarray}; do
for nu in ${nuarray}; do
for nflr in ${nflrarray}; do
  for L  in ${Larray}; do
       pretag=${code}.b${beta}.dtau${dtau}.t${rt}.plqU${plqu}.alpha${alpha}.theta${theta}.Nf${nflr}.nu${nu}.L${L}.proj${lprojplqu}${lproju}
       cd $WORKDIR
       rm ${pretag}*.dat
       for u in ${uarray}; do
         cd $datadir
         maindir=${code}.b${beta}.dtau${dtau}.t${rt}.U${u}.plqU${plqu}.alpha${alpha}.theta${theta}.Nf${nflr}.nu${nu}.L${L}.proj${lprojplqu}${lproju}
         if [ -f $maindir/energy.bin ]; then
            cd $maindir
	        lq=$( echo "$L" |awk '{print $1^2}' )
            echo "processing $maindir ... "
	        echo " lq = $lq "
            # analysis energy.bin

            # sgn
            echo "    get sgn ... "
            #nskip=$( wc -l energy.bin |awk '{print $1/4}' ) # set nskip to be nbin/4
            awk  '{if(NR>nskipv) print $1}' nskipv=$nskip energy.bin |sort -n > ener.tmp.tmp
            nb=$( wc -l ener.tmp.tmp |awk '{print $1}' )
            awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' ener.tmp.tmp > ener.tmp
            awk 'function abs(v) {return v < 0 ? -v : v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%24.12f %24.12f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
            ener.tmp > ener.tmp2
            awk '{if(NR==1) print uv, $0}' uv=$u ener.tmp2 >> $WORKDIR/${pretag}_sgn.dat

            # den
            echo "    get den ... "
            awk  '{if(NR>nskipv) print $3}' nskipv=$nskip energy.bin |sort -n > ener.tmp.tmp
            nb=$( wc -l ener.tmp.tmp |awk '{print $1}' )
            awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' ener.tmp.tmp > ener.tmp
            awk 'function abs(v) {return v < 0 ? -v : v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%24.12f %24.12f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
            ener.tmp > ener.tmp2
            awk '{if(NR==1) print uv, $1, $2}' uv=$u Lv=$L ener.tmp2 >> $WORKDIR/${pretag}_den.dat

            # ekint
            echo "    get ekint ... "
            awk  '{if(NR>nskipv) print $5}' nskipv=$nskip energy.bin |sort -n > ener.tmp.tmp
            nb=$( wc -l ener.tmp.tmp |awk '{print $1}' )
            awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' ener.tmp.tmp > ener.tmp
            awk 'function abs(v) {return v < 0 ? -v : v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%24.12f %24.12f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
            ener.tmp > ener.tmp2
            awk '{if(NR==1) print uv, $1/lqv, $2/lqv}' uv=$u lqv=$lq ener.tmp2 >> $WORKDIR/${pretag}_ekint.dat

            # eu
            echo "    get eu ... "
            awk  '{if(NR>nskipv) print $7}' nskipv=$nskip energy.bin |sort -n > ener.tmp.tmp
            nb=$( wc -l ener.tmp.tmp |awk '{print $1}' )
            awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' ener.tmp.tmp > ener.tmp
            awk 'function abs(v) {return v < 0 ? -v : v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%24.12f %24.12f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
            ener.tmp > ener.tmp2
            awk '{if(NR==1) print uv, $1/lqv, $2/lqv}' uv=$u lqv=$lq ener.tmp2 >> $WORKDIR/${pretag}_eu.dat

            rm ener.tmp2
            rm ener.tmp
            rm ener.tmp.tmp

            if [ -f spsm_orb1_k.bin ]; then
                # analysis spin correlation, from spsm_orb1_k.bin 
                # also calculate correlation ratio
                echo "    get spsm(pi,pi) ... "
                #nskip=$( grep -e "  0.00000000E+00  0.00000000E+00" spsm_orb1_k.bin |wc -l |awk '{print $1/6}' )
                grep -e "  0.00000000E+00  0.00000000E+00" spsm_orb1_k.bin |awk '{if(NR>nskipv) print $3-$5-$7+$9}' nskipv=$nskip |sort -n > spsm_pi.tmp.tmp
                nb=$( wc -l spsm_pi.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' spsm_pi.tmp.tmp > spsm_pi.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
                spsm_pi.tmp > spsm_pi.tmp2
                awk '{if(NR==1) print uv, $0}' uv=$u spsm_pi.tmp2 >> $WORKDIR/${pretag}_spsm_pi.dat

                grep -A 1 "  0.00000000E+00  0.00000000E+00" spsm_orb1_k.bin |grep -v "  0.00000000E+00  0.00000000E+00" |awk '{if(NF>3) print $3-$5-$7+$9}' \
                          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
                nb=$( wc -l spsm_pidq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' spsm_pidq.tmp.tmp > spsm_pidq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
                spsm_pidq.tmp > spsm_pidq.tmp2

                paste spsm_pi.tmp2 spsm_pidq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print uv, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print uv, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' uv=$u >> $WORKDIR/${pretag}_spsm_pi_r10.dat

                la1=$( echo "$L" |awk '{print $1+1}' )
                la2=$( echo "$L" |awk '{print $1+2}' )
                grep -A $la1 "  0.00000000E+00  0.00000000E+00" spsm_orb1_k.bin |grep -v "  0.00000000E+00  0.00000000E+00" |awk -v l1v=$la1 -v l2v=$la2 '{if(NR%l2v==l1v) print $3-$5-$7+$9}' \
                          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
                nb=$( wc -l spsm_pidq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' spsm_pidq.tmp.tmp > spsm_pidq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
                spsm_pidq.tmp > spsm_pidq.tmp2

                paste spsm_pi.tmp2 spsm_pidq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print uv, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print uv, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' uv=$u >> $WORKDIR/${pretag}_spsm_pi_r11.dat

                grep -A 2 "  0.00000000E+00  0.00000000E+00" spsm_orb1_k.bin |grep -v "  0.00000000E+00  0.00000000E+00" |awk -v l1v=$la1 -v l2v=$la2 '{if(NR%3==2) print $3-$5-$7+$9}' \
                          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
                nb=$( wc -l spsm_pidq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' spsm_pidq.tmp.tmp > spsm_pidq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
                spsm_pidq.tmp > spsm_pidq.tmp2

                paste spsm_pi.tmp2 spsm_pidq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print uv, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print uv, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' uv=$u >> $WORKDIR/${pretag}_spsm_pi_r20.dat

                rm spsm_pi.tmp2
                rm spsm_pi.tmp
                rm spsm_pidq.tmp2
                rm spsm_pidq.tmp
                rm spsm_pi.tmp.tmp spsm_pidq.tmp.tmp
            fi

            if [ -f bb_orb1_k.bin ]; then
                # analysis spin correlation, from bb_orb1_k.bin 
                # also calculate correlation ratio
                echo "    get bb(K) ... "
                #nskip=$( grep -e " -0.33333333E+00 -0.33333333E+00" bb_orb1_k.bin |wc -l |awk '{print $1/6}' )
                grep -e " -0.33333333E+00 -0.33333333E+00" bb_orb1_k.bin |awk '{if(NR>nskipv) print $3}' nskipv=$nskip |sort -n > bb_K.tmp.tmp
                nb=$( wc -l bb_K.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' bb_K.tmp.tmp > bb_K.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
                bb_K.tmp > bb_K.tmp2
                awk '{if(NR==1) print uv, $0}' uv=$u bb_K.tmp2 >> $WORKDIR/${pretag}_bb_K.dat

                grep -A 1 " -0.33333333E+00 -0.33333333E+00" bb_orb1_k.bin |grep -v " -0.33333333E+00 -0.33333333E+00" |awk '{if(NF>3) print $3}' \
                          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > bb_Kdq.tmp.tmp
                nb=$( wc -l bb_Kdq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' bb_Kdq.tmp.tmp > bb_Kdq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
                bb_Kdq.tmp > bb_Kdq.tmp2

                paste bb_K.tmp2 bb_Kdq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print uv, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print uv, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' uv=$u >> $WORKDIR/${pretag}_bb_K_r10.dat

                la1=$( echo "$L" |awk '{print $1+1}' )
                la2=$( echo "$L" |awk '{print $1+2}' )
                grep -A $la1 " -0.33333333E+00 -0.33333333E+00" bb_orb1_k.bin |grep -v " -0.33333333E+00 -0.33333333E+00" |awk -v l1v=$la1 -v l2v=$la2 '{if(NR%l2v==l1v) print $3}' \
                          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > bb_Kdq.tmp.tmp
                nb=$( wc -l bb_Kdq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' bb_Kdq.tmp.tmp > bb_Kdq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
                bb_Kdq.tmp > bb_Kdq.tmp2

                paste bb_K.tmp2 bb_Kdq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print uv, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print uv, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' uv=$u >> $WORKDIR/${pretag}_bb_K_r11.dat

                grep -A 2 " -0.33333333E+00 -0.33333333E+00" bb_orb1_k.bin |grep -v " -0.33333333E+00 -0.33333333E+00" |awk -v l1v=$la1 -v l2v=$la2 '{if(NR%3==2) print $3}' \
                          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > bb_Kdq.tmp.tmp
                nb=$( wc -l bb_Kdq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' bb_Kdq.tmp.tmp > bb_Kdq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR/NR) )} }'\
                bb_Kdq.tmp > bb_Kdq.tmp2

                paste bb_K.tmp2 bb_Kdq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print uv, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print uv, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' uv=$u >> $WORKDIR/${pretag}_bb_K_r20.dat

                rm bb_K.tmp2
                rm bb_K.tmp
                rm bb_Kdq.tmp2
                rm bb_Kdq.tmp
                rm bb_K.tmp.tmp bb_Kdq.tmp.tmp
            fi
         fi
       done
       paste $WORKDIR/${pretag}_ekint.dat $WORKDIR/${pretag}_eu.dat |awk -v uv=$u '{print $1, $2*2+$5*uv, sqrt(4*$3*$3+$6*$6*uv*uv)}' > $WORKDIR/${pretag}_etot.dat
  done
done
done
done
done
