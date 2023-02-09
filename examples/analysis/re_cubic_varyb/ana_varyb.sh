#!/bin/bash

source cal_para.sh
nskip=2
nskps=0
nskpl=0

WORKDIR="$PWD"
datadir=$WORKDIR/../../run_cubic/
echo $WORKDIR
cd $WORKDIR
for u in ${uarray}; do
for plqu in ${plquarray}; do
for nu in ${nuarray}; do
for nflr in ${nflrarray}; do
  for L  in ${Larray}; do
       pretag=${code}.dtau${dtau}.t${rt}.U${u}.plqU${plqu}.alpha${alpha}.theta${theta}.Nf${nflr}.nu${nu}.L${L}.proj${lprojplqu}${lproju}
       cd $WORKDIR
       rm ${pretag}*.dat
       for beta in ${betaarray}; do
         cd $datadir
         maindir=${code}.b${beta}.dtau${dtau}.t${rt}.U${u}.plqU${plqu}.alpha${alpha}.theta${theta}.Nf${nflr}.nu${nu}.L${L}.proj${lprojplqu}${lproju}
         if [ -f $maindir/energy.bin ]; then
            cd $maindir
	    lq=$( echo "$L" |awk '{print $1^3}' )
            echo "processing $maindir ... "
	    echo " lq = $lq "
            # analysis energy.bin

            # sgn
            echo "    get sgn ... "
            #nskip=$( wc -l energy.bin |awk '{print $1/4}' ) # set nskip to be nbin/4
            awk  '{if(NR>nskipv) print $1}' nskipv=$nskip energy.bin |sort -n > ener.tmp.tmp
            nb=$( wc -l ener.tmp.tmp |awk '{print $1}' )
            awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' ener.tmp.tmp > ener.tmp
            awk 'function abs(v) {return v < 0 ? -v : v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%24.12f %24.12f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
            ener.tmp > ener.tmp2
            awk '{if(NR==1) print betav, $0}' betav=$beta ener.tmp2 >> $WORKDIR/${pretag}_sgn.dat

            # den
            echo "    get den ... "
            awk  '{if(NR>nskipv) print $3}' nskipv=$nskip energy.bin |sort -n > ener.tmp.tmp
            nb=$( wc -l ener.tmp.tmp |awk '{print $1}' )
            awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' ener.tmp.tmp > ener.tmp
            awk 'function abs(v) {return v < 0 ? -v : v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%24.12f %24.12f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
            ener.tmp > ener.tmp2
            awk '{if(NR==1) print betav, $1, $2}' betav=$beta Lv=$L ener.tmp2 >> $WORKDIR/${pretag}_den.dat

            # ekint
            echo "    get ekint ... "
            awk  '{if(NR>nskipv) print $5}' nskipv=$nskip energy.bin |sort -n > ener.tmp.tmp
            nb=$( wc -l ener.tmp.tmp |awk '{print $1}' )
            awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' ener.tmp.tmp > ener.tmp
            awk 'function abs(v) {return v < 0 ? -v : v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%24.12f %24.12f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
            ener.tmp > ener.tmp2
            awk '{if(NR==1) print betav, $1/lqv, $2/lqv}' betav=$beta lqv=$lq ener.tmp2 >> $WORKDIR/${pretag}_ekint.dat

            # eu
            echo "    get eu ... "
            awk  '{if(NR>nskipv) print $7}' nskipv=$nskip energy.bin |sort -n > ener.tmp.tmp
            nb=$( wc -l ener.tmp.tmp |awk '{print $1}' )
            awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' ener.tmp.tmp > ener.tmp
            awk 'function abs(v) {return v < 0 ? -v : v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%24.12f %24.12f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
            ener.tmp > ener.tmp2
            awk '{if(NR==1) print betav, $1/lqv, $2/lqv}' betav=$beta lqv=$lq ener.tmp2 >> $WORKDIR/${pretag}_eu.dat

            rm ener.tmp2
            rm ener.tmp
            rm ener.tmp.tmp

            if [ -f spsm_orb1_k.bin ]; then
                # analysis spin correlation, from spsm_orb1_k.bin 
                # also calculate correlation ratio
	        lq=$( echo "$L" |awk '{print $1^3}' )
                echo "    get spsm(pi,pi,pi) ... "
                #nskip=$( grep -e "  0.50000000E+00  0.50000000E+00" spsm_orb1_k.bin |wc -l |awk '{print $1/6}' )
		awk -v lqv=$lq '{if(NR%lqv==0) print $0}' spsm_orb1_k.bin |awk '{if(NR>nskipv) print $4}' nskipv=$nskip |sort -n > spsm_pi.tmp.tmp
                nb=$( wc -l spsm_pi.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' spsm_pi.tmp.tmp > spsm_pi.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
                spsm_pi.tmp > spsm_pi.tmp2
                awk '{if(NR==1) print betav, $0}' betav=$beta spsm_pi.tmp2 >> $WORKDIR/${pretag}_spsm_pi.dat

                #grep -B 1 "  0.50000000E+00  0.50000000E+00  0.50000000E+00" spsm_orb1_k.bin |grep -v "  0.50000000E+00  0.50000000E+00  0.50000000E+00" |awk '{if(NF>3) print $4}' \
                #          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
		awk -v lqv=$lq '{if(NR%lqv==lqv-1) print $0}' spsm_orb1_k.bin |awk '{if(NR>nskipv) print $4}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
                nb=$( wc -l spsm_pidq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' spsm_pidq.tmp.tmp > spsm_pidq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
                spsm_pidq.tmp > spsm_pidq.tmp2

                paste spsm_pi.tmp2 spsm_pidq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print betav, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print betav, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' betav=$beta >> $WORKDIR/${pretag}_spsm_pi_r100.dat

                la1=$( echo "$L" |awk '{print $1+1}' )
                la2=$( echo "$L" |awk '{print $1+2}' )
                #grep -B $la1 "  0.50000000E+00  0.50000000E+00  0.50000000E+00" spsm_orb1_k.bin |grep -v "  0.50000000E+00  0.50000000E+00  0.50000000E+00" |awk -v l1v=$la1 -v l2v=$la2 '{if(NR%l2v==1) print $4}' \
                #          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
		awk -v lqv=$lq -v la1v=$la1 '{if(NR%lqv==lqv-la1v) print $0}' spsm_orb1_k.bin |awk '{if(NR>nskipv) print $4}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
                nb=$( wc -l spsm_pidq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' spsm_pidq.tmp.tmp > spsm_pidq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
                spsm_pidq.tmp > spsm_pidq.tmp2

                paste spsm_pi.tmp2 spsm_pidq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print betav, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print betav, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' betav=$beta >> $WORKDIR/${pretag}_spsm_pi_r110.dat

                la1=$( echo "$L" |awk '{print $1*$1+1}' )
                la2=$( echo "$L" |awk '{print $1*$1+2}' )
                #grep -B $la1 "  0.50000000E+00  0.50000000E+00  0.50000000E+00" spsm_orb1_k.bin |grep -v "  0.50000000E+00  0.50000000E+00  0.50000000E+00" |awk -v l1v=$la1 -v l2v=$la2 '{if(NR%l2v==1) print $4}' \
                #          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
		awk -v lqv=$lq -v la1v=$la1 '{if(NR%lqv==lqv-la1v) print $0}' spsm_orb1_k.bin |awk '{if(NR>nskipv) print $4}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
                nb=$( wc -l spsm_pidq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' spsm_pidq.tmp.tmp > spsm_pidq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
                spsm_pidq.tmp > spsm_pidq.tmp2

                paste spsm_pi.tmp2 spsm_pidq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print betav, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print betav, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' betav=$beta >> $WORKDIR/${pretag}_spsm_pi_r111.dat

                grep -B 2 "  0.50000000E+00  0.50000000E+00  0.50000000E+00" spsm_orb1_k.bin |grep -v "  0.50000000E+00  0.50000000E+00  0.50000000E+00" |awk -v l1v=$la1 -v l2v=$la2 '{if(NR%3==1) print $4}' \
                          |awk '{if(NR>nskipv) print $0}' nskipv=$nskip |sort -n > spsm_pidq.tmp.tmp
                nb=$( wc -l spsm_pidq.tmp.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' spsm_pidq.tmp.tmp > spsm_pidq.tmp
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
                spsm_pidq.tmp > spsm_pidq.tmp2

                paste spsm_pi.tmp2 spsm_pidq.tmp2 |awk 'function abs(v) {return v<0?-v:v} \
                {if(NR==1 && abs($3)<abs($1)) print betav, 1-$3/$1, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1); \
                 else if(NR==1) print betav, 0, sqrt($4*$1*$4*$1+$2*$3*$2*$3)/($1*$1)}' betav=$beta >> $WORKDIR/${pretag}_spsm_pi_r200.dat

                rm spsm_pi.tmp2
                rm spsm_pi.tmp
                rm spsm_pidq.tmp2
                rm spsm_pidq.tmp
                rm spsm_pi.tmp.tmp spsm_pidq.tmp.tmp
            fi
            if [ -f spsm_ktau.bin ]; then
		lta1=$( echo "$beta $dtau" |awk '{print int(($1+0.001)/$2)+1}' )
		lqta1=$( echo "$lq $lta1" |awk '{print $1*$2}')
		echo "lta1 = $lta1"
		echo "lqta1 = $lqta1"
		awk -v lta1v=$lta1 -v lqta1v=$lqta1 -v nskipv=$nskip '{if(NR>lqta1v*nskipv) {if(NR%lqta1v>=lqta1v-lta1v+2 || NR%lqta1v==0) sum +=$1; if(NR%lqta1v==0) {printf( "%12.6f \n", sum/(lta1v-1)); sum=0}} }' spsm_ktau.bin |sort -n > chi.tmp
                nb=$( wc -l chi.tmp |awk '{print $1}' )
                awk -v nbv=$nb -v nsv=$nskps -v nlv=$nskpl '{if(NR>nsv && NR<=nbv-nlv) print $0}' chi.tmp > chi.tmp1
                awk 'function abs(v) {return v<0?-v:v} {for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} END {for (i=1;i<=NF;i++)  { printf( "%12.6f %12.6f \n", sum[i]/NR, sqrt(abs(sumsq[i]-sum[i]^2/NR)/NR) )} }'\
                chi.tmp1 > chi.tmp2
                awk '{if(NR==1) print betav, $0}' betav=$beta chi.tmp2 >> $WORKDIR/${pretag}_chi.dat
		rm chi.tmp chi.tmp1 chi.tmp2
            fi

         fi
       done
       paste $WORKDIR/${pretag}_ekint.dat $WORKDIR/${pretag}_eu.dat |awk -v uv=$u '{print $1, $2*2+$5*uv, sqrt(4*$3*$3+$6*$6*uv*uv)}' > $WORKDIR/${pretag}_etot.dat
  done
done
done
done
done
