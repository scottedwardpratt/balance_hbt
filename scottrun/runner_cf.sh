#! /bin/bash
make balhbt
NPROC=24
BW_TAU=13.4
BW_RPERP=13.4
for ((i=0;i<${NPROC};i++))
do
	rm -f logfiles/balhbt${i}.txt;
	./balhbt ${BW_TAU} ${BW_RPERP} ${i} > logfiles/balhbt${i}.txt &
done
wait

./summer_cf
wait