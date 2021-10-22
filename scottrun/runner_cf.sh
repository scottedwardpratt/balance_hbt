#! /bin/bash
make balhbt
NMC=10000000
NPROC=12
BW_TAU=15.0
BW_R=12.0
for ((i=0;i<${NPROC};i++))
do
	rm -f logfiles/balhbt${i}.txt;
	./balhbt ${BW_TAU} ${BW_R} ${i} > logfiles/balhbt${i}.txt &
done
