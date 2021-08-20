#! /bin/bash
make balhbt
NMC=1000000000
NPROC=12
BW_TAU=$1
BW_RPERP=$2
for ((i=0;i<${NPROC};i++))
do
	rm -f logfiles/balhbt $1 $2 ${i}.txt;
	./balhbt $1 $2 ${i} &
done
wait

./summer.exe
wait