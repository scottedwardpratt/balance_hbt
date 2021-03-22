#! /bin/bash
make balance_hbt
for ((i=0;i<12;i++))
do
	rm -f logfiles/balhbt${i}.txt;
	./balance_hbt ${i} > logfiles/balhbt${i}.txt &
done
