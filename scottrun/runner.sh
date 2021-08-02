#! /bin/bash
make balhbt
rm logfiles/*.txt
for ((i=0;i<24;i++))
do
	rm -f logfiles/balhbt${i}.txt;
	echo 1000000000 | ./balhbt ${i} > logfiles/balhbt${i}.txt &
done
