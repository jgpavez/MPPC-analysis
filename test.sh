#!/bin/bash
a=0;
while read header line; do
	a=`expr $a + 1`
	if [ $a -eq 1 ]; then 
		box=$(printf %03d $line)
	fi
	if [ $a -eq 2 ]; then 
		directory=$line
	fi
	if [ $a -eq 3 ]; then 
		run=$line
	fi
		if [ $a -eq 4 ]; then 
		temp=$(printf %02d $line)
	fi
	if [ $a -eq 5 ]; then 
		maindir=$line
	fi
	if [ $a -eq 6 ]; then 
		vopth=$line
	fi
	if [ $a -eq 7 ]; then 
		begin=$line
	fi
	if [ $a -eq 8 ]; then 
		end=$line
	fi
done < "parameters.txt"

for i in $(seq 0 31); do 
	grep nan $maindir/USM$box/T$temp"_"$directory"_"run$run/results_MPPC_$i"_"* &>/dev/null;
	if [ $? -eq 0 ]; then 
		echo "mppc $i bad" 
	fi 
done

