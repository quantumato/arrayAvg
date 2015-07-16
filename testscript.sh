#!/bin/bash
echo "starting arrayAvg4"

for i in 4, 9, 16, 25, 64
do
	for j in 10, 50, 100, 1000
	do
		for k in {1..3}
		do
			mpirun -np $i arrayAvg4 $j $j
		done
	done
done

echo "starting arrayAvg5"

for i in {1..10}
do
mpirun -np 64 arrayAvg5
done

echo "starting arrayAvg6"
