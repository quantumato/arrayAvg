echo "starting arrayAvg4"

for i in {1..10}
do
mpirun -np 64 arrayAvg4
done

echo "starting arrayAvg5"

for i in {1..10}
do
mpirun -np 64 arrayAvg5
done

