#! /bin/bash 
num_iter=5
m_values=(10000 20000 30000 40000)

for ((iter=0; iter< $num_iter; iter++));
do
	for m in "${m_values[@]}"
	do

	echo "mpirun -n 256 --npernode 8 ./sys_3d_coded_summa -n 8 -k 8 -N 4 -K 2 -f 1 -m ${m} -l"
	mpirun -n 256 --npernode 8 ./sys_3d_coded_summa -n 8 -k 8 -N 4 -K 2 -f 1 -m ${m} -l       

	echo "mpirun -n 256 --npernode 8 ./replication_3d_summa -n 8 -k 8 -N 4 -K 2 -f 1 -m ${m}"
	mpirun -n 256 --npernode 8  ./replication_3d_summa -n 8 -k 8 -N 4 -K 2 -f 1 -m ${m}       

	echo "mpirun -n 256 --npernode 8 --hostfile ./hostfile ./3d_coded_summa -n 8 -k 8 -N 4 -K 2 -m ${m} -l"
	mpirun -n 256 --npernode 8  ./3d_coded_summa -n 8 -k 8 -N 4 -K 2  -m ${m} -l       

	done
done
