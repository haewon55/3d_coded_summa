#! /bin/bash 
num_iter=1
m_values=(10000 20000 30000 40000)

for ((iter=0; iter< $num_iter; iter++));
do
	for m in "${m_values[@]}"
	do

	echo "mpirun -n 256 --npernode 8 --hostfile ./hostfile ./coded_25d_summa_sys -n 8 -k 8 -N 4 -K 2 -f 1 -m ${m}"
	mpirun -n 256 --npernode 8 --hostfile ./hostfile ./coded_25d_summa_sys -n 8 -k 8 -N 4 -K 2 -f 1 -m ${m}       

	#echo "mpirun -n 256 --npernode 8 --hostfile ./hostfile ./coded_25d_summa_sys -n 8 -k 8 -N 4 -K 2 -f 3 -m ${m}"
	#mpirun -n 256 --npernode 16 --hostfile ./hostfile ./coded_25d_summa_sys -n 8 -k 8 -N 4 -K 2 -f 3 -m ${m}       

	echo "mpirun -n 256 --npernode 8 --hostfile ./hostfile ./replica_25d_summa -n 8 -k 8 -N 4 -K 2 -f 1 -m ${m}"
	mpirun -n 256 --npernode 8 --hostfile ./hostfile ./replica_25d_summa -n 8 -k 8 -N 4 -K 2 -f 1 -m ${m}       

	#echo "mpirun -n 256 --npernode 8 --hostfile ./hostfile ./replica_25d_summa -n 8 -k 8 -N 4 -K 2 -f 3 -m ${m}"
	#mpirun -n 256 --npernode 16 --hostfile ./hostfile ./replica_25d_summa -n 8 -k 8 -N 4 -K 2 -f 3 -m ${m}       

	echo "mpirun -n 256 --npernode 8 --hostfile ./hostfile ./coded_25d_summa_mem -n 8 -k 8 -N 4 -K 2 -m ${m}"
	mpirun -n 256 --npernode 8 --hostfile ./hostfile ./coded_25d_summa_mem -n 8 -k 8 -N 4 -K 2  -m ${m}       


#	echo "mpirun -n 640 --npernode 16  ./coded_25d_summa_mem -n 8 -k 8 -m ${m} -N 10 -K 5"
#        mpirun -n 640 --npernode 16 --hostfile ./hostfile ./coded_25d_summa_mem -n 8 -k 8  -N 10 -K 5 -m ${m}
	done
done
