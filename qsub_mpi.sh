#!/bin/sh
#PBS -N qsub_mpi
#PBS -e test.e
#PBS -o test.o
#PBS -l nodes=1:ppn=4

NODES=$(cat $PBS_NODEFILE | sort | uniq)
# 注意把所有的ntt换成你的选题

for node in $NODES; do
    scp master_ubss1:/home/${USER}/ntt/main ${node}:/home/${USER} 1>&2
    scp -r master_ubss1:/home/${USER}/ntt/files ${node}:/home/${USER}/ 1>&2
done

# 使用tee命令同时输出到终端和文件
/usr/local/bin/mpiexec -np 2 -machinefile $PBS_NODEFILE /home/${USER}/main 2>&1 | tee -a test.o

scp -r /home/${USER}/files/ master_ubss1:/home/${USER}/ntt/ 2>&1
