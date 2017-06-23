#!/bin/bash
#
# Run temper/grem with 16 replicas on 16 processors

export MPI_IMPLEMENTATION=openmpi

NSLOTS=16
NREPLICA=16 # Number of replicas

# Restart label
indx=$(cat last)

# Calculate partitions
NP_REP=$((${NSLOTS}/${NREPLICA}))
NPROCS=$((${NP_REP}*${NREPLICA}))
PROC=${NREPLICA}x${NP_REP} # M partitions, N procs per partition

# Run Simulation
if [ $indx -eq 0 ]; then
  mpirun -np $NPROCS ~/lmp_stable-latest -p $PROC -l log/log.lammps -pscreen none -in start.gREM-temper
else
  mpirun -np $NPROCS ~/lmp_stable-latest -p $PROC -l log/log.lammps -pscreen none -in temper.gREM 
fi
sleep 3

# Check sim, if finished advance counter
if [[ -s log/log.lammps ]]; then
  loglines=$(wc -l log/log.lammps | awk '{print $1}')
  if [ $loglines -lt 4 ]; then
    echo "Simulation failed! log/log.lammps has too few lines"
    exit 1
  else
    indx=$((indx+1))
    echo $indx > last
  fi
else
  echo "Simulation failed! log/log.lammps file empty"
  exit 1
fi

# Get new walker list
tail -1 log/log.lammps | awk '{for (i=2; i<=NF; i++) print $i}' | tr '\n' ' ' > lastwalker

# File 
CHK_WALKER=$(cat lastwalker | tr ' ' '\n' | wc -l)
if [ $CHK_WALKER -eq $(($NREPLICA)) ]; then
  cp log/log.lammps log/log.lammps-${indx}
  sed -e "s/ZZZ/$(sed 's:/:\\/:g' lastwalker)/" TEMPLATE.gREM-temper > temper.gREM
  for ((i=0; i<$NREPLICA; i++)); do
    cp log/log.lammps.$i $i/log.lammps.$i-${indx}
    cp $i/final_restart_file $i/final_restart_file-${indx}
    cp $i/final_restart_file $i/restart_file
  done
else
  echo "Restart failed! incorrect number of walker index"
  exit 1
fi
sleep 2

# uncomment to resubmit!
if [ "$indx" -lt 5 ]; then
  sleep 1
  ./submit-temper.bash
fi

exit 0
