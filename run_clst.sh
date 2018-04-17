#!/bin/bash
usage() {
cat <<EOF
run.sh
---------

DESCRIPTION:

Run the neolith code on <NPROCS> processes.

USAGE:

        run.sh [OPTIONS] <NPROCS>

OPTIONS:        

        -c|--compile
                    Recompile the code - useful if cparam.local been changed.

        -q|--quiet
                    Use nohup to run the code in the background. This means 
                    the job will continue to run even if the terminal is closed.
                    All output is directed to a log file 'out.log'.

        -f|--force
                    If the script detects the existence of process directories
                    in the current directory, new process directories are not
                    created and the script aborts.  Use -f to force removal of
                    the existing process directories and create new ones.

        -h|--help
                    Show usage information.
EOF
}
run() {
  echo Running job...
  pres_dir=`pwd`
  echo " cd '$pres_dir'" > ./temp.sh
  echo "/opt/mpiexec/bin/mpiexec -np $NPROCS $EXE" >> ./temp.sh
  qsub -l nodes=$NPROCSA:ppn=$NPROCSB,walltime=00:10:00,cput=00:10:00 ./temp.sh
  #qsub -l nodes=$NPROCSA:ppn=$NPROCSB ./temp.sh
}
run_quiet() {
  echo Running job quietly...
  nohup /usr/bin/time -p openmpiexec -n $NPROCS $EXE &> $LOGFILE &
}
COMPILE=0
QUIET=0
FORCE=0
RESTART=0
while getopts ":hcn:qfrp:" option;
do
 case $option in
  h)
   usage
   exit 1
   ;;
  f)
   FORCE=1
   ;;
  r)
   RESTART=1
   ;;
  c)
   COMPILE=1
   ;;
  q)
   QUIET=1
   ;;
  n)
   NPROCSA=$OPTARG
   ;;
  p)
   NPROCSB=$OPTARG
   ;;
  :)
   echo "option -$OPTARG needs an argument"
   ;;
  *)
   echo "invalid option -$OPTARG" 
   ;;
 esac
done
DIGITS="2"
EXE="./src/run.x"
LOGFILE="out.log"

echo Going to run on $NPROCSA nodes. With $NPROCSB cores
NPROCS=$(($NPROCSA * $NPROCSB))
echo Going to run on $NPROCS processes

for i in `seq -f "%0${DIGITS}g" 0 $(($NPROCS-1))`
do
  PROCDIR="./data/proc$i"
  if [ $RESTART -eq 1 ]; then
    if [ $i==0 ] ; then
      echo Attempting to restart not checking directory structure
    fi
  else
    if [ -d $PROCDIR ]; then
      if [ $FORCE -eq 0 ]; then
        echo WARNING: At least one process directory already exists - aborting.
        echo Use -f to force deletion of process directories.
        exit 1
      else
        rm -rf $PROCDIR
      fi
    fi
    mkdir $PROCDIR
  fi
done
if [ $COMPILE -eq 1 ]; then
  # Recompile
  echo Recompiling code...
  make clean ; make
fi
if [ $QUIET -eq 1 ]; then
  # Run quietly
  run_quiet
  exit 0
fi

# Run the job.
run

