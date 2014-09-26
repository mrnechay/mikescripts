#!/bin/bash
#
queue=" "
nofile=0
fileexist=0
toomanyargs=0
DIR=`pwd`
submit=0
cores="4"
node_type="4"
wall="24"
mem="0"
jobname="script.sh"
#
ECHO='echo -en'
#
usage() {
   $ECHO "\n\tAN ERROR WAS ENCOUNTERED:\n\n"
   $ECHO "\tUsage:\trunturbomole2.sh [-N cores] [-t time] [-q queue] [-sub] [-help] [-a jobname]  script_name \n\n"
   $ECHO "\t\t-N\t  :\tcores to use, \t\t\t\t[default:     8      \t]\n"
   $ECHO "\t\t-type\t :\tnode_typeto use, \t\t\t\t[default:     8      \t]\n"
   $ECHO "\t\t-mem\t  :\tmemory used,\t\t\t\t[default: 1024Mb*core\t] \n"
   $ECHO "\t\t-t\t  :\tWalltime to use, in hours.\t\t[default:    24   \t] \n"
   $ECHO "\t\t-sub\t  :\t(no argument) Submit the job.\t\t[default: don't submit\t]\n"
   $ECHO "\t\t-q\t  :\tSpecify the queue (highp or std).\t[default:    \t]\n"
   $ECHO "\t\t-a\t  :\tSpecify the job name.\t\t\t[default:    "script.sh"\t]\n"
   $ECHO "\t\t-help\t  :\tPrint this message.\n\n"
   if [ $fileexist -eq 1 ]; then
      $ECHO "\tERROR: The file you specified did not exist...\n\n"
   fi
   if [ $nofile -eq 1 ]; then
      $ECHO "\tERROR: You did not specify a file...\n\n"
   fi
   if [ $toomanyargs -eq 1 ]; then
      $ECHO "\tERROR: I expect only one argument...\n\n"
   fi
   exit 10
}
#
makefile() {

runfile=$1.tmole
name=$1


cat << %EOF% > $runfile
#!/bin/csh -f
#
#  SGE job for script built  $(date +%c)
#
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = $DIR/$1.joblog
#$ -o $DIR/$1.joblog.\$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  Parallelism:  $cores-way parallel

#  Resources requested
#$ -pe ${node_type}threads* $cores
##$ -pe dc_* $cores
#$ -l h_data=$((mem))G,h_rt=$wall:00:00${queue}

#  Name of application for log
#$ -v QQAPP=openmpi

#  Email address to notify
#$ -M gallup@chem.ucla.edu
#  Notify at beginning and end of job
#$ -m bea

#  Job name
#$ -N $jobname

#  Job is not rerunable
#$ -r n

#  Uncomment the next line to have your environment variables used by SGE
#$ -V

#
# Initialization for mpi parallel execution
#
  unalias *
  set qqversion = 
  set qqapp     = "openmpi parallel"
  set qqptasks  = $cores
  set qqidir    = $DIR
  set qqjob     = $1
  set qqodir    = $DIR

  ##set global_scratch = /u/scratch/$USER/job_\$JOB_ID
 set global_scratch = $SCRATCH/job_\$JOB_ID
 #  set global_scratch = $HOME/tmp/job_\$JOB_ID

  mkdir \$global_scratch
  echo "I created the scratch dir" \$global_scratch

  source /u/local/bin/qq.sge/qr.runtime
  if (\$status != 0) exit (1)
#
  echo "I sourced the runtime file"

  echo "  script directory: "\$qqidir
  echo "  Submitted to SGE: "\$qqsubmit
  echo "  'scratch' directory (on each node): "\$qqscratch
  echo "  script $cores-way parallel job configuration:"
  echo "    \$qqconfig" | tr "\\\\" "\n"
#
  echo ""
  echo "script started on:   "\` hostname -s \`
  echo "script started at:   "\` date \`
  echo ""
#
# Run the user program
#
   if !(\$?LD_LIBRARY_PATH) then
       setenv LD_LIBRARY_PATH  /u/local/intel/11.1/openmpi/1.4.2/lib:/u/local/compilers/intel/11.1/current/lib/intel64
   else
       setenv LD_LIBRARY_PATH /u/local/intel/11.1/openmpi/1.4.2/lib:/u/local/compilers/intel/11.1/current/lib/intel64:\${LD_LIBRARY_PATH}
   endif

   echo "    \$qqconfig" | tr "\\\\" "\n" | awk '{print \$2, \$4 }' > \$global_scratch/nodes.\$JOB_ID
   ~snedd/bin/list_hosts.sh \$global_scratch/nodes.\$JOB_ID \$global_scratch/nodes_long.\$JOB_ID

  setenv TURBODIR /u/home/ana/snedd/bin/TURBOMOLE
  source \$TURBODIR/Config_turbo_env.tcsh
  setenv PARA_ARCH SMP
  set path=(\$TURBODIR/bin/`sysname` \$path)
  setenv PARNODES  $cores
  setenv HOSTS_FILE \$global_scratch/nodes_long.\$JOB_ID

  setenv TURBOTMPDIR \$global_scratch

  cp -r \$qqidir/* \$global_scratch/
  cd \$global_scratch/
  rm $1.joblog.\$JOB_ID

  set time=(10 "System time: %S  User time: %U  Elapsed time: %E")
  time   $1  >& $1.output.\$JOB_ID

# clean up a little bit  
  rm slave*
  rm nodes*
  rm *dsk*

  cp -r \$global_scratch/* \$qqodir/
  rm -rf \$global_scratch/*

  echo ""
  echo "$1 finished at:  "` date `

#
# Cleanup after mpi parallel execution
  rm -rf \$global_scratch
#
  source /u/local/bin/qq.sge/qr.runtime
#
  exit (0)

%EOF%

}




#
while [ -n "`echo $1 | grep '^-'`" ]; do
   case $1 in
      -help)
         usage ;;
      -sub)
         submit=1 ;;
      -N)
         cores=$2
         shift ;;
      -a)           #Name addition
         jobname=$2
         shift ;;
      -type)
         node_type=$2
         shift ;;
      -q)
         queue=$2
         if [ "$queue" == "std" ]; then
            queue=" "
         fi
         if [ "$queue" == "highp" ]; then
            queue=",highp"
         fi
         shift ;;
      -mem)
         mem=$2
         shift ;;
      -t)
         wall=$2
         shift ;;
      -*)
         usage ;;
   esac
   shift
done
#
#echo "\$1 = $1"
if [ "$1" == "" ]; then
   nofile=1
   usage
fi
#
if [ ! -e "$1" ]; then
   fileexist=1
   usage
fi
#


if [ "$mem" == "0" ]; then
  mem=4
  (( mem*=$cores ))
  if [ "$mem" -gt "8" ]; then mem=4; fi 
#echo $mem
fi

if [ "$cores" -gt "1" ]; then
  export PARA_ARCH=SMP
fi

#
runfile="thisiscrap"
makefile $1 $runfile
chmod a+x $1
#
if [ $submit -eq 0 ]; then
   $ECHO "Now just do 'qsub $runfile'\n"
else
   $ECHO "Submitting '$runfile'\n"
   qsub $runfile
fi
