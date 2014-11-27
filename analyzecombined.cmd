#!/bin/csh -f
#  analyzecombined.cmd
#
#  UGE job for analyzecombined built Tue Nov  4 01:35:23 PST 2014
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/home/m/mnechay/mikescripts/analyzecombined.joblog.$JOB_ID
#$ -o /u/home/m/mnechay/mikescripts/analyzecombined.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/home/m/mnechay/mikescripts/analyzecombined
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#$ -pe shared 4
#$ -l h_data=1024M,h_rt=12:00:00
#
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M mnechay@mail
#  Notify at beginning and end of job
#$ -m bea
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "job serial"
  set qqidir    = /u/home/m/mnechay/mikescripts
  set qqjob     = analyzecombined
  set qqodir    = /u/home/m/mnechay/mikescripts
  cd     /u/home/m/mnechay/mikescripts
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for analyzecombined built Tue Nov  4 01:35:23 PST 2014"
  echo ""
  echo "  analyzecombined directory:"
  echo "    "/u/home/m/mnechay/mikescripts
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "analyzecombined started on:   "` hostname -s `
  echo "analyzecombined started at:   "` date `
  echo ""
#
# Run the user program
#
  source /u/local/Modules/default/init/modules.csh
  module load intel/13.cs
#
  echo analyzecombined "" \>\& analyzecombined.output.$JOB_ID
  echo ""
  time /u/home/m/mnechay/mikescripts/analyzecombined  >& /u/home/m/mnechay/mikescripts/analyzecombined.output.$JOB_ID
#
  echo ""
  echo "analyzecombined finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/home/m/mnechay/mikescripts/analyzecombined.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/job.log.serial
  if (`wc -l /u/home/m/mnechay/mikescripts/analyzecombined.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
	head -50 /u/home/m/mnechay/mikescripts/analyzecombined.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
	echo " "  >> /u/local/apps/queue.logs/job.log.serial
	tail -10 /u/home/m/mnechay/mikescripts/analyzecombined.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  else
	cat /u/home/m/mnechay/mikescripts/analyzecombined.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  endif
  exit (0)
