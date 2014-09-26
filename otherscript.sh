#!/bin/bash

# This is the base script I always will send to Hoffman. It sets up the MFILE
# necessary for certain parallelization in some TURBOMOLE scripts and then
# keeps track of where it is in the "script.sh" which contains the actual
# TURBOMOLE commands. Ttimer is a subroutine that activates when time is almost
# out. Once activated, Ttimer will send a STOP command to TURBOMOLE and submit
# the remainingscript.sh to start where it left off. One current limitation:
# all TURBOMOLE commands must be in ONE line, so if there are any if/for/while
# statements... they will have to be in one line. This is all very against
# standard style however not sure of another way to get functionality I want.

genMFILE
Ttimer $$

scriptLC=`wc -l < script.sh`

counter=1
while [ $counter -le $scriptLC ]; do
 cat script.sh | tail -$(($scriptLC-$counter+1)) > _remainingscript.sh
 mv _remainingscript.sh remainingscript.sh
 if [ -f runthis.sh ]; then rm runthis.sh; fi
 cat remainingscript.sh | head -n +1 > runthis.sh
 chmod +x runthis.sh
 if [ ! -f "$SGE_CWD_PATH/pleaseResubmit" ]; then
  runthis.sh
  if [ ! -f "$SGE_CWD_PATH/pleaseResubmit" ]; then
   counter=$(($counter+1))
  else
   break
  fi
 else
  break
 fi
done

if [ $counter -gt $scriptLC ]; then
 echo "ran every line of script.sh"
else
 echo "didn't get to every line of script.sh, will resubmit"
 sleep 500h
fi
