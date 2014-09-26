#!/bin/bash 

for use in ana cediev87 djurwin enolate gerardoc ha22m lushen mnechay ngallup oirving pcbao qsmith root snedd xmyway zhangjin iz0mbie  
do
	run=$(qstat -u $use | grep " r " | wc -l)
	que=$(qstat -u $use | grep " qw" | wc -l)
	eue=$(qstat -u $use | grep "Eqw" | wc -l)
	echo " $use has running: $run, queued: $que, and err'd: $eue " 
done 

echo "                _________________________          " 
echo "              DON'T TOUCH, I WILL FIND YOU         " 
echo "                _________________________          " 

exit 0 

