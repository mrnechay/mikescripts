#!/usr/bin/perl
#
# Script to remove duplicates of pairwise potentials from inConstr
# The original inConstr file is saved in inConstr.duplicate.sav
# The new inConstr file is saved in inConstr.temp
#
# by Sean Nedd
#

use strict;
use warnings;
#use diagnostics;
 
print "This script will remove duplicate pairwise potentials in inConstr in QMDMD.\n";
print "The original inConstr file is saved in inConstr.duplicate.sav.";
print "The new inConstr file is saved in inConstr.temp.";
#

my($a);
my($b);
my($c);
my($aa);
my($bb);
my($cc);
my($dup);
my($dupcounter);
my($INC_file); 
my($INC_SAV_file);
my($TEMP_file);
my(@pp);
my(@pp_label);
my(@non_pp);

$a=0;
$b=0;
$c=0;
$aa=0;
$bb=0;
$cc=0;
$dup=0;
$dupcounter=0;

$INC_file="inConstr";
$INC_SAV_file="inConstr.duplicate.sav";
$TEMP_file="inConstr.temp";

#Saving inConstr file
if ( -s $INC_file ) {
 print "Your inConstr file has been saved in inConstr.duplicate.sav\n";
 system("cat $INC_file>>$INC_SAV_file");
 #system("rm $INC_file");
}
open(INC, "<$INC_file") || die("Could not open INC_file!");
open(TMP, "+>$TEMP_file") || die("Could not open TEMP_file!");

#extracting pairwise-potentials, saved to @pp array
while (<INC>) {
 chomp($_);
 if ( /AtomPair / ) {
   push(@pp, $_);
  $a++;
 } else {
  push(@non_pp, $_);
  $b++;
 }
}

#removing duplicate; Compares a:x ... n-1 to b:x-1 ... n, saves a if no duplicate in b 
for ($aa = 0; $aa < @pp-1; $aa++) {
 for ($cc = $aa+1; $cc < @pp; $cc++) {
  if ( $pp[$aa] eq $pp[$cc] ) {
   print "\nDuplicate present for:\n$pp[$aa]";
   $dup = 1;
   $dupcounter++;
  }
 } 
 
 if ( $dup != 1) {
  print TMP "$pp[$aa]\n";
 }
 $dup=0;
}

#Adding the remainder of the inConstr file

for ($bb = 0; $bb <= $b-1; $bb++) {
 print TMP "$non_pp[$bb]\n";
}

system("cp $TEMP_file $INC_file");

print "\nYour inConstr file has been updated";
print "\nNumber of duplicates: $dupcounter";
print "\nAll done\n";

close(INC);
close(TMP);
