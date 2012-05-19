#!/bin/tcsh
source `which qaConfig.csh`

################################
#  
#  03-31-2008
#  Ann Zweig
#
#  Find the organism name given the assembly name
#
################################

set date=""

if ( $#argv < 1  | $#argv > 2 ) then
  echo
  echo " Finds the organism name given the assembly name" 
  echo "  usage: assemblyName [date]"
  echo "         will accept name with or without digit"
  echo "         use 'date' to also retrieve assembly date"
  echo "         (e.g. 'ornAna2' or 'ornAna')"
  echo
  exit 1
else
  set db=$argv[1]
endif

if ( $#argv == 2 ) then
  set date=$argv[2]
  if ( $date != "date" ) then
    echo "\n  error. use the word 'date' as second parameter\n"  
    exit 1
  endif
  set date=", description"
endif

if ( "$HOST" != "hgwdev" ) then
 echo "\n ERROR: you must run this script on dev!\n"
 exit 1
endif

hgsql -t -e "SELECT name, organism $date FROM dbDb WHERE NAME LIKE '$db%' \
  ORDER BY name" hgcentraltest  | tail -n+3

exit 0
