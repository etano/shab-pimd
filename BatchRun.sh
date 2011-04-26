#!/bin/sh
firstVar="1"
for var in "$@"
do
  if [ $firstVar -eq "1" ]
  then
    echo "Using" $1 "Processes"
    firstVar=0
  else
    nLines=$(wc -l < "inputs/$var" )
    echo "$var" "has" $nLines "input lines"

    line="1"
    echo "$line " > BatchRun.txt
    while [ $line -lt $nLines ]
    do
      line=$[$line+1]
      echo "$line " >> BatchRun.txt
    done

    cat BatchRun.txt | xargs -n1 -P$1 sh Run.sh "$var"
  fi
done
