#!/bin/sh
cat batchRun.txt | xargs -n1 -P6 sh run.sh $1
