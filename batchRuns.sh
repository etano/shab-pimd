#!/bin/sh
cat batchRuns.txt | xargs -n13 -P6 sh run.sh
