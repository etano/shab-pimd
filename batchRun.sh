#!/bin/sh
for i in {1..50}
do
  $i | xargs -n2 -P6 sh run.sh
done
