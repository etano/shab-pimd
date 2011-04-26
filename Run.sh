#!/bin/sh
echo "Running simulation" $2 "from" $1 "..."
./shab-pimd "inputs/$1" "$2" > "data/output/$1-$2.out"

