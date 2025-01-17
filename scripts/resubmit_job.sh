#!/bin/bash

rm -rf results/debug

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

test_dir="./results/$1"

make
EXIT_STATUS=$?

if [ $EXIT_STATUS -eq 0 ]; then
    cd "$test_dir"

    pjsub job.sh
	cd -
else
    echo "Make failed with status $EXIT_STATUS"
    exit $EXIT_STATUS
fi