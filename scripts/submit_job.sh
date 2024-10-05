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
    mkdir "$test_dir"
    EXIT_STATUS=$?
    if [ $EXIT_STATUS -ne 0 ]; then
        echo "finename $test_dir already exists"
        exit $EXIT_STATUS
    fi
    cp -r ./scripts/job.sh a.out "$test_dir"
    cd "$test_dir"

    pjsub job.sh
	cd -
else
    echo "Make failed with status $EXIT_STATUS"
    exit $EXIT_STATUS
fi