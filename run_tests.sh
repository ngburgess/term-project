#!/bin/bash

mkdir -p output
mkdir -p logs

# Compile
gcc adaptive_noise_filter.c -o adaptive_noise_filter -lm

if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

# Clear previous log
> logs/results.txt

for img in ./test_images/*/*.pgm; do
    if [[ "$img" == *clean* ]]; then
        continue
    fi

    base=$(basename "$img" .pgm)

    echo "===================================" | tee -a logs/results.txt
    echo "Processing: $base" | tee -a logs/results.txt
    echo "===================================" | tee -a logs/results.txt

    ./adaptive_noise_filter \
        "$img" \
        "./output/${base}_adaptive.pgm" \
        2>&1 | tee -a logs/results.txt

    echo "" | tee -a logs/results.txt

done

echo "Finished all tests."
