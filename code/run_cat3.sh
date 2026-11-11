#!/bin/sh
# exits with a non-zero status if any of the parallel jobs fails (xargs returns 123)
set -eu
cd "$(dirname "$0")"
cat night_files.txt | xargs -P 2 -I{} sh -c 'python3 catscan3.py ${WAVDIR:-../data/raw_wav/}{} >> cat3.log 2>&1'
echo DONE >> cat3.log
