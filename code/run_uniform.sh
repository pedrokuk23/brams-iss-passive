#!/bin/sh
# exits with a non-zero status if any of the parallel jobs fails (xargs returns 123)
set -eu
cd "$(dirname "$0")"
jobs="BEOTTI:20261103 BEOUDS:20261103 FRMAUB:20261103 BEHUMA:20261103 BEOTTI:20261104 BEOUDS:20261104 FRMAUB:20261104 BEHUMA:20261104 BEUCCE:20261104"
echo $jobs | tr ' ' '\n' | xargs -P 2 -I{} sh -c 'st=$(echo {} | cut -d: -f1); d=$(echo {} | cut -d: -f2); python3 air_uniform.py $st $d >> uniform.log 2>&1'
echo DONE >> uniform.log
