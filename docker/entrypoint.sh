#!/bin/sh

set -e

umask 000
mkdir /src/tmp && cd /src/tmp && git clone https://bitbucket.org/berkeleylab/upcxx.git

cd /src

bash "$@"
