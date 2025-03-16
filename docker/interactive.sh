#!/bin/sh

docker run --rm -it -v $(pwd)/..:/src cs267 
#docker run --gpus all --rm -it -v $(pwd)/..:/src cs267 
