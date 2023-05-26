#!/bin/bash

if (make 1> /dev/null); 
then
    ./main 4
    ./main 8
    ./main 16
    ./main 32
    ./main 64
    ./main 128
    ./main 256
    ./main 512

    make clean 1> /dev/null
fi
