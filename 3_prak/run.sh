#!/bin/bash
gcc main.c -lpthread -lm
        for (( counter=1; counter<29; counter++ ))
        do
               ./a.out $counter >> output
        done
