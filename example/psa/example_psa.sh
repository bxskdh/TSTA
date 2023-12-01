#!/bin/bash

##Pairwise sequence alignment##
../TSTA_psa -M 2 -X -3 -E -2 -O -4 -T 30 -S 30 seq/seqa1.fa seq/seqb1.fa > result.txt
echo done!
