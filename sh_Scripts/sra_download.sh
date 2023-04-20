#!/bin/bash

cat SRR_Acc_List1.txt | while read -r srr; do

prefetch $srr

fasterq-dump $srr --split-files

done


