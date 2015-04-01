#!/bin/bash
# This is simply part of a piped command from bedtools for use in BDS

awk '{print $1 "\t" $2 "\t" $3 "\n" $4 "\t" $5 "\t" $6}' $1
