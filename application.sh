#!/usr/bin/env bash
Input_file = $1
Params = $2
python generate_inputs.py -i $1 -p $2

module add gnu-parallel
parallel -j100 <inputs.txt