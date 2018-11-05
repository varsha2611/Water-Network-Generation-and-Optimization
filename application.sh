#!/usr/bin/env bash
$Input = $1
$params = $2

python generate_inputs.py -i $Input -p $params

module add gnu-parallel
parallel -j100 <inputs.txt
