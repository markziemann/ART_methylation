#!/bin/bash


ls *bam | parallel -j12 macs2 callpeak -t {} --outdir {}_peak -n {}_macs
