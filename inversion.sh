#!/bin/sh
export PATH=$HOME/PROGRAMS.330/bin:$PATH

echo "--->  Cleaning up"
surf96 39

echo "--->  Define damping"
surf96 32 1.

echo "--->  Differential smoothing"
surf96 36 1

echo "--->  Set up five iterations"
surf96 1 2 6 1 2 6 1 2
#surf96 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2

echo "--->  Make plots"
#srfphv96
#plotnps -EPS -K -F7 -W10 < SRFPHV96.PLT > figsrf1.eps

echo "---> Save"
surf96 28 modl.out
