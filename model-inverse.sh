#!/bin/sh
export PATH=$HOME/PROGRAMS.330/bin:$PATH

surf96 39
surf96 32 1.
surf96 36 0
surf96 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2 6 1 2
srfphv96
plotnps -EPS -K -F7 -W10 <SRFPHV96.PLT > figsrf1.eps
surf96 28 modl.out
surf96 39
