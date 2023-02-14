#!/bin/sh
export PATH=$HOME/PROGRAMS.330/bin:$PATH

sprep96 -M model-from-python.d -R -FARR freqdata -NMOD 2
sdisp96 
sdpsrf96 -R -TXT
