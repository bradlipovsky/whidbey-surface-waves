#!/bin/sh

cat > model.d << EOF 
MODEL 
TEST MODEL 
ISOTROPIC 
KGS 
FLAT EARTH 
1-D 
CONSTANT VELOCITY 
LINE08 
LINE09 
LINE10 
LINE11
HR VP VS RHO QP QS ETAP ETAS FREFP FREFS
0.1 1.5 0.0 1.0 0.0 0.0 0.0 0.0 1.0 1.0
0.1 1.3 0.05 1.5 0.0 0.0 0.0 0.0 1.0 1.0
0.1 1.6 0.1 1.5 0.0 0.0 0.0 0.0 1.0 1.0
1.0 2.0 1.0 2.0 0.0 0.0 0.0 0.0 1.0 1.0 
40. 6.0 3.5 2.8 0.0 0.0 0.0 0.0 1.0 1.0
00. 8.0 4.7 3.3 0.0 0.0 0.0 0.0 1.0 1.0 
EOF

cat > freqdata << EOF
0.1
0.2
0.3
0.4
0.5
0.6
0.7
0.8
0.9
1.0
1.1
1.2
1.3
1.4
1.5
1.6
1.7
1.8
1.9
2.0:
EOF

./sprep96 -M model.d -R -FARR freqdata
./sdisp96 
./sdpsrf96 -R -TXT
