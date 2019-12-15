#!/bin/sh

root -l -b -q 'plot.C("rhobpHist",5.0e6, 0.1, 0.16, 0.66)'
root -l -b -q 'plot.C("rerrHist",70.0e6, 0.1, 0.56, 0.66)'
root -l -b -q 'plot.C("numpcHist",70.0e6, 0.1, 0.56, 0.66)'
root -l -b -q 'plot.C("r1dwidecutWHist",9.0e6, 0.1, 0.56, 0.66)'
root -l -b -q 'plot.C("r1dwidecutWHist",0.9e6, 0.1, 0.56, 0.66)'

exit
