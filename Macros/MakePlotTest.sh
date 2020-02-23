#!/bin/sh

root -l -b -q 'plot.C("r1dwidecutPSHist",2.0e7, 3000.0, 0.56, 0.66, true, "R_{PS} (cm)", "Conversions per 0.1cm bin")'
root -l -b -q 'plot.C("r1dwidecutWHist",2.0e7, 3000.0, 0.56, 0.66, true, "R_{BPIX} (cm)", "Conversions per 0.1cm bin")'
root -l -b -q 'plot.C("rhobpHist",5.0e6, 0.1, 0.16, 0.66, false, "R_{BPIPE} (cm)", "Conversions per 0.05cm bin")'

exit
