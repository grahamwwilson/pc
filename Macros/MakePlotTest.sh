#!/bin/sh

root -l -b -q 'plot.C("r1dwidecutNomHist",2.0e7, 3000.0, 0.56, 0.66, true, "R_{nominal} (cm)", "Conversions per 0.1cm bin")'
root -l -b -q 'plot.C("r1dwidecutPSHist",2.0e7, 3000.0, 0.56, 0.66, true, "R_{PS} (cm)", "Conversions per 0.1cm bin")'
root -l -b -q 'plot.C("r1dwidecutWHist",2.0e7, 3000.0, 0.56, 0.66, true, "R_{BPIX} (cm)", "Conversions per 0.1cm bin")'
root -l -b -q 'plot.C("rbpHist",5.0e6, 0.1, 0.16, 0.66, false, "R_{BPIX} (cm)", "Conversions per 0.05cm bin")'
root -l -b -q 'plot.C("rhobpHist",5.0e6, 0.1, 0.16, 0.66, false, "R_{BPIPE} (cm)", "Conversions per 0.05cm bin")'
root -l -b -q 'plot.C("rnomHist",5.0e6, 0.1, 0.16, 0.66, false, "R_{nominal} (cm)", "Conversions per 0.05cm bin")'

root -l -b -q 'plotratio.C("r1dwidecutWHist",4.0, 0.25, 0.36, 0.66, false, "R_{BPIX} (cm)", "Data/MC per 0.1cm bin")'
root -l -b -q 'plotratio.C("r1dwidecutPSHist",4.0, 0.25, 0.36, 0.66, false, "R_{PS} (cm)", "Data/MC per 0.1cm bin")'
root -l -b -q 'plotratio.C("r1dwidecutNomHist",4.0, 0.25, 0.36, 0.66, false, "R_{nominal} (cm)", "Data/MC per 0.1cm bin")'
root -l -b -q 'plotratio.C("rhobpHist",4.0, 0.25, 0.40, 0.66, false, "R_{BPIPE} (cm)", "DATA/MC per 0.05cm bin")'
root -l -b -q 'plotratio.C("rbpHist",4.0, 0.25, 0.40, 0.66, false, "R_{BPIX} (cm)", "DATA/MC per 0.05cm bin")'
root -l -b -q 'plotratio.C("rnomHist",4.0, 0.25, 0.56, 0.66, false, "R_{nominal} (cm)", "DATA/MC per 0.05cm bin")'

exit
