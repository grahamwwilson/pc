#!/bin/sh

#root -l -b -q 'plot.C("asym1",1.2e6, 0.0, 0.36, 0.26, false, "Positron pT Fraction", "Conversions per bin")'
#root -l -b -q 'plot.C("asym2",0.6e6, 0.0, 0.36, 0.26, false, "Positron pT Fraction", "Conversions per bin")'
#root -l -b -q 'plot.C("asym4",0.22e6, 0.0, 0.36, 0.26, false, "Positron pT Fraction", "Conversions per bin")'
#root -l -b -q 'plot.C("asym8",0.15e6, 0.0, 0.36, 0.66, false, "Positron pT Fraction", "Conversions per bin")'
#root -l -b -q 'plot.C("asym16",0.13e6, 0.0, 0.36, 0.26, false, "Positron pT Fraction", "Conversions per bin")'

root -l -b -q 'plot2.C("r1dwidecutHist",4.0e6, 5.0e2, 0.6, 0.66, true, "Conversion Radius (cm)", "Conversions per bin")'


exit
