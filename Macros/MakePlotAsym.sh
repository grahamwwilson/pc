#!/bin/sh

#root -l -b -q 'plot.C("asym1",1.2e6, 0.0, 0.36, 0.26, false, "Positron pT Fraction", "Conversions per bin")'
#root -l -b -q 'plot.C("asym2",0.6e6, 0.0, 0.36, 0.26, false, "Positron pT Fraction", "Conversions per bin")'
#root -l -b -q 'plot.C("asym4",0.22e6, 0.0, 0.36, 0.26, false, "Positron pT Fraction", "Conversions per bin")'
#root -l -b -q 'plot.C("asym8",0.15e6, 0.0, 0.36, 0.66, false, "Positron pT Fraction", "Conversions per bin")'
#root -l -b -q 'plot.C("asym16",0.13e6, 0.0, 0.36, 0.26, false, "Positron pT Fraction", "Conversions per bin")'

#root -l -b -q 'plot2.C("r1dwidecutHist",4.0e6, 5.0e2, 0.6, 0.66, true, "Conversion Radius (cm)", "Conversions per bin")'

root -l -b -q 'plot2.C("conversionCandidateMassHist2",2.0e7, 5.0e2, 0.6, 0.66, true, "Pair Mass (GeV)", "Conversions per bin")'

root -l -b -q 'plot2.C("mpairHist",2.0e7, 5.0e2, 0.6, 0.66, true, "Pair Mass (GeV)", "Conversions per bin")'

root -l -b -q 'plot2.C("alpha16",10000.0, 0.0, 0.6, 0.66, false)'

root -l -b -q 'plot2.C("r0",1.5e6, 0.0, 0.6, 0.66, false, "Radius (cm)", "Conversions per bin")'
root -l -b -q 'plot2.C("r1p",0.7e6, 0.0, 0.6, 0.66, false, "Radius (cm)", "Conversions per bin")'
root -l -b -q 'plot2.C("r2p",0.2e6, 0.0, 0.6, 0.66, false, "Radius (cm)", "Conversions per bin")'
root -l -b -q 'plot2.C("r4p",0.07e6, 0.0, 0.6, 0.66, false, "Radius (cm)", "Conversions per bin")'
root -l -b -q 'plot2.C("r8p",0.025e6, 0.0, 0.6, 0.66, false, "Radius (cm)", "Conversions per bin")'
root -l -b -q 'plot2.C("r16",0.01e6, 0.0, 0.6, 0.66, false, "Radius (cm)", "Conversions per bin")'

root -l -b -q 'plot2.C("minptHist",2.0e6, 1.0, 0.6, 0.66, true, "minimum pT of pair (GeV)", "Conversions per bin")'
root -l -b -q 'plot2.C("maxptHist",2.0e6, 1.0, 0.6, 0.66, true, "maximum pT of pair (GeV)", "Conversions per bin")'

exit
