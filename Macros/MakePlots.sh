#!/bin/sh

root -l -b -q 'plot.C("rhobpHist",5.0e6, 0.1, 0.16, 0.66, false, "Rho (cm)", "Conversions per bin")'
root -l -b -q 'plot.C("mggCutHist",4000.0, 0.0, 0.56, 0.66, false, "Di-photon mass (GeV)", "Combinations per bin")'
root -l -b -q 'plot.C("rerrHist",20.0e7, 2000.0, 0.56, 0.66, true, "Radial uncertainty (cm)", "Conversions per bin")'
root -l -b -q 'plot.C("zerrHist",20.0e7, 100.0, 0.56, 0.66, true, "z uncertainty (cm)", "Conversions per bin")'
root -l -b -q 'plot.C("phierrHist",24.0e7, 1000.0, 0.56, 0.66, true, "phi uncertainty (rad)", "Conversions per bin")'


root -l -b -q 'plot.C("numpcHist",70.0e6, 0.1, 0.56, 0.66)'
root -l -b -q 'plot.C("r1dwidecutWHist",0.9e6, 0.1, 0.56, 0.66)'
root -l -b -q 'plot.C("r1dwidecutWHist",9.0e6, 0.1, 0.56, 0.66, false, "R (cm)", "Conversions per bin" )'

root -l -b -q 'plot.C("mpairHist",6.0e7, 0.1, 0.56, 0.66, false, "Pair mass (GeV)", "Conversions per bin")'
root -l -b -q 'plot.C("dminHist",6.0e7, 0.1, 0.56, 0.66, false, "distOfMinimumApproach (cm)", "Conversions per bin")'
root -l -b -q 'plot.C("dminHist",1.0e8, 1000.0, 0.56, 0.66, true, "distOfMinimumApproach (cm)", "Conversions per bin")'
root -l -b -q 'plot.C("dphiHist",3.5e7, 0.1, 0.56, 0.66, false, "dPhiTracksAtVtx (rad)","Conversions per bin")'
root -l -b -q 'plot.C("dcotthetaHist",8.0e7, 0.1, 0.56, 0.66, false, "pairCotThetaSeparation", "Conversions per bin")'
root -l -b -q 'plot.C("dcotthetaHist",1.2e8, 1.0, 0.56, 0.66, true, "pairCotThetaSeparation", "Conversions per bin")'

root -l -b -q 'plot.C("costhetaHist",0.22e7, 0.1, 0.36, 0.66, false, "cos(theta)", "Conversions per bin")'
root -l -b -q 'plot.C("zHist",3.0e6, 0.1, 0.61, 0.66, false, "z of PC vertex (cm)", "Conversions per bin")'
root -l -b -q 'plot.C("pfitHist",1.2e7, 0.1, 0.56, 0.66, false, "Fit probability", "Conversions per bin")'
root -l -b -q 'plot.C("phiHist2",0.9e7, 0.0, 0.56, 0.66, false, "Photon conversion phi (rad)", "Conversions per bin")'

root -l -b -q 'plot.C("r1dwidecutWHist",2.0e7, 3000.0, 0.56, 0.66, true, "R (cm)", "Conversions per bin")'
root -l -b -q 'plot.C("r1dwidecutDHist",3.0e7, 900.0, 0.56, 0.66, true)'
root -l -b -q 'plot.C("r1dwidecutDDHist",3.0e7, 900.0, 0.56, 0.66, true)'

root -l -b -q 'plot.C("r1dwidecutDDHist",5.0e6, 0.1, 0.56, 0.66, false)'

root -l -b -q 'plot.C("pTHist2",3.0e7, 30.0, 0.56, 0.66, true, "Photon conversion pT (GeV)", "Conversions per bin")'
root -l -b -q 'plot.C("numpvWHist",3.0e7, 0.1, 0.56, 0.66, true, "nPV", "Events per bin")'

root -l -b -q 'plot.C("numpcHist",1.0e8, 10.0, 0.56, 0.66, true, "nPC", "Conversions per bin")'

exit
