//       #include "histset2enums.h" 
// bookeeping enumeration: if we do this we don't need to worry about hist pointer copies and merging
       enum th1d_ids{id_ptHist, id_pzHist, id_numpcHist, id_numpvHist,id_numpvWHist,id_rerrHist, id_phierrHist, id_zerrHist,
            id_r1dHist2, id_r1dHist3, id_r1dHist4, id_r1dHist, id_r1dcutHist, id_r1dlowPUHist, id_r1dmedPUHist, id_r1dhiPUHist, 
            id_r1dlowPUcutHist, id_r1dmedPUcutHist, id_r1dhiPUcutHist, id_r1dwideHist, id_r1dwidecutHist, id_r1dwidecutWHist,
            id_r1dwidecutPSHist, id_r1dwidelowPUHist, id_r1dwidemedPUHist, id_r1dwidehiPUHist, id_r1dwidelowPUcutHist,
            id_r1dwidemedPUcutHist, id_r1dwidehiPUcutHist, id_rhobpHist, id_rbpHist, id_mggHist, id_mggCutHist,
            id_numnopcHist, id_numpvnopcHist, id_phiHist, id_mggallHist, id_pfitHist, id_zHist, id_costhetaHist,
            id_pTHist, id_EHist, id_pTHist2, id_EHist2, id_phiHist2, id_runHist, id_isdataHist, id_nPUHist, id_PUHist, id_wtHist,
            id_wwtHist, id_numpvUWHist, id_dminHist, id_dphiHist, id_mpairHist, id_dcotthetaHist, id_r1dwidecutDHist, id_r1dwidecutDDHist,
            id_r1dwidecutNomHist, id_rnomHist,
            id_xplus1Hist, id_xplus2Hist, id_xplus4Hist, id_xplus8Hist, id_xplus16Hist,
            id_alpha1Hist, id_alpha2Hist, id_alpha4Hist, id_alpha8Hist, id_alpha16Hist,
            id_q0Hist, id_q1Hist, id_qtotHist,
            id_zcutHist, id_zcutHist2, id_rendcapHist,
            id_conversionCandidateMassHist, id_conversionCandidateMassHist2,
            id_lambdaCandidateMassHist,id_lambdabarCandidateMassHist,id_lambdasCandidateMassHist,
            id_lambdasBkgdMassHist,id_lambdasSignalMassHist, id_lambdasBkgdMassHistR1,id_lambdasBkgdMassHistR2,id_lambdasBkgdMassHistR3,
            id_KShortMassHist,id_KShortBkgdMassHistR1,id_KShortBkgdMassHistR2, id_KShortBkgdMassHistR3,id_KShortBkgdMassHist,
            id_AP_pTminHist, id_AP_pTmaxHist, id_AP_pTaveHist, id_AP_alphaHist,
            id_alphaBkgdHist, id_alphaSignalHist, id_alphaBkgdHistR1, id_alphaBkgdHistR2, id_alphaBkgdHistR3,
            id_nconvHist, id_nassignedHist, id_nnonassignedHist,
            numTH1Hist};
       enum th2d_ids{id_pxpyHist, id_xyHist, id_xywideHist, id_rphiHist, 
            id_rzHist, id_rzHist2, id_rzHist3, id_rzHist4,
            id_xycutHist, id_xywidecutHist, id_xywidecutHist2,
            id_npv_rcutHist, id_mgg2Hist, id_mggRCutHist, id_npc_npvHist, id_rhophiHist, id_AP_pT_alphaHist,
            numTH2Hist};
