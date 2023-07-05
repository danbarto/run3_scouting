import ROOT
import os,sys
from DataFormats.FWLite import Events, Handle

def getMiniIsoDR(mvec,minDR=0.05,maxDR=0.2,kt=10.0):
    return min(maxDR, max(minDR, kt / mvec.Pt()))

def getPFIsolation(mvec,pfcands,maxDR=0.3,vetoDR=0.01,maxdz=0.1,domini=False,domindrmpf=True):
    mindrmpf=1e6
    chiso,nhiso,phiso,puiso=0.0,0.0,0.0,0.0
    if domini:
        maxDR = getMiniIsoDR(mvec)
    for p in pfcands:
        pfvec = ROOT.TLorentzVector()
        pfvec.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),0.0)
        dR = mvec.DeltaR(pfvec)
        if dR<mindrmpf:
            mindrmpf = dR
        if dR<vetoDR or dR>maxDR:
            continue
        frompv = not (abs(p.vertex())>0 or p.dz()>maxdz)
        isch,isnh,isph,ispu=False,False,False,False
        if abs(p.pdgId())==211:
            if frompv:
                isch  = True
                chiso = chiso+p.pt()
            else:
                ispu=True
                puiso = puiso+p.pt()
        if abs(p.pdgId())==130:
            isnh  = True
            nhiso = nhiso+p.pt()
        if abs(p.pdgId())==22:
            isph  = True
            phiso = phiso+p.pt()
    if domindrmpf:
        return chiso,nhiso,phiso,puiso,mindrmpf
    else:
        return chiso,nhiso,phiso,puiso
