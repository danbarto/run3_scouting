#!/usr/bin/env python3
import uproot
import pandas as pd
import numpy as np

from coffea.nanoevents.methods import vector
import awkward as ak
ak.behavior.update(vector.behavior)

import hist
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import mplhep as hep #matplotlib wrapper for easy plotting in HEP
plt.style.use(hep.style.CMS)


smaller_dtypes = [
    ["dimuon_mass","float32"],
    ["dimuon_pt","float32"],
    ["DV_chi2prob","float32"],
    ["DV_ndof","int8"],
    ["DV_redchi2","float32"],
    ["DV_layerPixel","int8"],
    ["Muon1_charge","int8"],
    ["Muon1_excesshits","int8"],
    ["Muon1_m","float32"],
    ["Muon1_mass","float32"],
    ["Muon1_nExcessPixelHits","int8"],
    ["Muon1_nExpectedPixelHits","int8"],
    ["Muon1_nMatchedStations","int8"],
    ["Muon1_nTrackerLayersWithMeasurement","int8"],
    ["Muon1_nValidMuonHits","int8"],
    ["Muon1_nValidPixelHits","int8"],
    ["Muon1_nValidStripHits","int8"],
    ["Muon2_charge","int8"],
    ["Muon2_excesshits","int8"],
    ["Muon2_m","float32"],
    ["Muon2_mass","float32"],
    ["Muon2_nExcessPixelHits","int8"],
    ["Muon2_nExpectedPixelHits","int8"],
    ["Muon2_nMatchedStations","int8"],
    ["Muon2_nTrackerLayersWithMeasurement","int8"],
    ["Muon2_nValidMuonHits","int8"],
    ["Muon2_nValidPixelHits","int8"],
    ["Muon2_nValidStripHits","int8"],
    ["categ","int8"],
    ["luminosityBlock","int32"],
    ["nDV","int8"],
    ["nDV_good","int8"],
    ["nDV_raw","int8"],
    ["nGenMuon","int8"],
    ["nGenPart","int16"],
    ["nJet","int8"],
    ["nMuon","int8"],
    ["nMuon_good","int8"],
    ["nMuon_raw","int8"],
    ["nPV","int8"],
    ["nPVM","int8"],
    ["run","int32"],
]

def make_df(fname,entrystart=None,entrystop=None,bigrho=True,iso=True):

    #events = uproot.open(f"{fname}:Events")
    #df = events.array(library="pd")
    tree = uproot.open(fname)["Events"]
    df = tree.arrays(library="pd", filter_name=["/n(DV|Jet|PV|PVM|Muon|GenMuon)$/","pass_*",
                     "/BS_(x|y|z)$/",
                     "/Muon1_n(Valid|Matched|Tracker|Expected).*/",
                     "/Muon1_(pt|eta|phi|m|trackIso|charge|dz.*|dxy.*|chi2|ndof|drjet|pass*)$/",
                     "/Muon2_n(Valid|Matched|Tracker|Expected).*/",
                     "/Muon2_(pt|eta|phi|m|trackIso|charge|dz.*|dxy.*|chi2|ndof|drjet|pass*)$/",
                     "/DV_(chi2|ndof|rho.*|inPixel.*|x|y|z|xError|yError|zError|pass*)$/",
                     "run","luminosityBlock","event",
                     "dimuon_*","cosphi*","absdphi*","minabs*","logabs*",
                     "L1_*", "lxy",
                    ])

    sel = df["pass_baseline"]

    # flatten into dataframe and require `sel`
    df = df[sel]
    #for k in arrs.fields:
    #    if any(k.startswith(y) for y in ["n","pass_","BS_","MET_","run","lumi","event","L1_",
    #                                    "dimuon","cosphi","absdphi","minabs","logabs"]):
    #        df[k] = arrs[k][sel]
    #    if k.startswith("DV_"):
    #        df[k] = arrs[k][sel][:,0]
    #    if k.startswith("Muon_"):
    #        df[k.replace("Muon_","Muon1_")] = arrs[k][sel][:,0]
    #        df[k.replace("Muon_","Muon2_")] = arrs[k][sel][:,1]

    mu1 = ak.zip({
        "pt": df["Muon1_pt"],
        "eta": df["Muon1_eta"],
        "phi": df["Muon1_phi"],
        "mass": df["Muon1_m"],
        }, with_name="PtEtaPhiMLorentzVector"
    )

    mu2 = ak.zip({
        "pt": df["Muon2_pt"],
        "eta": df["Muon2_eta"],
        "phi": df["Muon2_phi"],
        "mass": df["Muon2_m"],
        }, with_name="PtEtaPhiMLorentzVector"
    )

    df["dimuon_deltaeta"] = np.abs(mu1.eta - mu2.eta)
    df["dimuon_deltar"] = mu1.delta_r(mu2)
    df["DV_xyErrorMax"] = np.maximum(df["DV_xError"],df["DV_yError"])


    df["Muon1_excesshits"] = df.eval("Muon1_nValidPixelHits-Muon1_nExpectedPixelHits")
    df["Muon2_excesshits"] = df.eval("Muon2_nValidPixelHits-Muon2_nExpectedPixelHits")
    df["pass_excesshits"] = df.eval("DV_rhoCorr<3.5 or (Muon1_excesshits<=0 and Muon2_excesshits<=0)")


    for name,dtype in smaller_dtypes:
        if name not in df.columns: continue
        df[name] = df[name].astype(dtype, copy=False)

    return df


if __name__ == '__main__':


    df = make_df("../production/test_data.root")


    # muon dxy
    lxy_hist = hist.Hist(
        hist.axis.Regular(200, 0, 20, name="dxy", label=r"$L_{xy}\ (cm)$"),
    )

    lxy_hist.fill(dxy=abs(df['lxy']))
    fig, ax = plt.subplots(figsize=(8, 8))
    lxy_hist.plot1d(
        histtype="step",
        ax=ax,
    )

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi='X',
        com=13.6,
        loc=0,
        ax=ax,
        fontsize=15,
    )

    ax.set_yscale('log')
    fig.savefig(f'Lxy.png')


    # dimuon mass
    mass_hist = hist.Hist(
        hist.axis.Regular(1000, 0, 100, name="mass", label=r"$M(\mu\mu)\ (GeV)$"),
    )

    mass_hist.fill(mass=df['dimuon_mass'])
    fig, ax = plt.subplots(figsize=(8, 8))
    mass_hist.plot1d(
        histtype="step",
        ax=ax,
    )

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi='X',
        com=13.6,
        loc=0,
        ax=ax,
        fontsize=15,
    )

    ax.set_yscale('log')
    fig.savefig(f'dimuon_mass.png')



    # 2D vertex plot
    vert_hist = hist.Hist(
        hist.axis.Regular(1000, -10, 10, name="x", label=r"$x (cm)$"),
        hist.axis.Regular(1000, -10, 10, name="y", label=r"$y (cm)$"),
        hist.axis.Regular(100, -5, 5, name="z", label=r"$z (cm)$"),
    )

    vert_hist.fill(
        x = df['DV_x'],
        y = df['DV_y'],
        z = df['DV_z'],
    )

    fig, ax = plt.subplots(figsize=(8, 8))
    vert_hist[-10.0j:10.0j:2j, -10.0j:10.0j:2j, :][{'z':sum}].plot2d(
        ax=ax,
        norm=LogNorm(),
    )

    hep.cms.label(
        "Preliminary",
        data=True,
        lumi='X',
        com=13.6,
        loc=0,
        ax=ax,
        fontsize=15,
    )
    fig.savefig(f'vertices.png')
