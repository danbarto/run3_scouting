import ROOT
import os,sys,json
from datetime import date    
import numpy as np
import copy
sys.path.append('utils')
import plotUtils


indir = sys.argv[1]
if not os.path.isfile("%s/histograms_all.root"%indir):
    os.system('hadd '+indir+'/histograms_all.root $(find '+indir+' -name "histograms_*.root")')
    
user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")
outdir = 'plots_%s'%today
#outdir = indir
if not os.path.exists(outdir):
    os.makedirs(outdir)
os.system('cp '+os.environ.get("PWD")+'/utils/index.php '+outdir)

fin = ROOT.TFile.Open("%s/histograms_all.root"%indir,"r")
fin.cd()
listkeys = fin.GetListOfKeys()
size = listkeys.GetSize()
h1dn = []
h2dn = []
for i in range(0,size):
    if "TH1" in listkeys.At(i).GetClassName():
        h1dn.append(listkeys.At(i).GetName())
    elif "TH2" in listkeys.At(i).GetClassName():
        h2dn.append(listkeys.At(i).GetName())
h1d = []
h2d = []
for hn in h1dn:
    h1d.append(copy.deepcopy(fin.Get(hn)))
for hn in h2dn:
    h2d.append(copy.deepcopy(fin.Get(hn)))

# Draw histograms
unityArea = True
doLogy = True

# Labels
yearenergy = "Run2022D (13.6 TeV)"
cmsExtra = "Preliminary"
drawCMSOnTop = True
#
latex = ROOT.TLatex()
latex.SetTextFont(42)
latex.SetTextAlign(31)
latex.SetTextSize(0.04)
latex.SetNDC(True)
#
latexCMS = ROOT.TLatex()
latexCMS.SetTextFont(61)
latexCMS.SetTextSize(0.04)
latexCMS.SetNDC(True)
#
latexCMSExtra = ROOT.TLatex()
latexCMSExtra.SetTextFont(52)
latexCMSExtra.SetTextSize(0.04)
latexCMSExtra.SetNDC(True)

ROOT.gStyle.SetOptStat(0)
can = ROOT.TCanvas("can","",600,600)
for h in h1d:
    #if "lxy" in h.GetName() or "l3d" in h.GetName():
    #    h.Rebin(2)
    #    ytitle = h.GetYaxis().GetTitle()
    #    ytitle = ytitle.replace("0.1","0.2")
    #    h.GetYaxis().SetTitle(ytitle)
    if "_type" not in h.GetName():
        xmin=None
        xmax=None
        #if "lxy" in h.GetName():
        #    xmax=20.0
        plotUtils.PutUnderflowInFirstBin(h,xmin)
        plotUtils.PutOverflowInLastBin(h,xmax)
    if xmin!=None:
        h.GetXaxis().SetRangeUser(h.GetXaxis().GetBinLowEdge(h.GetXaxis().FindBin(xmin)),h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))
    if xmax!=None:
        h.GetXaxis().SetRangeUser(h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(h.GetXaxis().FindBin(xmax)))
    h.SetBinErrorOption(ROOT.TH1.kPoisson)
    h.GetYaxis().SetLabelSize(0.025)
    h.GetYaxis().SetMaxDigits(3)
    h.GetYaxis().SetRangeUser(0.0, 1.1*h.GetMaximum())
    h.SetLineWidth(2)
    ytitle = h.GetYaxis().GetTitle()
    if doLogy:
        h.GetYaxis().SetRangeUser(0.9, 2.0*h.GetMaximum())
    if unityArea:
        ytitle = ytitle.replace("Events", "Fraction of events")
        if h.Integral(0,-1)>0:
            h.Scale(1.0/h.Integral(0,-1))
        h.GetYaxis().SetRangeUser(0.0,1.1*h.GetMaximum())
        if doLogy:
            h.GetYaxis().SetRangeUser(0.9/h.GetEntries(), 2.0*h.GetMaximum())
    if "hsv" in h.GetName() and not "mind" in h.GetName() and not "maxd" in h.GetName():
        ytitle = ytitle.replace("Events", "Number of SVs")
        ytitle = ytitle.replace("events", "SVs")
    if "hmuon" in h.GetName() and h.GetName()!="hmuon_mindr" and h.GetName()!="hmuon_maxdr":
        ytitle = ytitle.replace("Events", "Number of muons")
        ytitle = ytitle.replace("events", "muons")
    h.GetYaxis().SetTitle(ytitle)
    can.cd()
    if doLogy:
        ROOT.gPad.SetLogy()
    h.Draw("E")
    can.Update()
    latex.DrawLatex(0.90, 0.91, yearenergy);
    if drawCMSOnTop:
        latexCMS.DrawLatex(0.1,0.91,"CMS");
        latexCMSExtra.DrawLatex(0.19,0.91, cmsExtra);
    else:
        latexCMS.DrawLatex(0.13,0.86,"CMS");
        latexCMSExtra.DrawLatex(0.13,0.815, cmsExtra);
    ROOT.gPad.RedrawAxis()
    outname = "%s"%h.GetName()
    if "h_" in outname:
        outname = outname.replace("h_","")
    elif "hsv" in outname or "hmuon" in outname:
        outname = outname.replace("hsv","sv")
        outname = outname.replace("hmuon","muon")
    can.SaveAs("%s/%s.png"%(outdir,outname))
    can.Clear()

for h in h2d:
    can.cd()
    ROOT.gPad.SetLogy(0)
    h.GetYaxis().SetLabelSize(0.025)
    h.GetYaxis().SetMaxDigits(3)
    ztitle = h.GetZaxis().GetTitle()
    if unityArea:
        ztitle = ztitle.replace("Events", "Fraction of events [%]")
        if h.Integral(0,-1,0,-1)>0:
            h.Scale(100.0/h.Integral(0,-1,0,-1))
        h.GetZaxis().SetRangeUser(0.0,1.1*h.GetMaximum())
    if "nmuons" in h.GetName():
        ztitle = ztitle.replace("Events", "Number of muons")
        ztitle = ztitle.replace("events", "muons")
    h.GetZaxis().SetTitle(ztitle)
    can.cd()
    if unityArea:
        ROOT.gStyle.SetPaintTextFormat(".2f");
        h.Draw("colz,text")
    can.Update()
    palette = h.GetListOfFunctions().FindObject("palette")
    if palette:
        palette.SetX2NDC(0.925)
    can.Update()
    h.GetZaxis().SetLabelSize(0.025)
    h.GetZaxis().SetMaxDigits(2)
    h.GetZaxis().SetRangeUser(0.0, 1.1*h.GetMaximum())
    latex.DrawLatex(0.90, 0.91, yearenergy);
    if drawCMSOnTop:
        latexCMS.DrawLatex(0.1,0.91,"CMS");
        latexCMSExtra.DrawLatex(0.19,0.91, cmsExtra);
    else:
        latexCMS.DrawLatex(0.13,0.86,"CMS");
        latexCMSExtra.DrawLatex(0.13,0.815, cmsExtra);
    ROOT.gPad.RedrawAxis()
    outname = "%s"%h.GetName()
    if "h_" in outname:
        outname = outname.replace("h_","")
    can.SaveAs("%s/%s.png"%(outdir,outname))
    can.Clear()

del can
