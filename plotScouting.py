import ROOT
import os,sys,json
from datetime import date    
import numpy as np
from DataFormats.FWLite import Events, Handle
sys.path.append('utils')
import plotUtils
import jsonUtils
import isolationUtils

oncondor = True
MUON_MASS = 0.10566

GlobalMuon = 1 << 1
TrackerMuon = 1 << 2
StandAloneMuon = 1 << 3
def isGlobalMuon(t):
    return t & GlobalMuon
def isTrackerMuon(t):
    return t & TrackerMuon
def isStandAloneMuon(t):
    return t & StandAloneMuon

doPFIso = True

isData = True
applyJSON = True
year   = "2022"
wjson  = "Golden"
injson = ""
if isData:
    if year=="2022":
        if wjson=="Golden":
            injson = "data/Cert_Collisions2022_355100_362760_Golden.json"
        elif wjson=="Muon":
            injson = "data/Cert_Collisions2022_355100_362760_Muon.json"
        else:
            print('ERROR: only Golden and Muon JSONs are available!')
    else:
        print('ERROR: only 2022 JSONs are available!')
        exit()
jsonruns = json.load(open(injson,"r"))

user = os.environ.get("USER")
today= date.today().strftime("%b-%d-%Y")
outdir = 'plots_%s'%today
if len(sys.argv)>4:
    outdir = sys.argv[4]
if not os.path.exists(outdir):
    os.makedirs(outdir)
#os.system('cp '+os.environ.get("PWD")+'/utils/index.php '+outdir)

files = []

if not oncondor:
    for f in os.listdir("/ceph/cms%s"%sys.argv[1]):
        if len(sys.argv)==3:
            thisfile="_%s.root"%sys.argv[2]
            if "output" in f and ".root" in f and thisfile in f and os.path.isfile("/ceph/cms%s/%s"%(sys.argv[1],f)):
                files.append("/ceph/cms%s/%s"%(sys.argv[1],f))
        else:
            if "output" in f and ".root" in f and os.path.isfile("/ceph/cms%s/%s"%(sys.argv[1],f)):
                files.append("/ceph/cms%s/%s"%(sys.argv[1],f))
else:
    infile = open('infiles.txt',"r")
    for l in infile.readlines():
        files.append("davs://redirector.t2.ucsd.edu:1095%s/%s"%(sys.argv[1],l.replace("\n","")))

if len(sys.argv)>3:
    index = int(sys.argv[2])
    pace  = int(sys.argv[3])
    files = files[index*pace:index*pace+pace]
print("Initialize Events class over %d input files, from %d to %d"%(len(files),index*pace,index*pace+pace))
events = Events(files)
# Muons
mhandle = Handle('std::vector<Run3ScoutingMuon>')
mlabel  = ("hltScoutingMuonPacker")
# Vertices
# displaced
svhandle = Handle('std::vector<Run3ScoutingVertex>')
svlabel  = ("hltScoutingMuonPacker","displacedVtx")
# primary
pvhandle = Handle('std::vector<Run3ScoutingVertex>')
pvlabel  = ("hltScoutingPrimaryVertexPacker","primaryVtx")
# Jets
jhandle = Handle('std::vector<Run3ScoutingPFJet>')
jlabel  = ("hltScoutingPFPacker")
# PF candidates
pfhandle = Handle('std::vector<Run3ScoutingParticle>')
pflabel  = ("hltScoutingPFPacker")
# HLT string
hlthandle = Handle('std::vector<std::string>')
hltlabel  = ("triggerMaker","hltname")
# L1 strings
l1shandle = Handle('std::vector<std::string>')
l1slabel  = ("triggerMaker","l1name")
# L1 prescales
l1phandle = Handle('std::vector<std::double>')
l1plabel  = ("triggerMaker","l1prescale")

# Histograms:
h1d = []
h2d = []

# Displaced vertices:
h_nsv = ROOT.TH1D("h_nsv","",10,0,10)
h_nsv.GetXaxis().SetTitle("Number of valid SVs")
h_nsv.GetYaxis().SetTitle("Events")
h_nsv.SetLineColor(1)
h1d.append(h_nsv)

h_nsvsel = ROOT.TH1D("h_nsvsel","",10,0,10)
h_nsvsel.GetXaxis().SetTitle("Number of selected valid SVs")
h_nsvsel.GetYaxis().SetTitle("Events")
h_nsvsel.SetLineColor(4)
h1d.append(h_nsvsel)

hsv_chi2ndof = ROOT.TH1D("hsv_chi2ndof","",100,0.0,10.0)
hsv_chi2ndof.GetXaxis().SetTitle("SV #chi^{2}/ndof")
hsv_chi2ndof.GetYaxis().SetTitle("Events / 0.1")
hsv_chi2ndof.SetLineColor(1)
h1d.append(hsv_chi2ndof)

hsv_chi2prob = ROOT.TH1D("hsv_chi2prob","",100,0.0,1.0)
hsv_chi2prob.GetXaxis().SetTitle("SV #chi^{2} probability")
hsv_chi2prob.GetYaxis().SetTitle("Events / 0.01")
hsv_chi2prob.SetLineColor(1)
h1d.append(hsv_chi2prob)

hsv_xerr = ROOT.TH1D("hsv_xerr","",50,0,0.5)
hsv_xerr.GetXaxis().SetTitle("SV x error [cm]")
hsv_xerr.GetYaxis().SetTitle("Events / 0.01 cm")
hsv_xerr.SetLineColor(1)
h1d.append(hsv_xerr)

hsv_yerr = ROOT.TH1D("hsv_yerr","",50,0,0.5)
hsv_yerr.GetXaxis().SetTitle("SV y error [cm]")
hsv_yerr.GetYaxis().SetTitle("Events / 0.01 cm")
hsv_yerr.SetLineColor(1)
h1d.append(hsv_yerr)

hsv_zerr = ROOT.TH1D("hsv_zerr","",50,0,0.5)
hsv_zerr.GetXaxis().SetTitle("SV z error [cm]")
hsv_zerr.GetYaxis().SetTitle("Events / 0.01 cm")
hsv_zerr.SetLineColor(1)
h1d.append(hsv_zerr)

hsv_lxy = ROOT.TH1D("hsv_lxy","",500,0,50)
hsv_lxy.GetXaxis().SetTitle("l_{xy} (from PV) [cm]")
hsv_lxy.GetYaxis().SetTitle("Events / 0.1 cm")
hsv_lxy.SetLineColor(1)
h1d.append(hsv_lxy)

hsvsel_lxy = ROOT.TH1D("hsvsel_lxy","",500,0,50)
hsvsel_lxy.GetXaxis().SetTitle("l_{xy} (from PV) [cm]")
hsvsel_lxy.GetYaxis().SetTitle("Events / 0.1 cm")
hsvsel_lxy.SetLineColor(4)
h1d.append(hsvsel_lxy)

hsv_l3d = ROOT.TH1D("hsv_l3d","",1000,0,100)
hsv_l3d.GetXaxis().SetTitle("l_{3D} (from PV) [cm]")
hsv_l3d.GetYaxis().SetTitle("Events / 0.1 cm")
hsv_l3d.SetLineColor(1)
h1d.append(hsv_l3d)

hsvsel_l3d = ROOT.TH1D("hsvsel_l3d","",1000,0,100)
hsvsel_l3d.GetXaxis().SetTitle("l_{3D} (from PV) [cm]")
hsvsel_l3d.GetYaxis().SetTitle("Events / 0.1 cm")
hsvsel_l3d.SetLineColor(4)
h1d.append(hsvsel_l3d)

hsvsel_mindx = ROOT.TH1D("hsvsel_mindx","",100,0,1)
hsvsel_mindx.GetXaxis().SetTitle("min D_{x} (SV_{i}, SV_{j}) [cm]")
hsvsel_mindx.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_mindx.SetLineColor(4)
h1d.append(hsvsel_mindx)

hsvsel_mindy = ROOT.TH1D("hsvsel_mindy","",100,0,1)
hsvsel_mindy.GetXaxis().SetTitle("min D_{y} (SV_{i}, SV_{j}) [cm]")
hsvsel_mindy.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_mindy.SetLineColor(4)
h1d.append(hsvsel_mindy)

hsvsel_mindxy = ROOT.TH1D("hsvsel_mindxy","",100,0,1)
hsvsel_mindxy.GetXaxis().SetTitle("min D_{xy} (SV_{i}, SV_{j}) [cm]")
hsvsel_mindxy.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_mindxy.SetLineColor(4)
h1d.append(hsvsel_mindxy)

hsvsel_mindz = ROOT.TH1D("hsvsel_mindz","",100,0,1)
hsvsel_mindz.GetXaxis().SetTitle("min d_{z} (SV_{i}, SV_{j}) [cm]")
hsvsel_mindz.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_mindz.SetLineColor(4)
h1d.append(hsvsel_mindz)

hsvsel_mind3d = ROOT.TH1D("hsvsel_mind3d","",100,0,1)
hsvsel_mind3d.GetXaxis().SetTitle("min D_{3D} (SV_{i}, SV_{j}) [cm]")
hsvsel_mind3d.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_mind3d.SetLineColor(4)
h1d.append(hsvsel_mind3d)

hsvsel_maxdx = ROOT.TH1D("hsvsel_maxdx","",100,0,1)
hsvsel_maxdx.GetXaxis().SetTitle("MAX D_{x} (SV_{i}, SV_{j}) [cm]")
hsvsel_maxdx.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_maxdx.SetLineColor(4)
h1d.append(hsvsel_maxdx)

hsvsel_maxdy = ROOT.TH1D("hsvsel_maxdy","",100,0,1)
hsvsel_maxdy.GetXaxis().SetTitle("MAX D_{y} (SV_{i}, SV_{j}) [cm]")
hsvsel_maxdy.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_maxdy.SetLineColor(4)
h1d.append(hsvsel_maxdy)

hsvsel_maxdxy = ROOT.TH1D("hsvsel_maxdxy","",100,0,1)
hsvsel_maxdxy.GetXaxis().SetTitle("MAX D_{xy} (SV_{i}, SV_{j}) [cm]")
hsvsel_maxdxy.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_maxdxy.SetLineColor(4)
h1d.append(hsvsel_maxdxy)

hsvsel_maxdz = ROOT.TH1D("hsvsel_maxdz","",100,0,1)
hsvsel_maxdz.GetXaxis().SetTitle("MAX D_{z} (SV_{i}, SV_{j}) [cm]")
hsvsel_maxdz.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_maxdz.SetLineColor(4)
h1d.append(hsvsel_maxdz)

hsvsel_maxd3d = ROOT.TH1D("hsvsel_maxd3d","",100,0,1)
hsvsel_maxd3d.GetXaxis().SetTitle("MAX D_{3D} (SV_{i}, SV_{j}) [cm]")
hsvsel_maxd3d.GetYaxis().SetTitle("Events / 0.01 cm")
hsvsel_maxd3d.SetLineColor(4)
h1d.append(hsvsel_maxd3d)

# Muons
h_nmuonsass = ROOT.TH1D("h_nmuonsass","",10,0,10)
h_nmuonsass.GetXaxis().SetTitle("Number of muons (from any SV)")
h_nmuonsass.GetYaxis().SetTitle("Events")
h_nmuonsass.SetLineColor(4)
h1d.append(h_nmuonsass)

h_nmuonssel = ROOT.TH1D("h_nmuonssel","",10,0,10)
h_nmuonssel.GetXaxis().SetTitle("Number of muons (p_{T}>3 GeV, |#eta|<2.4, from any SV)")
h_nmuonssel.GetYaxis().SetTitle("Events")
h_nmuonssel.SetLineColor(4)
h1d.append(h_nmuonssel)

hmuon_nSVs = ROOT.TH1D("hmuon_nSVs","",10,0,10)
hmuon_nSVs.GetXaxis().SetTitle("Number of muon-SV associations")
hmuon_nSVs.GetYaxis().SetTitle("Events")
hmuon_nSVs.SetLineColor(2)
h1d.append(hmuon_nSVs)

hmuon_pt = ROOT.TH1D("hmuon_pt","",50,0.0,150.0)
hmuon_pt.GetXaxis().SetTitle("Muon p_{T} [GeV]")
hmuon_pt.GetYaxis().SetTitle("Events / 3 GeV")
hmuon_pt.SetLineColor(2)
h1d.append(hmuon_pt)

hmuon_eta = ROOT.TH1D("hmuon_eta","",48,-2.4,2.4)
hmuon_eta.GetXaxis().SetTitle("Muon #eta")
hmuon_eta.GetYaxis().SetTitle("Events / 0.1")
hmuon_eta.SetLineColor(2)
h1d.append(hmuon_eta)

hmuon_phi = ROOT.TH1D("hmuon_phi","",64,-3.2,3.2)
hmuon_phi.GetXaxis().SetTitle("Muon #phi [rad]")
hmuon_phi.GetYaxis().SetTitle("Events / 0.1 rad")
hmuon_phi.SetLineColor(2)
h1d.append(hmuon_phi)

hmuon_ch = ROOT.TH1D("hmuon_ch","",3,-1.5,1.5)
hmuon_ch.GetXaxis().SetTitle("Muon charge")
hmuon_ch.GetYaxis().SetTitle("Events")
hmuon_ch.SetLineColor(2)
h1d.append(hmuon_ch)

hmuon_mindr = ROOT.TH1D("hmuon_mindr","",100,0.0,5.0)
hmuon_mindr.GetXaxis().SetTitle("min #DeltaR(#mu_{i}, #mu_{j})")
hmuon_mindr.GetYaxis().SetTitle("Events / 0.05")
hmuon_mindr.SetLineColor(2)
h1d.append(hmuon_mindr)

hmuon_maxdr = ROOT.TH1D("hmuon_maxdr","",100,0.0,5.0)
hmuon_maxdr.GetXaxis().SetTitle("MAX #DeltaR(#mu_{i}, #mu_{j})")
hmuon_maxdr.GetYaxis().SetTitle("Events / 0.05")
hmuon_maxdr.SetLineColor(2)
h1d.append(hmuon_maxdr)

hmuon_normchi2 = ROOT.TH1D("hmuon_normchi2","",50,0.0,5.0)
hmuon_normchi2.GetXaxis().SetTitle("Muon #chi^{2}/ndof")
hmuon_normchi2.GetYaxis().SetTitle("Events")
hmuon_normchi2.SetLineColor(2)
h1d.append(hmuon_normchi2)

mtypes = ["G & T","G & !T","!G & T","SA & !G & !T","!SA & !G & !T"]
hmuon_type = ROOT.TH1D("hmuon_type","",len(mtypes),0,len(mtypes))
hmuon_type.GetXaxis().SetTitle("Muon type")
hmuon_type.GetYaxis().SetTitle("Events")
hmuon_type.SetLineColor(2)
for b in range(1,hmuon_type.GetNbinsX()+1):
    hmuon_type.GetXaxis().SetBinLabel(b,mtypes[b-1])
h1d.append(hmuon_type)

hmuon_ecaliso = ROOT.TH1D("hmuon_ecaliso","",40,0,20)
hmuon_ecaliso.GetXaxis().SetTitle("Muon ECAL isolation [GeV]")
hmuon_ecaliso.GetYaxis().SetTitle("Events / 0.5 GeV")
hmuon_ecaliso.SetLineColor(2)
h1d.append(hmuon_ecaliso)

hmuon_ecalreliso = ROOT.TH1D("hmuon_ecalreliso","",40,0,2.0)
hmuon_ecalreliso.GetXaxis().SetTitle("Muon ECAL isolation / p_{T}")
hmuon_ecalreliso.GetYaxis().SetTitle("Events / 0.05")
hmuon_ecalreliso.SetLineColor(2)
h1d.append(hmuon_ecalreliso)

hmuon_hcaliso = ROOT.TH1D("hmuon_hcaliso","",40,0,20)
hmuon_hcaliso.GetXaxis().SetTitle("Muon HCAL isolation [GeV]")
hmuon_hcaliso.GetYaxis().SetTitle("Events / 0.5 GeV")
hmuon_hcaliso.SetLineColor(2)
h1d.append(hmuon_hcaliso)

hmuon_hcalreliso = ROOT.TH1D("hmuon_hcalreliso","",40,0,2.0)
hmuon_hcalreliso.GetXaxis().SetTitle("Muon HCAL isolation / p_{T}")
hmuon_hcalreliso.GetYaxis().SetTitle("Events / 0.05")
hmuon_hcalreliso.SetLineColor(2)
h1d.append(hmuon_hcalreliso)

hmuon_trackiso = ROOT.TH1D("hmuon_trackiso","",20,0,10)
hmuon_trackiso.GetXaxis().SetTitle("Muon track isolation [GeV]")
hmuon_trackiso.GetYaxis().SetTitle("Events / 0.5 GeV")
hmuon_trackiso.SetLineColor(2)
h1d.append(hmuon_trackiso)

hmuon_trackreliso = ROOT.TH1D("hmuon_trackreliso","",20,0,1.0)
hmuon_trackreliso.GetXaxis().SetTitle("Muon track isolation / p_{T}")
hmuon_trackreliso.GetYaxis().SetTitle("Events / 0.05")
hmuon_trackreliso.SetLineColor(2)
h1d.append(hmuon_trackreliso)

hmuon_mindrjet = ROOT.TH1D("hmuon_mindrjet","",100,0,5.0)
hmuon_mindrjet.GetXaxis().SetTitle("min #DeltaR(#mu, PF jet)")
hmuon_mindrjet.GetYaxis().SetTitle("Events / 0.05")
hmuon_mindrjet.SetLineColor(2)
h1d.append(hmuon_mindrjet)

hmuon_dxy = ROOT.TH1D("hmuon_dxy","",100,0,5.0)
hmuon_dxy.GetXaxis().SetTitle("Muon |d_{xy}| [cm]")
hmuon_dxy.GetYaxis().SetTitle("Events / 0.05 cm")
hmuon_dxy.SetLineColor(2)
h1d.append(hmuon_dxy)

hmuon_dxysig = ROOT.TH1D("hmuon_dxysig","",100,0,50)
hmuon_dxysig.GetXaxis().SetTitle("Muon |d_{xy}|/#sigma_{xy}")
hmuon_dxysig.GetYaxis().SetTitle("Events / 0.5")
hmuon_dxysig.SetLineColor(2)
h1d.append(hmuon_dxysig)

hmuon_dz = ROOT.TH1D("hmuon_dz","",300,0,30)
hmuon_dz.GetXaxis().SetTitle("Muon |d_{z}| [cm]")
hmuon_dz.GetYaxis().SetTitle("Events / 0.1 cm")
hmuon_dz.SetLineColor(2)
h1d.append(hmuon_dz)

hmuon_dzsig = ROOT.TH1D("hmuon_dzsig","",100,0,50)
hmuon_dzsig.GetXaxis().SetTitle("Muon |d_{z}|/#sigma_{z}")
hmuon_dzsig.GetYaxis().SetTitle("Events / 0.5")
hmuon_dzsig.SetLineColor(2)
h1d.append(hmuon_dzsig)

hmuon_nsahits = ROOT.TH1D("hmuon_nsahits","",50,0,50)
hmuon_nsahits.GetXaxis().SetTitle("Number of valid SA muon hits")
hmuon_nsahits.GetYaxis().SetTitle("Events")
hmuon_nsahits.SetLineColor(2)
h1d.append(hmuon_nsahits)

hmuon_nsamatchedstats = ROOT.TH1D("hmuon_nsamatchedstats","",10,0,10)
hmuon_nsamatchedstats.GetXaxis().SetTitle("Number of SA muon matched stations")
hmuon_nsamatchedstats.GetYaxis().SetTitle("Events")
hmuon_nsamatchedstats.SetLineColor(2)
h1d.append(hmuon_nsamatchedstats)

hmuon_nmuhits = ROOT.TH1D("hmuon_nmuhits","",50,0,50)
hmuon_nmuhits.GetXaxis().SetTitle("Number of valid muon hits")
hmuon_nmuhits.GetYaxis().SetTitle("Events")
hmuon_nmuhits.SetLineColor(2)
h1d.append(hmuon_nmuhits)

hmuon_nmuchambs = ROOT.TH1D("hmuon_nmuchambs","",25,0,25)
hmuon_nmuchambs.GetXaxis().SetTitle("Number of muon chambers")
hmuon_nmuchambs.GetYaxis().SetTitle("Events")
hmuon_nmuchambs.SetLineColor(2)
h1d.append(hmuon_nmuchambs)

hmuon_nmuchambsCSCorDT = ROOT.TH1D("hmuon_nmuchambsCSCorDT","",20,0,20)
hmuon_nmuchambsCSCorDT.GetXaxis().SetTitle("Number of muon chambers (CSC or DT)")
hmuon_nmuchambsCSCorDT.GetYaxis().SetTitle("Events")
hmuon_nmuchambsCSCorDT.SetLineColor(2)
h1d.append(hmuon_nmuchambsCSCorDT)

hmuon_nmumatches = ROOT.TH1D("hmuon_nmumatches","",10,0,10)
hmuon_nmumatches.GetXaxis().SetTitle("Number of muon matches")
hmuon_nmumatches.GetYaxis().SetTitle("Events")
hmuon_nmumatches.SetLineColor(2)
h1d.append(hmuon_nmumatches)

hmuon_nmumatchedstats = ROOT.TH1D("hmuon_nmumatchedstats","",10,0,10)
hmuon_nmumatchedstats.GetXaxis().SetTitle("Number of muon matched stations")
hmuon_nmumatchedstats.GetYaxis().SetTitle("Events")
hmuon_nmumatchedstats.SetLineColor(2)
h1d.append(hmuon_nmumatchedstats)

hmuon_nmuexpmatchedstats = ROOT.TH1D("hmuon_nmuexpmatchedstats","",10,0,10)
hmuon_nmuexpmatchedstats.GetXaxis().SetTitle("Number of muon expected matched stations")
hmuon_nmuexpmatchedstats.GetYaxis().SetTitle("Events")
hmuon_nmuexpmatchedstats.SetLineColor(2)
h1d.append(hmuon_nmuexpmatchedstats)

hmuon_nmumatchedstatsmexp = ROOT.TH1D("hmuon_nmumatchedstatsmexp","",20,-10,10)
hmuon_nmumatchedstatsmexp.GetXaxis().SetTitle("Number of muon matched stations - expected")
hmuon_nmumatchedstatsmexp.GetYaxis().SetTitle("Events")
hmuon_nmumatchedstatsmexp.SetLineColor(2)
h1d.append(hmuon_nmumatchedstatsmexp)

hmuon_nmumatchedRPClayers = ROOT.TH1D("hmuon_nmumatchedRPClayers","",10,0,10)
hmuon_nmumatchedRPClayers.GetXaxis().SetTitle("Number of muon matched RPC layers")
hmuon_nmumatchedRPClayers.GetYaxis().SetTitle("Events")
hmuon_nmumatchedRPClayers.SetLineColor(2)
h1d.append(hmuon_nmumatchedRPClayers)

hmuon_npixelhits = ROOT.TH1D("hmuon_npixelhits","",10,0,10)
hmuon_npixelhits.GetXaxis().SetTitle("Number of pixel hits")
hmuon_npixelhits.GetYaxis().SetTitle("Events")
hmuon_npixelhits.SetLineColor(2)
h1d.append(hmuon_npixelhits)

hmuon_npixellayers = ROOT.TH1D("hmuon_npixellayers","",10,0,10)
hmuon_npixellayers.GetXaxis().SetTitle("Number of pixel layers")
hmuon_npixellayers.GetYaxis().SetTitle("Events")
hmuon_npixellayers.SetLineColor(2)
h1d.append(hmuon_npixellayers)

hmuon_nstriphits = ROOT.TH1D("hmuon_nstriphits","",30,0,30)
hmuon_nstriphits.GetXaxis().SetTitle("Number of strip hits")
hmuon_nstriphits.GetYaxis().SetTitle("Events")
hmuon_nstriphits.SetLineColor(2)
h1d.append(hmuon_nstriphits)

hmuon_ntrackerlayers = ROOT.TH1D("hmuon_ntrackerlayers","",25,0,25)
hmuon_ntrackerlayers.GetXaxis().SetTitle("Number of tracker layers")
hmuon_ntrackerlayers.GetYaxis().SetTitle("Events")
hmuon_ntrackerlayers.SetLineColor(2)
h1d.append(hmuon_ntrackerlayers)

if doPFIso:
    hmuon_mindrpfc = ROOT.TH1D("hmuon_mindrpfc","",100,0,0.1)
    hmuon_mindrpfc.GetXaxis().SetTitle("min #DeltaR(#mu, PF candidate)")
    hmuon_mindrpfc.GetYaxis().SetTitle("Events / 0.001")
    hmuon_mindrpfc.SetLineColor(2)
    h1d.append(hmuon_mindrpfc)

    hmuon_pfiso0p3chg = ROOT.TH1D("hmuon_pfiso0p3chg","",40,0,20)
    hmuon_pfiso0p3chg.GetXaxis().SetTitle("Muon PF charged isolation [GeV]")
    hmuon_pfiso0p3chg.GetYaxis().SetTitle("Events / 0.5 GeV")
    hmuon_pfiso0p3chg.SetLineColor(2)
    h1d.append(hmuon_pfiso0p3chg)

    hmuon_pfreliso0p3chg = ROOT.TH1D("hmuon_pfreliso0p3chg","",40,0,2.0)
    hmuon_pfreliso0p3chg.GetXaxis().SetTitle("Muon PF charged isolation / p_{T}")
    hmuon_pfreliso0p3chg.GetYaxis().SetTitle("Events / 0.05")
    hmuon_pfreliso0p3chg.SetLineColor(2)
    h1d.append(hmuon_pfreliso0p3chg)

    hmuon_pfiso0p3all = ROOT.TH1D("hmuon_pfiso0p3all","",40,0,20)
    hmuon_pfiso0p3all.GetXaxis().SetTitle("Muon PF all isolation (#delta#Beta) [GeV]")
    hmuon_pfiso0p3all.GetYaxis().SetTitle("Events / 0.5 GeV")
    hmuon_pfiso0p3all.SetLineColor(2)
    h1d.append(hmuon_pfiso0p3all)

    hmuon_pfreliso0p3all = ROOT.TH1D("hmuon_pfreliso0p3all","",40,0,2.0)
    hmuon_pfreliso0p3all.GetXaxis().SetTitle("Muon PF all isolation (#delta#Beta) / p_{T}")
    hmuon_pfreliso0p3all.GetYaxis().SetTitle("Events / 0.05")
    hmuon_pfreliso0p3all.SetLineColor(2)
    h1d.append(hmuon_pfreliso0p3all)

drt = [0.05,0.1,0.2,0.3,0.4,0.5,0.8,1.1,1.5,1e10]
drtLabels = []
for drtn in range(len(drt)-1):
    drtLabels.append("<%.2f"%drt[drtn])
drtLabels.append(">%.2f"%drt[len(drt)-2])
h_nmuons_vs_dr = ROOT.TH2D("h_nmuons_vs_dr","",len(drt),0,len(drt),10,0,10)
h_nmuons_vs_dr.GetXaxis().SetTitle("#DeltaR(#mu_{i}, #mu_{j})")
h_nmuons_vs_dr.GetYaxis().SetTitle("Number of other muons within (exclusive) cone")
h_nmuons_vs_dr.GetZaxis().SetTitle("Events")
for b in range(1,h_nmuons_vs_dr.GetNbinsX()+1):
    h_nmuons_vs_dr.GetXaxis().SetBinLabel(b,drtLabels[b-1])
h2d.append(h_nmuons_vs_dr)

###

for h in h1d:
    h.Sumw2(ROOT.kFALSE)
for h in h2d:
    h.Sumw2(ROOT.kFALSE)

###

print("Starting loop over %d events"%events.size())
for en,e in enumerate(events):
    ###if en>=1000:
    ###    break
    # Run & lumisection
    eventAux = e.eventAuxiliary()
    eid = eventAux.id()
    run  = int(eid.run())
    lumi = int(eid.luminosityBlock())
    evtn = int(eid.event())
    # JSON
    if isData and applyJSON:
        if not jsonUtils.isgoodrun(run,lumi,jsonruns):
            continue
    # L1t seed
    e.getByLabel(l1slabel,l1shandle)
    l1sStrings = l1shandle.product()
    e.getByLabel(l1plabel,l1phandle)
    l1sPrescales = l1phandle.product()
    passl1t = False
    l1seeds = []
    l1seeds.append("L1_DoubleMu_15_7")
    l1seeds.append("L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7")
    l1seeds.append("L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18")
    l1seeds.append("L1_DoubleMu4_SQ_OS_dR_Max1p2")
    l1seeds.append("L1_DoubleMu4p5_SQ_OS_dR_Max1p2")
    for l1sn,l1s in enumerate(l1sStrings):
        if l1sStrings[l1sn] in l1seeds and l1sPrescales[l1sn]==1:
            passl1t=True
            break
    if not passl1t:
        continue
    # HLT path
    e.getByLabel(hltlabel,hlthandle)
    hltString = hlthandle.product()[0]
    if "DoubleMu3_noVtx" not in hltString:
        continue

    # Muons
    nmuonsass=0
    nmuonssel=0
    e.getByLabel(mlabel,mhandle)
    muons = mhandle.product()

    # Jets
    e.getByLabel(jlabel,jhandle)
    pfjets = jhandle.product()

    # PF candidates
    e.getByLabel(pflabel,pfhandle)
    pfcands = pfhandle.product()

    # Displaced vertices
    e.getByLabel(svlabel,svhandle)
    sverts = svhandle.product()
    dtype = [('index',int), ('x',float), ('y',float), ('z',float), ('xe',float), ('ye',float), ('ze',float),
             ('chi2',float), ('ndof',int), ('prob',float), ('isValid', bool)]
    svvalues = []
    for svn,sv in enumerate(sverts):
        if sv.isValidVtx():
            svvalues.append((svn, sv.x(), sv.y(), sv.z(), sv.xError(), sv.yError(), sv.zError(),
                             sv.chi2(), sv.ndof(), ROOT.TMath.Prob(sv.chi2(), sv.ndof()), sv.isValidVtx()))
    svs = np.array(svvalues,dtype=dtype)
    if len(svs)<1:
        continue
    # Sort displaced vertices according to chi2 probability
    np.sort(svs, order='prob')

    # Primary vertices
    e.getByLabel(pvlabel,pvhandle)
    pverts = pvhandle.product()
    pvx,pvy,pvz=0.0,0.0,0.0
    if len(pverts)>0 and pverts[0].isValidVtx():
        pvx = pverts[0].x()
        pvy = pverts[0].y()
        pvz = pverts[0].z()

    h_nsv.Fill(len(svs))
    svvalues_sel = []
    for v in svs:
        hsv_chi2ndof.Fill(v['chi2']/v['ndof'])
        hsv_chi2prob.Fill(v['prob'])
        hsv_xerr.Fill(v['xe'])
        hsv_yerr.Fill(v['ye'])
        hsv_zerr.Fill(v['ze'])
        lxy = ROOT.TMath.Sqrt((v['x']-pvx)*(v['x']-pvx)+(v['y']-pvy)*(v['y']-pvy))
        hsv_lxy.Fill(lxy)
        l3d = ROOT.TMath.Sqrt(lxy*lxy+(v['z']-pvz)*(v['z']-pvz))
        hsv_l3d.Fill(l3d)
        if v['chi2']/v['ndof']>5.0 or v['xe']>0.05 or v['ye']>0.05 or v['ze']>0.1:
            continue
        hsvsel_lxy.Fill(lxy)
        hsvsel_l3d.Fill(l3d)
        svvalues_sel.append((v['index'],v['x'],v['y'],v['z'],v['xe'],v['ye'],v['ze'],
                             v['chi2'],v['ndof'],v['prob'],v['isValid']))
    selsvs = np.array(svvalues_sel,dtype=dtype)
    h_nsvsel.Fill(len(selsvs))

    if len(selsvs)>1:
        mindx,mindy,mindxy,mindz,mind3d = 1e6,1e6,1e6,1e6,1e6
        maxdx,maxdy,maxdxy,maxdz,maxd3d = 0.0,0.0,0.0,0.0,0.0
        for vn in range(0,len(selsvs)):
            for vvn in range(vn+1,len(selsvs)):
                thisdx = (selsvs['x'][vn]-selsvs['x'][vvn])*(selsvs['x'][vn]-selsvs['x'][vvn])
                thisdy = (selsvs['y'][vn]-selsvs['y'][vvn])*(selsvs['y'][vn]-selsvs['y'][vvn])
                thisdxy = ROOT.TMath.Sqrt(thisdx+thisdy)
                thisdz = (selsvs['z'][vn]-selsvs['z'][vvn])*(selsvs['z'][vn]-selsvs['z'][vvn])
                thisd3d = ROOT.TMath.Sqrt(thisdx+thisdy+thisdz)
                thisdx = ROOT.TMath.Sqrt(thisdx)
                thisdy = ROOT.TMath.Sqrt(thisdy)
                thisdz = ROOT.TMath.Sqrt(thisdz)
                if thisdx < mindx:
                    mindx = thisdx
                if thisdx > maxdx:
                    maxdx = thisdx
                if thisdy < mindy:
                    mindy = thisdy
                if thisdy > maxdy:
                    maxdy = thisdy
                if thisdxy < mindxy:
                    mindxy = thisdxy
                if thisdxy > maxdxy:
                    maxdxy = thisdxy
                if thisdz < mindz:
                    mindz = thisdz
                if thisdz > maxdz:
                    maxdz = thisdz
                if thisd3d < mind3d:
                    mind3d = thisd3d
                if thisd3d > maxd3d:
                    maxd3d = thisd3d
        hsvsel_mindx.Fill(mindx)
        hsvsel_maxdx.Fill(maxdx)
        hsvsel_mindy.Fill(mindy)
        hsvsel_maxdy.Fill(maxdy)
        hsvsel_mindxy.Fill(mindxy)
        hsvsel_maxdxy.Fill(maxdxy)
        hsvsel_mindz.Fill(mindz)
        hsvsel_maxdz.Fill(maxdz)
        hsvsel_mind3d.Fill(mind3d)
        hsvsel_maxd3d.Fill(maxd3d)

    vtxidxs = []
    nass    = []
    muchi2ndof = []
    muvec  = []  
    much   = []
    mutype = [] #mtypes = ["G & T","G & !T","!G & T","SA & !G & !T","!SA & !G & !T"]
    mueiso = []
    muhiso = []
    mutiso = []
    mudxy  = []
    mudxye = []
    mudz   = []
    mudze  = []
    musahits = []
    musamatchedstats = []
    mumuhits = []
    mumuchambs = []
    mumucscdt = []
    mumumatch = []
    mumumatchedstats = []
    mumuexpmatchedstats = []
    mumumatchedrpc = []
    mupixhits = []
    mustriphits = []
    mupixlayers = []
    mutrklayers = []
    mupfisochg = []
    mupfisoall = []
    mumindrpfc = []
    for mn,m in enumerate(muons):
        tvtxidxs = np.array(m.vtxIndx())
        tnass = 0
        for v in tvtxidxs:
            if v in selsvs['index']:
                tnass = tnass+1
        if tnass<1:
            continue
        nmuonsass = nmuonsass+1
        if not (m.pt()>3.0 and abs(m.eta())<2.4):
            continue
        nmuonssel = nmuonssel+1
        vtxidxs.append(tvtxidxs)
        nass.append(tnass)
        tvec = ROOT.TLorentzVector()
        tvec.SetPtEtaPhiM(m.pt(), m.eta(), m.phi(), MUON_MASS)
        muvec.append(tvec)
        much.append(m.charge())
        muchi2ndof.append(m.normalizedChi2())
        if isGlobalMuon(m.type()) and isTrackerMuon(m.type()):
            mutype.append(0.5)
        elif isGlobalMuon(m.type()) and not isTrackerMuon(m.type()):
            mutype.append(1.5)
        elif not isGlobalMuon(m.type()) and isTrackerMuon(m.type()):
            mutype.append(2.5)
        elif not isGlobalMuon(m.type()) and not isTrackerMuon(m.type()) and isStandAloneMuon(m.type()):
            mutype.append(3.5)
        elif not isGlobalMuon(m.type()) and not isTrackerMuon(m.type()) and not isStandAloneMuon(m.type()):
            mutype.append(4.5)
        mueiso.append(m.ecalIso())
        muhiso.append(m.hcalIso())
        mutiso.append(m.trackIso())
        mudxy.append(m.trk_dxy())
        mudxye.append(m.trk_dxyError())
        mudz.append(m.trk_dz())
        mudze.append(m.trk_dzError())
        musahits.append(m.nValidStandAloneMuonHits())
        musamatchedstats.append(m.nStandAloneMuonMatchedStations())
        mumuhits.append(m.nValidRecoMuonHits())
        mumuchambs.append(m.nRecoMuonChambers())
        mumucscdt.append(m.nRecoMuonChambersCSCorDT())
        mumumatch.append(m.nRecoMuonMatches())
        mumumatchedstats.append(m.nRecoMuonMatchedStations())
        mumuexpmatchedstats.append(m.nRecoMuonExpectedMatchedStations())
        mumumatchedrpc.append(m.nRecoMuonMatchedRPCLayers())
        mupixhits.append(m.nValidPixelHits())
        mustriphits.append(m.nValidStripHits())
        mupixlayers.append(m.nPixelLayersWithMeasurement())
        mutrklayers.append(m.nTrackerLayersWithMeasurement())
        if doPFIso:
            chiso,nhiso,phiso,puiso,mindrpfcmu = isolationUtils.getPFIsolation(muvec[nmuonssel-1],pfcands)
            mupfisochg.append(chiso)
            mupfisoall.append(chiso+max(0.0,nhiso+phiso-0.5*puiso))
            mumindrpfc.append(mindrpfcmu)
    h_nmuonsass.Fill(nmuonsass)
    h_nmuonssel.Fill(nmuonssel)
    if nmuonssel<2:
        continue

    mindr = 1e6
    maxdr = 0.0
    ndr = []
    for drtn in range(len(drt)):
        ndr.append(0)
    for m in range(0,nmuonssel):
        hmuon_nSVs.Fill(nass[m])
        hmuon_pt.Fill(muvec[m].Pt())
        hmuon_eta.Fill(muvec[m].Eta())
        hmuon_phi.Fill(muvec[m].Phi())
        hmuon_ch.Fill(much[m])
        hmuon_normchi2.Fill(muchi2ndof[m])
        hmuon_type.Fill(mutype[m])
        for mm in range(m+1,nmuonssel):
            tdr = muvec[m].DeltaR(muvec[mm])
            drtn = 0
            while drtn < len(drt):
                if tdr < drt[drtn]:
                    ndr[drtn] = ndr[drtn]+1
                drtn = drtn+1
            if tdr < mindr:
                mindr = tdr
            if tdr > maxdr:
                maxdr = tdr
        mindrj= 1e6
        for j in pfjets:
            if j.pt()<20.0 or abs(j.eta())>2.5:
                continue
            jvec = ROOT.TLorentzVector()
            jvec.SetPtEtaPhiM(j.pt(),j.eta(),j.phi(),j.m())
            tdrj = muvec[m].DeltaR(jvec)
            if tdrj < mindrj:
                mindrj = tdrj
        hmuon_ecaliso.Fill(mueiso[m])
        hmuon_ecalreliso.Fill(mueiso[m]/muvec[m].Pt())
        hmuon_hcaliso.Fill(muhiso[m])
        hmuon_hcalreliso.Fill(muhiso[m]/muvec[m].Pt())
        hmuon_trackiso.Fill(mutiso[m])
        hmuon_trackreliso.Fill(mutiso[m]/muvec[m].Pt())
        hmuon_dxy.Fill(abs(mudxy[m]))
        hmuon_dxysig.Fill(abs(mudxy[m]/mudxye[m]))
        hmuon_dz.Fill(abs(mudz[m]))
        hmuon_dzsig.Fill(abs(mudz[m]/mudze[m]))
        hmuon_nsahits.Fill(musahits[m])
        hmuon_nsamatchedstats.Fill(musamatchedstats[m])
        hmuon_nmuhits.Fill(mumuhits[m])
        hmuon_nmuchambs.Fill(mumuchambs[m])
        hmuon_nmuchambsCSCorDT.Fill(mumucscdt[m])
        hmuon_nmumatches.Fill(mumumatch[m])
        hmuon_nmumatchedstats.Fill(mumumatchedstats[m])
        hmuon_nmuexpmatchedstats.Fill(mumuexpmatchedstats[m])
        hmuon_nmumatchedstatsmexp.Fill(mumumatchedstats[m]-mumuexpmatchedstats[m])
        hmuon_nmumatchedRPClayers.Fill(mumumatchedrpc[m])
        hmuon_npixelhits.Fill(mupixhits[m])
        hmuon_npixellayers.Fill(mupixlayers[m])
        hmuon_nstriphits.Fill(mustriphits[m])
        hmuon_ntrackerlayers.Fill(mutrklayers[m])
        hmuon_mindrjet.Fill(mindrj)
        if doPFIso:
            hmuon_mindrpfc.Fill(mumindrpfc[m])
            hmuon_pfiso0p3chg.Fill(mupfisochg[m])
            hmuon_pfreliso0p3chg.Fill(mupfisochg[m]/muvec[m].Pt())
            hmuon_pfiso0p3all.Fill(mupfisoall[m])
            hmuon_pfreliso0p3all.Fill(mupfisoall[m]/muvec[m].Pt())
    hmuon_mindr.Fill(mindr)
    hmuon_maxdr.Fill(maxdr)
    for drtn in range(len(drt)):
        h_nmuons_vs_dr.Fill(drtn,ndr[drtn])

doPlots = True
if len(sys.argv)>2:
    doPlots = False
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
    if not doPlots:
        break
    if "_type" not in h.GetName():
        plotUtils.PutUnderflowInFirstBin(h)
        plotUtils.PutOverflowInLastBin(h)
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
    if "hsv" in h.GetName():
        ytitle = ytitle.replace("Events", "Number of SVs")
        ytitle = ytitle.replace("events", "SVs")
    if "hmuon" in h.GetName():
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
    if not doPlots:
        break
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

foname = "%s/histograms_all.root"%outdir
if len(sys.argv)>2:
    if len(sys.argv)>3:
        foname = "%s/histograms_%s.root"%(outdir,sys.argv[2])
    else:
        foname = "%s/histograms_file%s.root"%(outdir,sys.argv[2])
fout = ROOT.TFile(foname,"RECREATE")
fout.cd()
for h in h1d:
    h.Write()
for h in h2d:
    h.Write()
fout.Close()
