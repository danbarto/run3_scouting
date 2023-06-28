#include <algorithm>
#include <filesystem>
#include <numeric>
#include <string>
#include <tuple>
namespace fs = std::filesystem;

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
using namespace fwlite;

#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"

#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"

#include "tools/dorky.h"
#include "tools/goodrun.h"
#include "tools/tqdm.h"


template<class T>
const T getObject(Event &ev, const char* prodLabel, const char* outLabel = "") {
  Handle<T> obj;
  obj.getByLabel(ev,prodLabel,outLabel);
  return *obj.product();
}

template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(const std::vector<T>& vec, Compare compare) {
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
  return p;
}

template <typename T>
void apply_permutation_in_place(std::vector<T>& vec, const std::vector<std::size_t>& p) {
  std::vector<bool> done(vec.size());
  for (std::size_t i=0; i<vec.size(); ++i) {
    if (done[i])
      continue;
    done[i] = true;
    std::size_t prev_j = i;
    std::size_t j = p[i];
    while (i != j) {
      std::swap(vec[prev_j], vec[j]);
      done[j] = true;
      prev_j = j;
      j = p[j];
    }
  }
}


float getMiniIsoDR(TLorentzVector muVec, float minDR=0.05, float maxDR=0.2,float kt=10.0) {
  return std::min(maxDR, std::max(minDR, kt/(float)muVec.Pt()));
}

std::tuple<float,float,float,float,float> getPFIsolation(TLorentzVector muVec, std::vector<Run3ScoutingParticle> pfs, float maxDR=0.3, float vetoDR=0.01, float maxdz=0.1, bool doMini=false, bool doMindrMuPF=true) {
  float maxdr, mindrMuPF=1e6;
  float chIso=0, nhIso=0, phIso=0, puIso=0;
  if (doMini)
    maxdr = getMiniIsoDR(muVec);
  for (auto pf : pfs) {
    TLorentzVector pfvec;
    pfvec.SetPtEtaPhiM(pf.pt(),pf.eta(),pf.phi(),0.0);
    float dr = muVec.DeltaR(pfvec);
    if (dr<mindrMuPF)
      mindrMuPF = dr;
    if (dr<vetoDR || dr>maxDR)
      continue;
    bool fromPV = !(abs(pf.vertex())>0 || pf.dz()>maxdz);
    bool isCh=false, isNh=false, isPh=false, isPU=false;
    if (fabs(pf.pdgId())==211) {
      if (fromPV) {
        isCh = true;
        chIso = chIso + pf.pt();
      }
      else {
        isPU = true;
        puIso = puIso + pf.pt();
      }
    }
    if (fabs(pf.pdgId())==130) {
      isNh  = true;
      nhIso = nhIso + pf.pt();
    }
    if (fabs(pf.pdgId())==22) {
      isPh  = true;
      phIso = phIso + pf.pt();
    }
  }
  if (doMindrMuPF)
    return {chIso,nhIso,phIso,puIso,mindrMuPF};
  else
    return {chIso,nhIso,phIso,puIso,1e6};
}


struct SV {
  std::vector<unsigned int> index, ndof;
  std::vector<float> x, y, z;
  std::vector<float> xe, ye, ze;
  std::vector<float> chi2, prob, chi2Ndof;
  std::vector<float> lxy, l3d;
  std::vector<float> mindx, mindy, mindz, mindxy, mind3d;
  std::vector<float> maxdx, maxdy, maxdz, maxdxy, maxd3d;
  std::vector<bool> isValid;

  void clear() {
    index.clear(); ndof.clear();
    x.clear(); y.clear(); z.clear();
    xe.clear(); ye.clear(); ze.clear();
    chi2.clear(); prob.clear(); chi2Ndof.clear();
    lxy.clear(); l3d.clear();
    mindx.clear(); mindy.clear(); mindz.clear(); mindxy.clear(); mind3d.clear();
    maxdx.clear(); maxdy.clear(); maxdz.clear(); maxdxy.clear(); maxd3d.clear();
    isValid.clear();
  }

  void sort() {
    auto comp = sort_permutation(prob, [](float const& a, float const& b){ return a < b; }); // Define permutations based on the prob vector
    // Apply the permutation to all vectors
    apply_permutation_in_place(index, comp);
    apply_permutation_in_place(ndof, comp);
    apply_permutation_in_place(x, comp);
    apply_permutation_in_place(y, comp);
    apply_permutation_in_place(z, comp);
    apply_permutation_in_place(xe, comp);
    apply_permutation_in_place(ye, comp);
    apply_permutation_in_place(ze, comp);
    apply_permutation_in_place(chi2, comp);
    apply_permutation_in_place(prob, comp);
    apply_permutation_in_place(chi2Ndof, comp);
    apply_permutation_in_place(lxy, comp);
    apply_permutation_in_place(l3d, comp);
    apply_permutation_in_place(mindx, comp);
    apply_permutation_in_place(mindy, comp);
    apply_permutation_in_place(mindz, comp);
    apply_permutation_in_place(mindxy, comp);
    apply_permutation_in_place(mind3d, comp);
    apply_permutation_in_place(maxdx, comp);
    apply_permutation_in_place(maxdy, comp);
    apply_permutation_in_place(maxdz, comp);
    apply_permutation_in_place(maxdxy, comp);
    apply_permutation_in_place(maxd3d, comp);
    apply_permutation_in_place(isValid, comp);
  }
};

float MUON_MASS = 0.10566;

bool isGlobalMuon(unsigned int type) {
  return type & (1<<1);
}

bool isTrackerMuon(unsigned int type) {
  return type & (1<<2);
}

bool isStandAloneMuon(unsigned int type) {
  return type & (1<<3);
}

struct Muon {
  std::vector<std::vector<int>> vtxIdxs;
  std::vector<unsigned int> saHits, saMatchedStats;
  std::vector<unsigned int> muHits, muChambs, muCSCDT, muMatch, muMatchedStats, muExpMatchedStats, muMatchedRPC;
  std::vector<unsigned int> pixHits, stripHits;
  std::vector<unsigned int> pixLayers, trkLayers;
  std::vector<int> ch;
  std::vector<float> pt, eta, phi;
  std::vector<float> chi2ndof;
  std::vector<float> ecalIso, hcalIso, trackIso;
  std::vector<float> ecalRelIso, hcalRelIso, trackRelIso;
  std::vector<float> dxy, dxye, dz, dze;
  std::vector<float> dxysig, dzsig;
  std::vector<float> PFIsoChg, PFIsoAll;
  std::vector<float> PFRelIsoChg, PFRelIsoAll;
  std::vector<float> mindr, maxdr;
  std::vector<float> mindrPF, mindrJet;
  std::vector<TLorentzVector> vec;
  std::vector<TString> type; // types = ["G & T","G & !T","!G & T","SA & !G & !T","!SA & !G & !T"]

  void clear() {
  vtxIdxs.clear();
  saHits.clear(); saMatchedStats.clear();
  muHits.clear(); muChambs.clear(); muCSCDT.clear(); muMatch.clear(); muMatchedStats.clear(); muExpMatchedStats.clear(); muMatchedRPC.clear();
  pixHits.clear(); stripHits.clear();
  pixLayers.clear(); trkLayers.clear();
  ch.clear();
  pt.clear(); eta.clear(); phi.clear();
  chi2ndof.clear();
  ecalIso.clear(); hcalIso.clear(); trackIso.clear();
  ecalRelIso.clear(); hcalRelIso.clear(); trackRelIso.clear();
  dxy.clear(); dxye.clear(); dz.clear(); dze.clear();
  dxysig.clear(); dzsig.clear();
  PFIsoChg.clear(); PFIsoAll.clear();
  PFRelIsoChg.clear(); PFRelIsoAll.clear();
  mindr.clear(); maxdr.clear();
  mindrPF.clear(); mindrJet.clear();
  vec.clear();
  type.clear();
  }
};


std::vector<std::string> selL1Seeds = {"L1_DoubleMu_15_7",
                                       "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7",
                                       "L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18",
                                       "L1_DoubleMu4_SQ_OS_dR_Max1p2",
                                       "L1_DoubleMu4p5_SQ_OS_dR_Max1p2"};


void run3ScoutingLooper(std::vector<TString> inputFiles, TString year, TString process, const char* outdir="temp_data") {
  // Output folders and files
  fs::create_directory(outdir);
  fs::permissions(outdir,fs::perms::owner_all | fs::perms::group_read | fs::perms::group_exec | fs::perms::others_read | fs::perms::others_exec);
  TFile* fout = new TFile(TString(outdir)+"/output_"+process+"_"+year+".root", "RECREATE");
  TTree* tout = new TTree("tout","Run3ScoutingTree");

  // Branch variables
  float PV_x, PV_y, PV_z;
  SV SVs;
  int nMuonAssoc;
  Muon Muons;

  // Branch definition
  tout->Branch("PV_x", &PV_x);
  tout->Branch("PV_y", &PV_y);
  tout->Branch("PV_z", &PV_z);

  tout->Branch("SV_index", &SVs.index);
  tout->Branch("SV_ndof", &SVs.ndof);
  tout->Branch("SV_x", &SVs.x);
  tout->Branch("SV_y", &SVs.y);
  tout->Branch("SV_z", &SVs.z);
  tout->Branch("SV_xe", &SVs.xe);
  tout->Branch("SV_ye", &SVs.ye);
  tout->Branch("SV_ze", &SVs.ze);
  tout->Branch("SV_chi2", &SVs.chi2);
  tout->Branch("SV_prob", &SVs.prob);
  tout->Branch("SV_chi2OverNdof", &SVs.chi2Ndof);
  tout->Branch("SV_lxy", &SVs.lxy);
  tout->Branch("SV_l3d", &SVs.l3d);
  tout->Branch("SV_mindx", &SVs.mindx);
  tout->Branch("SV_mindy", &SVs.mindy);
  tout->Branch("SV_mindz", &SVs.mindz);
  tout->Branch("SV_maxdx", &SVs.maxdx);
  tout->Branch("SV_maxdy", &SVs.maxdy);
  tout->Branch("SV_maxdz", &SVs.maxdz);
  tout->Branch("SV_isValid", &SVs.isValid);

  tout->Branch("nMuonAssoc", &nMuonAssoc);
  tout->Branch("Muon_vtxIdx", &Muons.vtxIdxs);
  tout->Branch("Muon_saHits", &Muons.saHits);
  tout->Branch("Muon_saMatchedStats", &Muons.saMatchedStats);
  tout->Branch("Muon_muHits", &Muons.muHits);
  tout->Branch("Muon_muChambs", &Muons.muChambs);
  tout->Branch("Muon_muCSCDT", &Muons.muCSCDT);
  tout->Branch("Muon_muMatch", &Muons.muMatch);
  tout->Branch("Muon_muMatchedStats", &Muons.muMatchedStats);
  tout->Branch("Muon_muExpMatchedStats", &Muons.muExpMatchedStats);
  tout->Branch("Muon_muMatchedRPC", &Muons.muMatchedRPC);
  tout->Branch("Muon_pixHits", &Muons.pixHits);
  tout->Branch("Muon_stripHits", &Muons.stripHits);
  tout->Branch("Muon_pixLayers", &Muons.pixLayers);
  tout->Branch("Muon_trkLayers", &Muons.trkLayers);
  tout->Branch("Muon_pt", &Muons.pt);
  tout->Branch("Muon_eta", &Muons.eta);
  tout->Branch("Muon_phi", &Muons.phi);
  tout->Branch("Muon_ch", &Muons.ch);
  tout->Branch("Muon_chi2ndof", &Muons.chi2ndof);
  tout->Branch("Muon_ecalIso", &Muons.ecalIso);
  tout->Branch("Muon_hcalIso", &Muons.hcalIso);
  tout->Branch("Muon_trackIso", &Muons.trackIso);
  tout->Branch("Muon_ecalRelIso", &Muons.ecalRelIso);
  tout->Branch("Muon_hcalRelIso", &Muons.hcalRelIso);
  tout->Branch("Muon_trackRelIso", &Muons.trackRelIso);
  tout->Branch("Muon_dxy", &Muons.dxy);
  tout->Branch("Muon_dxye", &Muons.dxye);
  tout->Branch("Muon_dz", &Muons.dz);
  tout->Branch("Muon_dze", &Muons.dze);
  tout->Branch("Muon_dxysig", &Muons.dxysig);
  tout->Branch("Muon_dzsig", &Muons.dzsig);
  tout->Branch("Muon_PFIsoChg", &Muons.PFIsoChg);
  tout->Branch("Muon_PFIsoAll", &Muons.PFIsoAll);
  tout->Branch("Muon_PFRelIsoChg", &Muons.PFRelIsoChg);
  tout->Branch("Muon_PFRelIsoAll", &Muons.PFRelIsoAll);
  tout->Branch("Muon_mindr", &Muons.mindr);
  tout->Branch("Muon_maxdr", &Muons.maxdr);
  tout->Branch("Muon_mindrJet", &Muons.mindrJet);
  tout->Branch("Muon_mindrPF", &Muons.mindrPF);
  tout->Branch("Muon_vec", &Muons.vec);
  tout->Branch("Muon_type", &Muons.type);
  

  // Event setup
  tqdm bar;

  bool isMC = true;
  if (process == "Data")
    isMC = false;

  if ( !isMC ) {
    if ( year == "2022" )
      set_goodrun_file_json("../data/Cert_Collisions2022_355100_362760_Golden.json");
  }

  
  // File loop
  unsigned int iFile = 1;
  for (auto inputFile : inputFiles) {
    std::cout << "File number: " << iFile << "\n";
    TFile file(inputFile);
    auto nEventsFile = ((TTree*)file.Get("Events"))->GetEntries();
    Event ev(&file);

    // Event loop
    unsigned int iEv = 0;
    for (ev.toBegin(); ! ev.atEnd(); ++ev) {
      iEv++;
      //if (iEv > 10)
      //  break;
      bar.progress(iEv, nEventsFile);

      auto evAux = ev.eventAuxiliary();
      auto eID = evAux.id();
      auto run  = eID.run();
      auto lumi = eID.luminosityBlock();
      auto evtn = eID.event();

      // JSON and duplicate removal
      if ( !isMC ) {
        if ( !(goodrun(run, lumi)) )
          continue;
        duplicate_removal::DorkyEventIdentifier id(run, evtn, lumi);
        if ( is_duplicate(id) )
          continue;
      }


      // L1 selection
      auto l1s = getObject<std::vector<std::string>>(ev, "triggerMaker", "l1name");
      auto l1Prescales = getObject<std::vector<double>>(ev, "triggerMaker", "l1prescale");

      bool passL1 = false;
      for (unsigned int iL1=0; iL1<l1s.size(); ++iL1) {
        for (auto selL1Seed : selL1Seeds) {
          if (l1s[iL1] == selL1Seed) {
            passL1 = true;
            break;
          }
        }
        if (passL1 && l1Prescales[iL1])
          break;
      }
      if (!passL1)
        continue;

      // HLT selection
      auto hlts = getObject<std::vector<std::string>>(ev, "triggerMaker", "hltname");
      for (auto hlt : hlts) { 
        if (!TString(hlt).Contains("DoubleMu3_noVtx"))
          continue;
      }


      // PV selection
      auto pvs = getObject<std::vector<Run3ScoutingVertex>>(ev, "hltScoutingPrimaryVertexPacker", "primaryVtx");
      PV_x = 0; PV_y = 0; PV_z = 0;

      if (pvs.size() > 0 && pvs[0].isValidVtx()) {
        PV_x = pvs[0].x();
        PV_y = pvs[0].y();
        PV_z = pvs[0].z();
      }

      // SV selection
      auto svs = getObject<std::vector<Run3ScoutingVertex>>(ev, "hltScoutingMuonPacker", "displacedVtx");
      SVs.clear();

      unsigned int nSVs = svs.size();
      if (nSVs < 1)
        continue;

      for (unsigned int iSV=0; iSV<nSVs; ++iSV) {
        auto sv = svs[iSV];
        if (!(sv.isValidVtx()))
          continue;

        SVs.index.push_back(iSV);
        SVs.ndof.push_back(sv.ndof());
        SVs.x.push_back(sv.x());
        SVs.y.push_back(sv.y());
        SVs.z.push_back(sv.z());
        SVs.xe.push_back(sv.xError());
        SVs.ye.push_back(sv.yError());
        SVs.ze.push_back(sv.zError());
        SVs.chi2.push_back(sv.chi2());
        SVs.prob.push_back(TMath::Prob(sv.chi2(), sv.ndof()));
        SVs.chi2Ndof.push_back(sv.chi2()/sv.ndof());
        float lxy = TMath::Sqrt((sv.x()-PV_x)*(sv.x()-PV_x)+(sv.y()-PV_y)*(sv.y()-PV_y));
        SVs.lxy.push_back(lxy);
        SVs.l3d.push_back(TMath::Sqrt(lxy*lxy+(sv.z()-PV_z)*(sv.z()-PV_z)));

        float mindx=1e6, mindy=1e6, mindz=1e6, mindxy=1e6, mind3d=1e6;
        float maxdx=-1, maxdy=-1, maxdz=-1, maxdxy=-1, maxd3d=-1;
        for (unsigned int jSV=0; jSV<nSVs; ++jSV) {
          if (jSV==iSV)
            continue;
          auto svOther = svs[jSV];
          if (!(svOther.isValidVtx()))
            continue;

          float dx = (sv.x()-svOther.x())*(sv.x()-svOther.x()); // Squared for the moment
          float dy = (sv.y()-svOther.y())*(sv.y()-svOther.y()); // Squared for the moment
          float dz = (sv.z()-svOther.z())*(sv.z()-svOther.z()); // Squared for the moment
          float dxy = TMath::Sqrt(dx+dy);
          float d3d = TMath::Sqrt(dx+dy+dz);
          dx = TMath::Sqrt(dx); // Back to the proper value
          dy = TMath::Sqrt(dy); // Back to the proper value
          dz = TMath::Sqrt(dz); // Back to the proper value

          if (dx<mindx) mindx = dx;
          if (dx>maxdx) maxdx = dx;
          if (dy<mindy) mindy = dy;
          if (dy>maxdy) maxdy = dy;
          if (dz<mindz) mindz = dz;
          if (dz>maxdz) maxdz = dz;
          if (dxy<mindxy) mindxy = dxy;
          if (dxy>maxdxy) maxdxy = dxy;
          if (d3d<mind3d) mind3d = d3d;
          if (d3d>maxd3d) maxd3d = d3d;
        }
        SVs.mindx.push_back(mindx);
        SVs.mindy.push_back(mindy);
        SVs.mindz.push_back(mindz);
        SVs.mindxy.push_back(mindxy);
        SVs.mind3d.push_back(mind3d);
        SVs.maxdx.push_back(maxdx);
        SVs.maxdy.push_back(maxdy);
        SVs.maxdz.push_back(maxdz);
        SVs.maxdxy.push_back(maxdxy);
        SVs.maxd3d.push_back(maxd3d);
      }
      SVs.sort();

      // Muon selection
      auto mus = getObject<std::vector<Run3ScoutingMuon>>(ev, "hltScoutingMuonPacker");
      auto jets = getObject<std::vector<Run3ScoutingPFJet>>(ev, "hltScoutingPFPacker");
      auto pfs = getObject<std::vector<Run3ScoutingParticle>>(ev, "hltScoutingPFPacker");
      Muons.clear();
      unsigned int nMus = mus.size();
      nMuonAssoc=0;

      for (unsigned int iMu=0; iMu<nMus; ++iMu) {
        auto mu = mus[iMu];
        std::vector<int> matchedAndSelVtxIdxs;
        for (auto matchedVtxIdx : mu.vtxIndx()) {
          for (auto selVtxIdx : SVs.index) {
            if (matchedVtxIdx==selVtxIdx)
              matchedAndSelVtxIdxs.push_back(matchedVtxIdx);
          }
        }
        if (matchedAndSelVtxIdxs.size() < 1) // Require muon to be associated to at least one of the selected SVs
          continue;
        nMuonAssoc++; // Count muon associated to at least one selected SV
        if (!(mu.pt()>3.0 && fabs(mu.eta())<2.4))
          continue;
        Muons.vtxIdxs.push_back(matchedAndSelVtxIdxs);
        Muons.saHits.push_back(mu.nValidStandAloneMuonHits());
        Muons.saMatchedStats.push_back(mu.nStandAloneMuonMatchedStations());
        Muons.muHits.push_back(mu.nValidRecoMuonHits());
        Muons.muChambs.push_back(mu.nRecoMuonChambers());
        Muons.muCSCDT.push_back(mu.nRecoMuonChambersCSCorDT());
        Muons.muMatch.push_back(mu.nRecoMuonMatches());
        Muons.muMatchedStats.push_back(mu.nRecoMuonMatchedStations());
        Muons.muExpMatchedStats.push_back(mu.nRecoMuonExpectedMatchedStations());
        Muons.muMatchedRPC.push_back(mu.nRecoMuonMatchedRPCLayers());
        Muons.pixHits.push_back(mu.nValidPixelHits());
        Muons.stripHits.push_back(mu.nValidStripHits());
        Muons.pixLayers.push_back(mu.nPixelLayersWithMeasurement());
        Muons.trkLayers.push_back(mu.nTrackerLayersWithMeasurement());
        Muons.pt.push_back(mu.pt());
        Muons.eta.push_back(mu.eta());
        Muons.phi.push_back(mu.phi());
        Muons.ch.push_back(mu.charge());
        Muons.chi2ndof.push_back(mu.normalizedChi2());
        Muons.ecalIso.push_back(mu.ecalIso());
        Muons.hcalIso.push_back(mu.hcalIso());
        Muons.trackIso.push_back(mu.trackIso());
        Muons.ecalRelIso.push_back(mu.ecalIso()/mu.pt());
        Muons.hcalRelIso.push_back(mu.hcalIso()/mu.pt());
        Muons.trackRelIso.push_back(mu.trackIso()/mu.pt());
        Muons.dxy.push_back(mu.trk_dxy());
        Muons.dxye.push_back(mu.trk_dxyError());
        Muons.dz.push_back(mu.trk_dz());
        Muons.dze.push_back(mu.trk_dzError());
        Muons.dxysig.push_back(mu.trk_dxy()/mu.trk_dxyError());
        Muons.dzsig.push_back(mu.trk_dz()/mu.trk_dzError());

        TLorentzVector muVec; muVec.SetPtEtaPhiM(mu.pt(), mu.eta(), mu.phi(), MUON_MASS);
        Muons.vec.push_back(muVec);

        if (isGlobalMuon(mu.type()) && isTrackerMuon(mu.type()))
          Muons.type.push_back(0.5); // "G & T"
        if (isGlobalMuon(mu.type()) && !isTrackerMuon(mu.type()))
          Muons.type.push_back(1.5); // "G & !T"
        if (!isGlobalMuon(mu.type()) && isTrackerMuon(mu.type()))
          Muons.type.push_back(2.5); // "!G & T"
        if (!isGlobalMuon(mu.type()) && !isTrackerMuon(mu.type()) && isStandAloneMuon(mu.type()))
          Muons.type.push_back(3.5); // "SA & !G & !T"
        if (!isGlobalMuon(mu.type()) && !isTrackerMuon(mu.type()) && !isStandAloneMuon(mu.type()))
          Muons.type.push_back(4.5); // "!SA & !G & !T"

        const auto pfIsos = getPFIsolation(muVec, pfs);
        const auto pfIsoChg = std::get<0>(pfIsos);
        const auto pfIsoAll = pfIsoChg + std::max(0.0, std::get<1>(pfIsos)+std::get<2>(pfIsos)-0.5*std::get<3>(pfIsos));
        Muons.PFIsoChg.push_back(pfIsoChg);
        Muons.PFIsoAll.push_back(pfIsoAll);;
        Muons.PFRelIsoChg.push_back(pfIsoChg/mu.pt());
        Muons.PFRelIsoAll.push_back(pfIsoAll/mu.pt());
        Muons.mindrPF.push_back(std::get<4>(pfIsos));

        float mindr=1e6;
        float maxdr=-1;
        for (unsigned int jMu=0; jMu<nMus; ++jMu) {
          if (jMu==iMu)
            continue;
          auto muOther = mus[jMu];
          bool matchedAndSelVtx = false;
          for (auto matchedVtxIdx : muOther.vtxIndx()) {
            for (auto selVtxIdx : SVs.index) {
              if (matchedVtxIdx==selVtxIdx) {
                matchedAndSelVtx = true;
                break;
              }
              if (matchedAndSelVtx)
                break;
            }
          }
          if (!matchedAndSelVtx)
            continue;
          if (!(muOther.pt()>3.0 && fabs(muOther.eta())<2.4))
            continue;
          TLorentzVector muVecOther; muVecOther.SetPtEtaPhiM(muOther.pt(), muOther.eta(), muOther.phi(), MUON_MASS);
          float dr = muVec.DeltaR(muVecOther);
          if (dr<mindr) mindr = dr;
          if (dr>maxdr) maxdr = dr;
        }
        Muons.mindr.push_back(mindr);
        Muons.maxdr.push_back(maxdr);

        float mindrJet=1e6;
        for (auto jet : jets) {
          if (!(jet.pt()>20.0 && fabs(jet.eta())<2.5))
            continue;
          TLorentzVector jetVec; jetVec.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.m());
          float dr = muVec.DeltaR(jetVec);
          if (dr<mindrJet) mindrJet = dr;
        }
        Muons.mindrJet.push_back(mindrJet);
      }
      if (Muons.pt.size() < 2)
        continue;

      tout->Fill();
    }
    iFile++;
    std::cout<<"\n\n";
  }

  bar.finish();
  fout->Write();
  fout->Close();

  return;
}

