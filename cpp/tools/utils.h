#include <algorithm>
#include <numeric>

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

//import ROOT
//import os,sys
//from DataFormats.FWLite import Events, Handle
//
//def getMiniIsoDR(mvec,minDR=0.05,maxDR=0.2,kt=10.0):
//    return min(maxDR, max(minDR, kt / mvec.Pt()))
//
//def getPFIsolation(mvec,pfcands,maxDR=0.3,vetoDR=0.01,maxdz=0.1,domini=False,domindrmpf=True):
//    mindrmpf=1e6
//    chiso,nhiso,phiso,puiso=0.0,0.0,0.0,0.0
//    if domini:
//        maxDR = getMiniIsoDR(mvec)
//    for p in pfcands:
//        pfvec = ROOT.TLorentzVector()
//        pfvec.SetPtEtaPhiM(p.pt(),p.eta(),p.phi(),0.0)
//        dR = mvec.DeltaR(pfvec)
//        if dR<mindrmpf:
//            mindrmpf = dR
//        if dR<vetoDR or dR>maxDR:
//            continue
//        frompv = not (abs(p.vertex())>0 or p.dz()>maxdz)
//        isch,isnh,isph,ispu=False,False,False,False
//        if abs(p.pdgId())==211:
//            if frompv:
//                isch  = True
//                chiso = chiso+p.pt()
//            else:
//                ispu=True
//                puiso = puiso+p.pt()
//        if abs(p.pdgId())==130:
//            isnh  = True
//            nhiso = nhiso+p.pt()
//        if abs(p.pdgId())==22:
//            isph  = True
//            phiso = phiso+p.pt()
//    if domindrmpf:
//        return chiso,nhiso,phiso,puiso,mindrmpf
//    else:
//        return chiso,nhiso,phiso,puiso
//
