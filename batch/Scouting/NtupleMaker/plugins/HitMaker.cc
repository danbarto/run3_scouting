#include "Scouting/NtupleMaker/plugins/HitMaker.h" 
#include "TMath.h"

using namespace edm;
using namespace std;

namespace {

  Surface::RotationType rotation(const GlobalVector& zDir) {
    GlobalVector zAxis = zDir.unit();
    GlobalVector yAxis(zAxis.y(), -zAxis.x(), 0);
    GlobalVector xAxis = yAxis.cross(zAxis);
    return Surface::RotationType(xAxis, yAxis, zAxis);
  }


}

HitMaker::HitMaker(const edm::ParameterSet& iConfig):
magFieldToken_(esConsumes()),
trackingGeometryToken_(esConsumes()),
measurementTrackerToken_(esConsumes()),
propagatorToken_(esConsumes(edm::ESInputTag("", "PropagatorWithMaterial"))),
geomToken_(esConsumes())
{
    muonToken_ = consumes<Run3ScoutingMuonCollection>(iConfig.getParameter<InputTag>("muonInputTag"));
    dvToken_ = consumes<Run3ScoutingVertexCollection>(iConfig.getParameter<InputTag>("dvInputTag"));
    measurementTrackerEventToken_ = consumes<MeasurementTrackerEvent>(iConfig.getParameter<InputTag>("measurementTrackerEventInputTag"));

    produces<vector<vector<bool> > >("isbarrel").setBranchAlias("Muon_hit_barrel");
    produces<vector<vector<bool> > >("ispixel").setBranchAlias("Muon_hit_pixel");
    produces<vector<vector<bool> > >("isactive").setBranchAlias("Muon_hit_active");
    produces<vector<vector<int> > >("layernum").setBranchAlias("Muon_hit_layer");
    produces<vector<vector<int> > >("ndet").setBranchAlias("Muon_hit_ndet");
    produces<vector<vector<float> > >("x").setBranchAlias("Muon_hit_x");
    produces<vector<vector<float> > >("y").setBranchAlias("Muon_hit_y");
    produces<vector<vector<float> > >("z").setBranchAlias("Muon_hit_z");
    produces<vector<int> >("nexpectedhits").setBranchAlias("Muon_nExpectedPixelHits");
    produces<vector<int> >("ncompatible").setBranchAlias("Muon_nCompatiblePixelLayers");
    produces<vector<int> >("nexpectedhitsmultiple").setBranchAlias("Muon_nExpectedPixelHitsMultiple");
    produces<vector<int> >("ncompatibletotal").setBranchAlias("Muon_nCompatibleTrackerLayers");
    produces<vector<int> >("nexpectedhitsmultipletotal").setBranchAlias("Muon_nExpectedTrackerHitsMultiple");
    produces<vector<int> >("nexpectedhitstotal").setBranchAlias("Muon_nExpectedTrackerHits");
    produces<vector<float> >("pxatdv").setBranchAlias("Muon_pxatdv");
    produces<vector<float> >("pyatdv").setBranchAlias("Muon_pyatdv");
    produces<vector<float> >("pzatdv").setBranchAlias("Muon_pzatdv");
    produces<vector<bool> >("onmodule").setBranchAlias("Muon_patdvOnModule");
}

HitMaker::~HitMaker(){
}

void HitMaker::beginJob(){}

void HitMaker::endJob(){}

void HitMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup){}


bool HitMaker::checkDet(const GeomDet *det, float pxatdv, float pyatdv, float pzatdv){

    bool debug = 0;
    bool ret=0;
    float xy[4][2];
    float dz;
    const Bounds *b = &((det->surface()).bounds());

    if (const TrapezoidalPlaneBounds *b2 = dynamic_cast<const TrapezoidalPlaneBounds *>(b)) {
        // See sec. "TrapezoidalPlaneBounds parameters" in doc/reco-geom-notes.txt
        std::array<const float, 4> const &par = b2->parameters();
        xy[0][0] = -par[0];
        xy[0][1] = -par[3];
        xy[1][0] = -par[1];
        xy[1][1] = par[3];
        xy[2][0] = par[1];
        xy[2][1] = par[3];
        xy[3][0] = par[0];
        xy[3][1] = -par[3];
        dz = par[2];
    }
    else if (const RectangularPlaneBounds *b2 = dynamic_cast<const RectangularPlaneBounds *>(b)) {
        // Rectangular
        float dx = b2->width() * 0.5;   // half width
        float dy = b2->length() * 0.5;  // half length
        xy[0][0] = -dx;
        xy[0][1] = -dy;
        xy[1][0] = -dx;
        xy[1][1] = dy;
        xy[2][0] = dx;
        xy[2][1] = dy;
        xy[3][0] = dx;
        xy[3][1] = -dy;
        dz = b2->thickness() * 0.5;  // half thickness
    }
    //check compatibility
    for (int i = 0; i < 4; ++i) {
        Local3DPoint lp1(xy[i][0], xy[i][1], -dz);
        Local3DPoint lp2(xy[i][0], xy[i][1], dz);
        GlobalPoint gp1 = det->surface().toGlobal(lp1);
        GlobalPoint gp2 = det->surface().toGlobal(lp2);
        if (debug) {
           std::cout << gp1.z() << " z1 " << pzatdv  << std::endl;
           std::cout << gp2.z() << " z2 " << pzatdv  << std::endl;
           std::cout << gp1.y() << " y1 " << pyatdv  << std::endl;
           std::cout << gp2.y() << " y2 " << pyatdv  << std::endl;
           std::cout << gp1.x() << " x1 " << pxatdv  << std::endl;
           std::cout << gp2.x() << " x2 " << pxatdv  << std::endl;
        }
       //TO DO : check if position is compatible  
    }
   

    return ret;
}




bool HitMaker::positionOnModule(const TrackerGeometry trackerGeom, float pxatdv, float pyatdv, float pzatdv){
    bool ret = false;
    for (auto &det : trackerGeom.detsPXB()){
        ret = checkDet(det, pxatdv, pyatdv, pzatdv);
        if (ret) return ret;
    }
    for (auto &det : trackerGeom.detsPXF()){
        ret = checkDet(det, pxatdv, pyatdv, pzatdv);
        if (ret) return ret;
    }
    return ret;
}

void HitMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){

    bool debug = false;

    propagatorHandle_ = iSetup.getHandle(propagatorToken_);
    magfield_ = iSetup.getHandle(magFieldToken_);
    theGeo_ = iSetup.getHandle(trackingGeometryToken_);
    measurementTracker_ = iSetup.getHandle(measurementTrackerToken_);

    auto const& searchGeom = *(*measurementTracker_).geometricSearchTracker();
    auto const& prop = *propagatorHandle_;

    edm::Handle<MeasurementTrackerEvent> measurementTrackerEvent;
    iEvent.getByToken(measurementTrackerEventToken_, measurementTrackerEvent);

    edm::Handle<Run3ScoutingMuonCollection> muonHandle;
    iEvent.getByToken(muonToken_, muonHandle);

    edm::Handle<Run3ScoutingVertexCollection> dvHandle;
    iEvent.getByToken(dvToken_, dvHandle);

    const TrackerGeometry trackerGeom_ = iSetup.get<TrackerRecoGeometryRecord>().get(geomToken_);

    if (debug) {
        std::cout << std::endl;
        std::cout << "------- Run " << iEvent.id().run() << " Lumi " << iEvent.luminosityBlock() << " Event " << iEvent.id().event() << " -------" << std::endl;
    }

    // nlohmann::json j;

    unique_ptr<vector<vector<bool> > > v_isbarrel(new vector<vector<bool> >);
    unique_ptr<vector<vector<bool> > > v_ispixel(new vector<vector<bool> >);
    unique_ptr<vector<vector<bool> > > v_isactive(new vector<vector<bool> >);
    unique_ptr<vector<vector<int> > > v_layernum(new vector<vector<int> >);
    unique_ptr<vector<vector<int> > > v_ndet(new vector<vector<int> >);
    unique_ptr<vector<vector<float> > > v_hitx(new vector<vector<float> >);
    unique_ptr<vector<vector<float> > > v_hity(new vector<vector<float> >);
    unique_ptr<vector<vector<float> > > v_hitz(new vector<vector<float> >);
    unique_ptr<vector<int> > v_nexpectedhits(new vector<int>);
    unique_ptr<vector<int> > v_ncompatible(new vector<int>);
    unique_ptr<vector<int> > v_nexpectedhitsmultiple(new vector<int>);
    unique_ptr<vector<int> > v_nexpectedhitstotal(new vector<int>);
    unique_ptr<vector<int> > v_ncompatibletotal(new vector<int>);
    unique_ptr<vector<int> > v_nexpectedhitsmultipletotal(new vector<int>);
    unique_ptr<vector<float> > v_pxatdv(new vector<float>);
    unique_ptr<vector<float> > v_pyatdv(new vector<float>);
    unique_ptr<vector<float> > v_pzatdv(new vector<float>);
    unique_ptr<vector<bool> > v_onmodule(new vector<bool>);

    for (auto const& muon : *muonHandle) {
        vector<int> vertex_indices = muon.vtxIndx();
        int first_good_index = 0;
        for (auto idx : vertex_indices) {
            if (idx >= 0) {
                first_good_index = idx;
                break;
            }
        }
        int nDV = (*dvHandle).size();
        float dv_x = 0;
        float dv_y = 0;
        float dv_z = 0;
        if (first_good_index < nDV) {
            Run3ScoutingVertex dv = (*dvHandle).at(first_good_index);
            dv_x = dv.x();
            dv_y = dv.y();
            dv_z = dv.z();
        }
        TLorentzVector lv;
        lv.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), 0.10566);

        float track_px = lv.Px();
        float track_py = lv.Py();
        float track_pz = lv.Pz();
        float track_phi = muon.trk_phi();
        float track_dz = muon.trk_dz();
        float track_dxy = muon.trk_dxy();
        float track_dsz = muon.trk_dsz();
        float track_lambda = muon.trk_lambda();
        float track_qoverpError = muon.trk_qoverpError();
        float track_lambdaError = muon.trk_lambdaError();
        float track_phiError = muon.trk_phiError();
        float track_dxyError = muon.trk_dxyError();
        float track_dszError = muon.trk_dszError();
        int track_charge = muon.charge();
        int nvalidpixelhits = muon.nValidPixelHits();

        /* Compute track reference point position. This is where the momentum is reported.
         * Invert this to solve:
         * / dxy \   / -sin(phi)              cos(phi)              0           \ / vx \
         * | dsz | = | -cos(phi)*sin(lambda)  -sin(phi)*sin(lambda) cos(lambda) | | vy |
         * \  dz /   \ 0                      0                     1           / \ vz /
         * Based on DataFormats/TrackReco/interface/TrackBase.h
         */
        float sinphi = sin(track_phi);
        float cosphi = cos(track_phi);
        float sinlmb = sin(track_lambda);
        float coslmb = cos(track_lambda);
        float tanlmb = sinlmb/coslmb;
        float track_vz = track_dz;
        float track_vx = -sinphi*track_dxy - (cosphi/sinlmb)*track_dsz + (cosphi/tanlmb)*track_vz;
        float track_vy =  cosphi*track_dxy - (sinphi/sinlmb)*track_dsz + (sinphi/tanlmb)*track_vz;

        // Is track reference point inside a cylinder with the DV? This should always be true from what I've seen.
        bool track_ref_inside_dv = (track_vx*track_vx+track_vy*track_vy) < (dv_x*dv_x+dv_y*dv_y);

        if (debug) {
            std::cout << "=== Muon "
                << " nDV=" << vertex_indices.size()
                << " (pt,eta,phi)=(" << muon.pt() << "," << muon.eta() << "," << muon.phi() << ")"
                << " (px,py,pz)=(" << track_px << "," << track_py << "," << track_pz << ")"
                << " (dvx,dvy,dvz)=(" << dv_x << "," << dv_y << "," << dv_z << ")"
                << " (tvx,tvy,tvz)=(" << track_vx << "," << track_vy << "," << track_vz << ")"
                << " trkindv=" << track_ref_inside_dv
                << " q=" << track_charge
                << " nvalid=" << nvalidpixelhits
                << std::endl;

            // j["event"] = iEvent.id().event();
            // j["nMuon"] = (*muonHandle).size();
            // j["nDV"] = (*dvHandle).size();
            // j["Muon_pt"] = muon.pt();
            // j["Muon_eta"] = muon.eta();
            // j["Muon_phi"] = muon.phi();
            // j["Muon_px"] = track_px;
            // j["Muon_py"] = track_py;
            // j["Muon_pz"] = track_pz;
            // j["Muon_nDV"] = vertex_indices.size();
            // j["Muon_dvx"] = dv_x;
            // j["Muon_dvy"] = dv_y;
            // j["Muon_dvz"] = dv_z;
            // j["Muon_tvx"] = track_vx;
            // j["Muon_tvy"] = track_vy;
            // j["Muon_tvz"] = track_vz;
            // j["Muon_trkindv"] = track_ref_inside_dv;
            // j["Muon_charge"] = track_charge;
            // j["Muon_nvalid"] = nvalidpixelhits;
        }

        reco::TrackBase::CovarianceMatrix track_cov;
        track_cov(0,0) = pow(track_qoverpError,2);
        track_cov(1,1) = pow(track_lambdaError,2);
        track_cov(2,2) = pow(track_phiError,2);
        track_cov(3,3) = pow(track_dxyError,2);
        track_cov(4,4) = pow(track_dszError,2);

        CurvilinearTrajectoryError err(track_cov);
        // Default parameters according to https://github.com/cms-sw/cmssw/blob/master/TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorParams.h
        Chi2MeasurementEstimator estimator(30., 3., 0.5, 2.0, 0.5, 1.e12);



        GlobalVector startingMomentum(track_px, track_py, track_pz);
        GlobalPoint startingPosition;
        if (track_ref_inside_dv) {
            startingPosition = GlobalPoint(track_vx, track_vy, track_vz);
        } else {
            // If the ref point is outside the DV cylinder (due to rounding issues or DV not corresponding to the right muon), just use the DV cylinder
            startingPosition = GlobalPoint(dv_x, dv_y, dv_z);
        }

        PlaneBuilder pb;
        auto startingPlane = pb.plane(startingPosition, rotation(startingMomentum));

        TrajectoryStateOnSurface startingStateP(
                GlobalTrajectoryParameters(startingPosition, startingMomentum, track_charge, magfield_.product()),
                err, *startingPlane
                );

        // float cylx = dv_x;
        // float cyly = dv_y;
        // float cylz = dv_z;

        if (track_ref_inside_dv) {
            float dv_rho = sqrt(dv_x*dv_x+dv_y*dv_y);
            if (debug) {
                std::cout << "   Before propagating trk ref to DV cyl POS/MOM: " << startingStateP.globalPosition() << "/" << startingStateP.globalMomentum() << std::endl;
            }
            // https://github.com/cms-sw/cmssw/blob/949a7b9d2c1bfde1458e01da1c14da0cd53a0ccf/HLTriggerOffline/Muon/src/PropagateToMuon.cc#L159
            // https://github.com/cms-sw/cmssw/blob/c9b012f3388a39f64eb05980e3732d0484539f14/DataFormats/GeometrySurface/interface/Cylinder.h
            double valid_distance = prop.propagateWithPath(startingStateP, Cylinder(dv_rho)).second;
            //very seldom (1/500000 probably) the propagation seems to give seg faults - since this is very rare I just check the distance
            if(valid_distance<=0) continue;
            startingStateP = prop.propagate(startingStateP, Cylinder(dv_rho));
            // cylx = startingStateP.globalPosition().x();
            // cyly = startingStateP.globalPosition().y();
            // cylz = startingStateP.globalPosition().z();
            if (debug) {
                std::cout << "   After propagating trk ref to DV cyl POS/MOM: " << startingStateP.globalPosition() << "/" << startingStateP.globalMomentum() << std::endl;
                // j["Muon_dvcylx"] = cylx;
                // j["Muon_dvcyly"] = cyly;
                // j["Muon_dvcylz"] = cylz;
            }
        }

        float pxatdv = startingStateP.globalMomentum().x();
        float pyatdv = startingStateP.globalMomentum().y();
        float pzatdv = startingStateP.globalMomentum().z();
        bool onmodule = false;//positionOnModule(trackerGeom_, pxatdv, pyatdv, pzatdv);
       
        // or could get searchGeom.allLayers() and require layer->subDetector() enum is PixelBarrel/PixelEndcap 
        //vector<DetLayer const*> layers_pixel;
        //for (auto layer : searchGeom.pixelBarrelLayers()) layers_pixel.push_back(layer);
        //for (auto layer : searchGeom.negPixelForwardLayers()) layers_pixel.push_back(layer);
        //for (auto layer : searchGeom.posPixelForwardLayers()) layers_pixel.push_back(layer);
        vector<bool> isbarrel;
        vector<bool> ispixel;
        vector<bool> isactive;
        vector<int> layernum;
        vector<int> ndet;
        vector<float> hitx;
        vector<float> hity;
        vector<float> hitz;
        int nexpectedhits = 0;
        int nexpectedhitsmultiple = 0;
        int nexpectedhitsmultipleraw = 0;
        int nexpectedhitstotal = 0;
        int nexpectedhitsmultipletotal = 0;
        auto tsos = startingStateP;
        int ncompatible = 0;
        int ncompatibletotal = 0;
        for (auto const& layer : searchGeom.allLayers()) {
            if(debug) std::cout << (int)layer->subDetector() << " " << layer->subDetector() << " layers subdet" << std::endl;
            bool pixel=((int)layer->subDetector())<2; // pixel barrel = 0, endcap = 1
            bool compatible = layer->compatible(tsos, prop, estimator).first;
            if (debug) std::cout << "compatible=" << compatible << std::endl;
            if (pixel) ncompatible += compatible;
            ncompatibletotal += compatible;
            // auto tsos = startingStateP;
            // /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_5/src/TrackingTools/DetLayers/src/BarrelDetLayer.cc
            auto const& detWithState = layer->compatibleDets(tsos, prop, estimator);
            if (debug) std::cout << "detWithState.size()=" << detWithState.size() << std::endl;
            if (!detWithState.size()) continue;
            tsos = detWithState.front().second;
            DetId did = detWithState.front().first->geographicalId();
            MeasurementDetWithData measDet = measurementTracker_->idToDet(did, *measurementTrackerEvent);
            bool active = measDet.isActive() && measDet.isValid(); // From what I see, isValid is always true, but just be safe.
            bool barrel = layer->isBarrel();
            int seq = layer->seqNum();
            int sdet = detWithState.size();
            auto pos = tsos.globalPosition();
            if (debug) {
                std::cout << "HIT subdet=" << layer->subDetector()
                    << " layer=" << seq 
                    << " detSize=" << sdet
                    << " pos=" << pos
                    << " active=" << active 
                    << std::endl;
            }
            for (auto ds : detWithState) {
                auto did2 = ds.first->geographicalId();
                auto md = measurementTracker_->idToDet(did2, *measurementTrackerEvent);
                if (debug) std::cout << "   subhit active=" << md.isActive() << " valid=" << md.isValid() << std::endl;
                if (pixel) nexpectedhitsmultiple += md.isActive()*md.isValid();
                nexpectedhitsmultipletotal += md.isActive()*md.isValid();
                nexpectedhitsmultipleraw += 1;
            }
            isbarrel.push_back(barrel);
            ispixel.push_back(pixel);
            isactive.push_back(active);
            layernum.push_back(seq);
            ndet.push_back(sdet);
            hitx.push_back(pos.x());
            hity.push_back(pos.y());
            hitz.push_back(pos.z());
            if (pixel) nexpectedhits += active;
            nexpectedhitstotal += active;
        }
        v_isbarrel->push_back(isbarrel);
        v_ispixel->push_back(ispixel);
        v_isactive->push_back(isactive);
        v_layernum->push_back(layernum);
        v_ndet->push_back(ndet);
        v_hitx->push_back(hitx);
        v_hity->push_back(hity);
        v_hitz->push_back(hitz);
        v_nexpectedhits->push_back(nexpectedhits);
        v_ncompatible->push_back(ncompatible);
        v_nexpectedhitsmultiple->push_back(nexpectedhitsmultiple);
        v_nexpectedhitstotal->push_back(nexpectedhitstotal);
        v_ncompatibletotal->push_back(ncompatibletotal);
        v_nexpectedhitsmultipletotal->push_back(nexpectedhitsmultipletotal);
        v_pxatdv->push_back(pxatdv);
        v_pyatdv->push_back(pyatdv);
        v_pzatdv->push_back(pzatdv);
        v_onmodule->push_back(onmodule);

        if (debug) {
            std::cout <<  " valid: " << nvalidpixelhits <<  " exp: " << nexpectedhits <<  " expmultiple: " << nexpectedhitsmultiple 
                      <<  " valid-exp: " << nvalidpixelhits-nexpectedhits <<  " valid-expmultiple: " << nvalidpixelhits-nexpectedhitsmultiple <<  std::endl;
            // j["Muon_expected"] = nexpectedhits;
            // j["Muon_expectedmultiple"] = nexpectedhitsmultiple;
            // j["Muon_expectedmultipleraw"] = nexpectedhitsmultipleraw;

        }

        if (debug) {
            // std::cout << "JSON: " << j.dump(-1) << std::endl;
        }


    }

    iEvent.put(std::move(v_isbarrel), "isbarrel");
    iEvent.put(std::move(v_ispixel), "ispixel");
    iEvent.put(std::move(v_isactive), "isactive");
    iEvent.put(std::move(v_layernum), "layernum");
    iEvent.put(std::move(v_ndet), "ndet");
    iEvent.put(std::move(v_hitx), "x");
    iEvent.put(std::move(v_hity), "y");
    iEvent.put(std::move(v_hitz), "z");
    iEvent.put(std::move(v_nexpectedhits), "nexpectedhits");
    iEvent.put(std::move(v_ncompatible), "ncompatible");
    iEvent.put(std::move(v_nexpectedhitsmultiple), "nexpectedhitsmultiple");
    iEvent.put(std::move(v_nexpectedhitstotal), "nexpectedhitstotal");
    iEvent.put(std::move(v_ncompatibletotal), "ncompatibletotal");
    iEvent.put(std::move(v_nexpectedhitsmultipletotal), "nexpectedhitsmultipletotal");
    iEvent.put(std::move(v_pxatdv), "pxatdv");
    iEvent.put(std::move(v_pyatdv), "pyatdv");
    iEvent.put(std::move(v_pzatdv), "pzatdv");
    iEvent.put(std::move(v_onmodule), "onmodule");
}

DEFINE_FWK_MODULE(HitMaker);
