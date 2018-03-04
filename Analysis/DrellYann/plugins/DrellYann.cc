// -*- C++ -*-
//
// Package:    Analysis/DrellYann
// Class:      DrellYann
// 
/**\class DrellYann DrellYann.cc Analysis/DrellYann/plugins/DrellYann.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gourab Saha
//         Created:  Mon, 19 Feb 2018 13:45:26 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
//header for Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TLorentzVector.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DrellYann : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit DrellYann(const edm::ParameterSet&);
  ~DrellYann();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------//
  
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  TFileDirectory* histoDir;
  bool verbosity_;
  bool isMC_;
  TH1I *h_Size;
  TH1I *h_Charge;
  TH1F *h_P;
  TH1F *h_PtAll;
  TH1F *h_PtBarrel;
  TH1F *h_PtEnd;
  TH1F *h_invM;
  TH1F *h_Eta;
  TH1F *h_Phi;
  TH1F *h_E;
  TH1F *h_Dxy;
  TH1F *h_Dz;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DrellYann::DrellYann(const edm::ParameterSet& iConfig):
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexSrc"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  verbosity_(iConfig.getParameter<bool>("verbosity")),
  isMC_(iConfig.getParameter<bool>("isMC"))
{
   //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  histoDir = new TFileDirectory(fs->mkdir("muon"));
  h_Size           = histoDir->make<TH1I>("h_Size" , "Muon_size" , 20 , 0 , 20);
  h_Charge         = histoDir->make<TH1I>("h_Charge" , "Muon_charge" , 5 , -2 , 2);
  h_P              = histoDir->make<TH1F>("h_P" , "Muon_momentum" , 300 , 0. , 150.);
  h_PtAll          = histoDir->make<TH1F>("h_PtAll" , "Muon_pt_All" , 300 , 0. , 150.);
  h_PtBarrel       = histoDir->make<TH1F>("h_PtBarrel" , "Muon_pt_Barrel" , 300 , 0. , 150.);
  h_PtEnd          = histoDir->make<TH1F>("h_PtEndCap" , "Muon_pt_Endcap" , 300 , 0. , 150.);
  h_Eta            = histoDir->make<TH1F>("h_Eta" , "Muon_Eta" , 200 , -5. , 5.);
  h_Phi            = histoDir->make<TH1F>("h_Phi" , "Muon_Phi" , 200 , -5. , 5.);
  h_E              = histoDir->make<TH1F>("h_E" , "Muon_energy" , 400 , 0. , 200.);
  h_Dxy            = histoDir->make<TH1F>("h_Dxy" , "Muon_Dxy" , 200 , -5. , 5.);
  h_Dz             = histoDir->make<TH1F>("h_Dz" , "Muon_Dz" , 200 , -5. , 5.);
  h_invM           = histoDir->make<TH1F>("h_invM" , "Invariant_mass" , 400 , 10. , 210.);
} 


DrellYann::~DrellYann()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DrellYann::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  using namespace edm;
  edm::Handle<reco::VertexCollection> primaryVertices; 
  iEvent.getByToken(vertexToken_, primaryVertices);

  const reco::Vertex& vit = primaryVertices->front();
  //std::cout << "Number of Primary Vertices in the event=" 
  //				<< primaryVertices->size() << std::endl;
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  //std::cout << "Number of Muons in the event=" 
  //		      << muons->size() << std::endl;
  h_Size -> Fill(muons->size());
  //if (muons->size() < 2) std::continue;
  std::vector<pat::Muon> mu_vec;
  mu_vec.clear();
  TLorentzVector L1, L2, L;
  for ( unsigned int i = 0; i < muons->size(); ++i ) { 
    const pat::Muon& mu = muons->at(i);
    h_Charge ->Fill(mu.charge());
    h_P      ->Fill(mu.p());
    h_PtAll  ->Fill(mu.pt());
    if(std::abs(mu.eta()) < 1.2) h_PtBarrel->Fill(mu.pt());
    else h_PtEnd->Fill(mu.pt());

    h_Eta    ->Fill(mu.eta());
    h_Phi    ->Fill(mu.phi());
    h_E      ->Fill(mu.energy());
    //How to get the dxy and dz w.r.t a vertex                                                                                                                            
    reco::TrackRef tk = mu.muonBestTrack();
    double dxyWrtPV = -99.;
    double dzWrtPV = -99.;
    if(tk.isNonnull()) {
      dxyWrtPV = tk->dxy(vit.position());
      dzWrtPV = tk->dz(vit.position());
    }
    //    std::cout << "Dxy=" << dxyWrtPV << "\tdzWrtPV=" << dzWrtPV << std::endl;

    h_Dxy    ->Fill(dxyWrtPV);
    h_Dz     ->Fill(dzWrtPV);

    if (mu.isGlobalMuon() == 1 && fabs(mu.pt()) > 5 && fabs(mu.eta()) < 2.4 && fabs(dxyWrtPV) < 0.5 && fabs(dzWrtPV) < 1) mu_vec.push_back(mu);
  }
  //  std::cout<<"SSSSSSSSSSSIze:    "<<mu_vec.size()<<std::endl;
  if (mu_vec.size() >= 2) {
    for (unsigned int i = 0; i < mu_vec.size(); ++i){
      const pat::Muon& imu = mu_vec.at(i);
      int i_chrg = imu.charge();
      L1.SetPtEtaPhiE(imu.pt(), imu.eta(), imu.phi(), imu.energy());
      for (unsigned int j = i+1; j < mu_vec.size(); ++j){
	const pat::Muon& jmu = mu_vec.at(j);
	L2.SetPtEtaPhiE(jmu.pt(), jmu.eta(), jmu.phi(), jmu.energy());
	int j_chrg = jmu.charge();
	if ((i_chrg + j_chrg) == 0){
	  L = L1 + L2;
	  std::cout<<"inv_mass: "<<L.M()<<std::endl;
	  h_invM->Fill(L.M());
	}
      }
    }
  }
    
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
DrellYann::beginJob(){
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DrellYann::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DrellYann::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DrellYann);
