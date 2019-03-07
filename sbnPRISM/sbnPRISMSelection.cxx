#include <iostream>
#include <math.h>
#include <vector>
#include <TH2D.h>
#include <TH1D.h>
#include <TH3D.h>
#include <THStack.h>
#include <TLorentzVector.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "sbnPRISMSelection.h"
#include "sbnPRISMTools.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"

namespace ana {
  namespace sbnPRISM {

sbnPRISMSelection::sbnPRISMSelection() : SelectionBase(), fEventCounter(0), nu(0), muon(0), problem(0) {}


void sbnPRISMSelection::Initialize(Json::Value* config) {

  // Make a histogram
  h_reco_ene = new TH1D("h_reco_ene", "", 200, 0.0, 10.0); 
  h_nu_ene = new TH1D("h_nu_ene", "", 200, 0.0, 10.0);
  h_true_loop = new TH1D("h_true_loop", "", 200, 0.0, 10.0);

//Load configuration parameters
  fMyParam = 0;
  fTruthTag = { "generator" };
  fFluxTag = { "generator" };
  fPartTag = { "largeant" };  

  if (config) {
    fMyParam = (*config)["sbnPRISM"].get("parameter", 0).asInt();
    fTruthTag = { (*config)["sbnPRISM"].get("MCTruthTag", "generator").asString() };
    fPartTag = { (*config)["sbnPRISM"].get("MCPartTag", "largeant").asString() };
  }

  // Add custom branches
   AddBranch("nu_int_vertex", &bNuVertex); //mctruth nu interaction vertex
   AddBranch("oaa", &bAngle); //angle between the neutrino vector and the center of the beam 
   AddBranch("nu_ene", &bEne); //mctruth nu energy
   AddBranch("nu_reco_ene", &bNuRecoEne); //reconstructed neutrino energy from mcparticles using numuCC formula
   AddBranch("nu_type", &bNuType);
   AddBranch("nu_mode", &bNuMode);
   AddBranch("nu_x", &bNuX);
   AddBranch("nu_y", &bNuY);
   AddBranch("nu_z", &bNuZ);
   AddBranch("nu_ene_loop", &bNuEneLoop);
   AddBranch("ccnc", &bNuCCNC);
}


void sbnPRISMSelection::Finalize() {
  // Output our histograms to the ROOT file
  fOutputFile->cd();
  
  h_reco_ene->Write();
  h_nu_ene->Write();
  h_true_loop->Write();

for(int i=0; i<check.size(); i++) {
std::cout << "reco energy = " << check[i] << " || numer = " << numer[i] << " || denomi = " << denomi[i] << " || mom_mag = " << mom_mag[i] << " || costheta = " << costheta[i] << std::endl;
}
std::cout << "Total Neutrinos = " << nu << " || muon (mcparticle) = " << muon << std::endl;
std::cout << "Problem = " << problem << std::endl;
}


bool sbnPRISMSelection::ProcessEvent(gallery::Event& ev) {

//clearing all branch vectors for the new event
//branches are vectors to get each mctruth instead of writing the value for only the first mctruth in the event
  bAngle.clear();
  bEne.clear();
  bNuVertex.clear(); 
  bNuRecoEne.clear();
  bNuType.clear();
  bNuMode.clear();
  bNuX.clear();
  bNuY.clear();
  bNuZ.clear();
  bNuEneLoop.clear();
  bNuCCNC.clear();

  if (fEventCounter % 10 == 0) {
    std::cout << "sbnPRISMSelection: Processing event " << fEventCounter << std::endl;
      }
  fEventCounter++;

  // Grab a data product from the event
  auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);
  auto const& mctruth_handle = ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);
 
  //Association between MCTruth and MCParticle
  const art::FindManyP <simb::MCParticle, sim::GeneratedParticleInfo> find_many_mcpart(mctruth_handle, ev, fPartTag);
   
  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) 
  {
    nu++; //number of neutrino interactions

    //grab each neutrino interaction in the event
    auto const& mctruth = mctruths.at(i);

    const simb::MCNeutrino& nu = mctruth.GetNeutrino();

    //calculating and storing the off axis angle
    TLorentzVector nuVertex = nu.Nu().Position();
    double nuX = nuVertex.X(), nuY = nuVertex.Y(), nuZ = nuVertex.Z();

    h_nu_ene->Fill(nu.Nu().E());

    //fill branches
    bNuVertex.push_back(nuVertex);
    bNuX.push_back(nuX);
    bNuY.push_back(nuY);
    bNuZ.push_back(nuZ);
    bEne.push_back(nu.Nu().E()); //neutrino energy
    bNuType.push_back(nu.Nu().PdgCode());
    bNuMode.push_back(nu.Mode());
    bNuCCNC.push_back(nu.CCNC());

    //calcuatinf off-axis angle with respect to center of neutrino beam (known)
    double nuVerMag = sqrt( pow((nuX-cX),2) + pow((nuY-cY),2) + pow((nuZ-cZ+10950),2) );
    double cosine = ((nuZ-cZ+10950)/nuVerMag);
    double oaa = acos(cosine)*180/pi;
 
    bAngle.push_back(oaa);

   //getting the reconstructed energies using numuCC formula on the mcparticles associated with this neutrino
    std::vector <art::Ptr<simb::MCParticle>> const& mcparticle = find_many_mcpart.at(i);

    for (auto const& mcpart : mcparticle) {
      if (mcpart->PdgCode() == 13) {
        muon++;

    //calculating reconstructed energy using the qe interaction formula   
        double px = mcpart->Px(), py = mcpart->Py(), pz = mcpart->Pz(), El = mcpart->E(); //muon momentum 4-vector
        double pmag = TMath::Sqrt( px*px + py*py + pz*pz ), pxy = TMath::Sqrt( px*px + py*py ), theta = TMath::ACos(pz/pmag);

        double num =  ( mp*mp - (mn - Eb)*(mn - Eb) - ml*ml + 2*(mn - Eb)*El );
        double deno = 2.0*((mn - Eb) - El + pmag*TMath::Cos(theta));
        nuRecoEne = num/deno;
        
        if((nuRecoEne > 0.0)&&(nuRecoEne < 10.0)) {
          bNuRecoEne.push_back(nuRecoEne);
          h_reco_ene->Fill(nuRecoEne);
          bNuEneLoop.push_back(nu.Nu().E());
          h_true_loop->Fill(nu.Nu().E());
        }
        else { 
        problem++;
        numer.push_back(num);
        denomi.push_back(deno);
        mom_mag.push_back(pmag);
        costheta.push_back(cos(theta));     
        check.push_back(nuRecoEne);
        }

      } //muon if loop
    } //mcparticle 
  } //mctruth

 return true;
}

  }  // namespace sbnPRISM
}  // namespace ana

// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::sbnPRISM::sbnPRISMSelection)

