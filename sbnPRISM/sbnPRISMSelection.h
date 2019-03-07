#ifndef __sbnanalysis_ana_sbnPRISM_sbnPRISMSelection__
#define __sbnanalysis_ana_sbnPRISM_sbnPRISMSelection__

#define _USE_MATH_DEFINES
/**
 * \file sbnPRISMSelection.h
 *
 * This is an implementation of the core::SelectionBase class. We define
 * the methods called for initialization, finalization, and event-by-event
 * processing. The selection uses MCTruth and MCParticle objects to get the truth and reconstructed energies
 * in addition to other variables of interest
 *
 * Author:
 */

// Includes
#include <iostream>
#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"

// Forward declarations
class TH2D;

// All analysis code is defined in namespace "ana" 
namespace ana {

  // Code specific to the sbnPRISM
  namespace sbnPRISM {

class sbnPRISMSelection : public core::SelectionBase {
public:
  /** Constructor. */
  sbnPRISMSelection();

  /**
   * Initialization.
   *
   * Here we load configuration parameters, set up histograms for output, and
   * add our own branches to the output tree.
   *
   * \param config A configuration, as a JSON object
   */
  void Initialize(Json::Value* config=NULL);

  /** Finalize and write objects to the output file. */
  void Finalize();

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \return True to keep event
   */

  bool ProcessEvent(gallery::Event& ev);

protected:
  unsigned fEventCounter;  //!< Count processed events
  double pi = 3.141592653589793238463;

  /** Configuration parameters */
  art::InputTag fTruthTag; 
  art::InputTag fFluxTag;
  art::InputTag fTrackTag;
  art::InputTag fPartTag;

  int fMyParam;  //!< A parameter from the configuration file
  /** Custom data branches */
  std::vector<double> bEne;   
  std::vector<double> bNuRecoEne;
  std::vector<double> bAngle;
  std::vector<TLorentzVector> bNuVertex;
  std::vector<double> bNuX;
  std::vector<double> bNuY;
  std::vector<double> bNuZ;
  std::vector<double> bNuEneLoop;
  std::vector<int> bNuType;
  std::vector<int> bNuMode;
  std::vector<int> bNuCCNC;

  std::vector<double> check;
  std::vector<double> numer;
  std::vector<double> denomi;
  std::vector<double> mom_mag;
  std::vector<double> costheta;
 
  unsigned nu, muon, problem;
  double nuRecoEne;
  double nuX, nuY, nuZ;
//center of the neutrino beam at SBND 
  double cX = 130.1 , cY = -0.516 , cZ = 0.0; 
//constants to calculate reconstructed neutrino energy
  double mn = 0.93956541 /*mn = neutron mass */, mp = 0.93827208 /*mp = proton mass*/, ml = 0.10565837 /*ml = lepton mass*/, Eb = 0.030 /*Eb = binding energy*/; 

  /** HISTOGRAMS **/
  TH1D *h_reco_ene, *h_nu_ene, *h_true_loop;  
};

  }  // namespace sbnPRISM
}  // namespace ana

#endif  // __sbnanalysis_ana_sbnPRISM_sbnPRISMSelection__

