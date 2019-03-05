#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include "TTree.h"
#include <TH1D.h>

using namespace std;

//TMinuit requires the fit objects to be global
vector<TH1D*> h; 
TH1D *h_gaus;
unsigned step = 1;

//*********************************************************************************************
//creates a function to then feed to fcn to fit (TMinuit - ROOT)
//the number of histograms in the linear combination depends on the oaa interval and hence needs to be changed manually here if the total number of OAA cuts changes
TH1D* myfunc(vector<TH1D*> h, Double_t *par) {
  TH1D *sbnd = new TH1D("sbnd", "", 80, 0.0, 4.0); //create a histogram of linear combination of sbnd off-axis flux histograms with par as coefficients
  *sbnd = ( par[0]*(*h[0]) + par[1]*(*h[1]) + par[2]*(*h[2]) + par[3]*(*h[3]) + par[4]*(*h[4]) + par[5]*(*h[5]) + par[6]*(*h[6]) + par[7]*(*h[7]) + par[8]*(*h[8]) + par[9]*(*h[9]) + par[10]*(*h[10]) + par[11]*(*h[11]) + par[12]*(*h[12]) + par[13]*(*h[13]) + par[14]*(*h[14]) + par[15]*(*h[15]) + par[16]*(*h[16]) + par[17]*(*h[17]) );
  return sbnd;
}

//************************************************************************************************
//Feed in fake-gaussian and the OAA bin histograms to fit 
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  //defining chi-square to minimize (standard chi-square definition summed over each bin) 
  double chisq=0;
  for(Int_t i=1; i<=(h_gaus->GetNbinsX()); i++) 
  {
    if(h_gaus->GetBinContent(i) != 0) 
    {
      chisq = chisq + (pow( (myfunc(h, par)->GetBinContent(i) - h_gaus->GetBinContent(i)), 2) / h_gaus->GetBinContent(i) );
    }
  }
  cout << "  step = " << step << "  ||  chisquare = " << chisq << endl; ;
  f = chisq;
  step++;
}

//*************************************************************************************************
//declarations
TH1D* monoHist(double mean, double sigma); //function for making a "fake" gaussian histogram to fit
void intervalHists(vector<TH1D*> hist, vector<string> offaxis); //function for displaying all the off-axis flux histograms in the same canvas
void recoFit(double mean, double sigma, TH1D* reco_fit, TH1D* qe_fit, TH1D* nonqe_fit); //function to plot and compare the fit to the fake gaussian
void fitCheck(TH1D* h_gaus, TH1D* h_fit, double mean, double sigma); //function plot reco, qe, and nonqe fit histogram (result of the fit applied to reco)

//************************************************************************************************
//************************************************************************************************
void ana_prism() { 

  //input the parameters of the desired gaussian
  double mean = 0.9; 
  double sigma = 0.2; 

  vector<TCut> cut; //defining cuts to for each off-axis interval
  vector<string> offaxis; //string of cut intervals
  vector<TH1D*> reco; //dividing reco histograms according to oaa
  vector<TH1D*> h_nonqe; //------------||---------------
  vector<TH1D*> h_qe; //-------------||-------------

  TCanvas *c = new TCanvas();
  TFile *f = new TFile("output_sbndPRISM.root"); //created using sbnana
  TTree *sbnana = (TTree*)f->Get("sbnana");

  vector<double> nu_ene, oaa, nu_x, nu_y, nu_z, nu_reco_ene;
  vector<int> nu_type, nu_mode, ccnc;

//GET BRANCHES AND REDEFINE LOCALLY
  sbnana->SetBranchAddress("nu_type", &nu_type); //true - neutrino species
  sbnana->SetBranchAddress("nu_mode", &nu_mode); //true - interaction mode
  sbnana->SetBranchAddress("nu_x", &nu_x); //int vertex x
  sbnana->SetBranchAddress("nu_y", &nu_y); //int vetex y
  sbnana->SetBranchAddress("nu_z", &nu_z); //int vertex z
  sbnana->SetBranchAddress("nu_ene", &nu_ene); //true neutrino energy
  sbnana->SetBranchAddress("nu_reco_ene", &nu_reco_ene); //reco neutrino energy
  sbnana->SetBranchAddress("oaa", &oaa); //off-axis-angle measured from the center of neutrino beam at sbnd
  sbnana->SetBranchAddress("ccnc", &ccnc); //true CC or NC information

//DEFINE CUTS
  TCut nu_mu("nu_type == 14"); //ONLY MUON-NEUTRINOS
  TCut cc("ccnc == 0"); //ONLY CC INTERACTIONS
  TCut int_qe("nu_mode == 0"); //SPECIFICALLY QE INTERACTIONS
  TCut int_nonqe("(nu_mode != 0)&&(nu_mode != -1)"); //ALL INTERACTIONS THAT ARE NOT-QE AND NOT-UNKNOWN
  TCut nomec("nu_mode != 10"); //IN CASE ONE WANTS TO PLAY WITH THE FIT (REMOVING MEC SAMPLE)
 
  //sbnd ACTIVE VOLUME CUT -- MODIFY AS PER FIDUCAIL VOLUME DEFINITION
  TCut x_active("(nu_x >= -200)&&(nu_x <= 200)"); 
  TCut y_active("(nu_y >= -200)&&(nu_y <= 200)"); 
  TCut z_active("(nu_z >= 0)&&(nu_z <= 500)"); 

  double interval = 0.1; //THE INTERVAL OF OAA GROUPS (eg. 1.5-1.6, 1.6, 1.7)

//DIVING EVENTS ACC. TO OAA INTERVALS (USUAL INTERVAL: 0.1 OR 0.2 DEGREES)
  cut.push_back("(oaa >= 0.0)&&(oaa < 0.2)");
  offaxis.push_back("0.0 - 0.2");

  for(Int_t i=0; i<16; i++) 
  {
    cut.push_back(Form("(oaa >= %.2f)&&(oaa < %.2f)", (i+2)*0.1, (i+3)*0.1));
    offaxis.push_back(Form("%.1f - %.1f", (i+2)*0.1, (i+3)*0.1));
  }

  cut.push_back("(oaa >= 1.80)&&(oaa <= 2.00)");
  offaxis.push_back("1.8 - 2.0");

cout << "----------------------------------------------------------------" << endl;

//USE THE CUTS ABOVE TO DIVIDE EVENTS INTO HISTOGRAMS OF EACH CUT-INTERVAL (THESE WILL BE USED TO DO THE LINEAR COMBINATION FIT)
  for(Int_t i = 0; i<cut.size(); i++) 
  {
    h.push_back(new TH1D(Form("h%d", i), "", 80, 0.0, 4.0));
    sbnana->Draw(Form("nu_ene>>h%d", i), nu_mu && cc && cut.at(i) && x_active && y_active && z_active );

    reco.push_back(new TH1D(Form("reco%d", i), "", 80, 0.0, 4.0));
    sbnana->Draw(Form("nu_reco_ene>>reco%d", i), nu_mu && cut.at(i) && x_active && y_active && z_active );

    h_nonqe.push_back(new TH1D(Form("h_nonqe%d", i), "", 80, 0.0, 4.0));
    sbnana->Draw(Form("nu_reco_ene>>h_nonqe%d", i), nu_mu && cut.at(i) && x_active && y_active && z_active && int_nonqe );

    h_qe.push_back(new TH1D(Form("h_qe%d", i), "", 80, 0.0, 4.0));
    sbnana->Draw(Form("nu_reco_ene>>h_qe%d", i), nu_mu && cut.at(i) && x_active && y_active && z_active && int_qe );

  cout << "True, Reco, QE-reco and NON-QE reco histograms created for OAA bin = " << offaxis.at(i) << endl;
  }

cout << "----------------------------------------------------------------" << endl;  
  cout << "Total size of histogram array = " << h.size() << endl;
  cout << "Histograms created for each off-axis interval" << endl;

  int para = h.size(); //NO. OF HISTOGRAMS EQUALS THE NUMBER OF PARAMETERS FOR THE FIR (IMP VARIABLE SINCE IT WILL BE USED IN MOST LOOPS HENCEFORTH)

//FUNCTION CALL TO CREATE "FAKE" GAUSSIAN HISTOGRAM WITH THE GIVEN VALUES OF MEAN AND SIGMA
  h_gaus = monoHist(mean, sigma); 

//FUNCTION TO DISPLAY ALL HISTOGRAMS IN THE SAME CANVAS
//  intervalHists(h, offaxis);

//MAIN TMINUIT SCRIPT (TO CARRY OUT THE FIT)
  Int_t ierflg = 0;
  Double_t arglist[10];

  arglist[0] = 1;

  TMinuit *gMinuit = new TMinuit(para); //DEFINITION
  gMinuit->SetFCN(fcn); //SET THE FUNCTION TO MINIMIZE (USER DEFINED AS IN FCN ABOVE)
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

//SET THE STARTING VALUES AND INTERVAL STEPS FOR EACH PARAMTER OF THE FIT
  static double start = 0.0;
  static double step = 0.05;

  for(Int_t i=0; i<para; i++) 
  {  
    gMinuit->mnparm(i, Form("a%d",i), start, step, 0, 0, ierflg);
  } 

  arglist[0] = 500; //MAX NUMBER OF STEPS  
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg); //ALGORITHM USED FOR THE FIR (REFER TO TMINUIT FOR OTHER OPTIONS)

//PRINTING OUT THE FIR RESULTS
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);

cout << "----------------------------------------------------------------" << endl;
  if(gMinuit->GetStatus() == 0) cout << "******************** FIT CONVERGED ****************" << endl;
  else cout << "******************* FIT DID NOT; I REPEAT; DID NOT CONVERGE ****************" << endl;

  gMinuit->mnimpr(); //TRY TO IMPROVE THE MINIMUM (WILL USUALLY GET BACK TO THE INITIAL MINIMUM AND HENCE COMMENT OUT IF IT TAKES TOO LONG)
cout << "----------------------------------------------------------------" << endl;

//TMINUIT IS COMPLETE, GETTING FIT RESULTS AND USING THEM TO MAKE RESULT HISTOGRAMS 
  vector<double> coeff;
  vector<double> error;

  for(Int_t i=0; i<para; i++) 
  {
    double w0, e0;
    gMinuit->GetParameter(i, w0, e0);
    coeff.push_back(w0);
    error.push_back(e0);
    cout << "Paramter " << i << " : Value = " << w0 << " ; Error = " << e0 << endl;
  }

cout << "\n----------------------------------------------------------------" << endl;
  cout << "Number of Parameters = " << para << "\n" << endl;
  
//DEFINING HISTOGRAMS TO DRAW THE FIT RESULTS
  TH1D *h_fit = new TH1D("h_fit", "", 80, 0.0, 4.0); //TRUE ENERGY
  TH1D *reco_fit = new TH1D("reco_fit", "", 80, 0.0, 4.0); //RECO ENERGY
  TH1D *nonqe_fit = new TH1D("nonqe_fit", "", 80, 0.0, 4.0); //NON-QE (MEC) RECO INTERACTIONS
  TH1D *qe_fit = new TH1D("qe_fit", "", 80, 0.0, 4.0); //QE RECO INTERACTIONS

//SCALING AND SUMMING EACH OF THE HSITOGRAMS USING THE FIT PARAMETERS AS COEFFICIENTS
  for(Int_t i=0; i<para; i++) 
  {
    *h_fit = *h_fit + coeff[i]*(*h[i]); //TRUE ENERGY
    *reco_fit = *reco_fit + coeff[i]*(*reco[i]); //RECO ENERGY
    *nonqe_fit = *nonqe_fit + coeff[i]*(*h_nonqe[i]); //NON-QE (MEC) RECO  
    *qe_fit = *qe_fit + coeff[i]*(*h_qe[i]); //QE RECO
  }
  cout << "Linear Combinations histograms made from the histogram arrays and the weights from the fit" << endl;

//FUNCTION CALL TO PLOT THE RESULT OF LINEAR COMB. OF FIT FOR TRUE ENERGY HISTOGRAMS (OG CHECK FOR HOW GOOD THE FIT IS)
  fitCheck(h_gaus, h_fit, mean, sigma);

//FUNCTION CALL TO PLOT THE RESULT OF LINEAR COMB. OF RECO, NON-QE AND QE HISTOGRAMS (VISIBLE SEPERATION BETWEEN NON-QE AND QE PEAKS IS THE GOAL) 
  recoFit(mean, sigma, reco_fit, qe_fit, nonqe_fit);

}

//************************************************************************************************************
//************************************************************************************************************

//PLOTS THE RESULT FROM TMINUIT FIT WITH THE ORIGIANL FAKE HISTOGRAM TO VISIBLY SHOW THE GOODNESS OF FIT
void fitCheck(TH1D* h_gaus, TH1D* h_fit, double mean, double sigma) {
  
  TCanvas *c_check = new TCanvas();
  h_gaus->SetTitle(";Neutrino Energy (GeV); Counts / 0.5 GeV");
  h_gaus->GetXaxis()->CenterTitle();
  h_gaus->GetYaxis()->CenterTitle();
  h_gaus->SetLineColor(kBlue);
  h_gaus->SetStats(kFALSE);
  h_gaus->Draw("hist");
  h_fit->SetLineColor(kBlack);
  h_fit->Draw("ep same");

  TLegend *leg = new TLegend(0.60, 0.70, 0.90, 0.85);
  leg->AddEntry(h_gaus, "(Fake) Mono Energy Gaussian", "l");
  leg->AddEntry(h_fit, "Fit (SBND offaxis Fluxes)", "lep");
  leg->SetBorderSize(0);
  leg->Draw();

  TLatex *text = new TLatex(0.6, 0.6, Form("#splitline{Fake Histogram:}{Mean = %.2f; Sigma = %.2f}", mean, sigma));
  text->SetNDC(kTRUE);
  text->SetTextFont(43);
  text->SetTextSize(14);
  text->Draw();

  TCanvas *c_fit = new TCanvas();
  h_fit->SetStats(kFALSE);
  h_fit->GetXaxis()->SetRangeUser(0., 3.);
  h_fit->Draw("ep");
  h_fit->Fit("gaus");
  h_fit->SetTitle(";Neutrino Energy (GeV); Counts / 50 MeV");

  TF1 *fit_func = (TF1*)h_fit->GetFunction("gaus");

  TLegend *leg1 = new TLegend(0.63, 0.75, 0.90, 0.85);
  leg1->AddEntry(h_fit, "Fit (SBND offaxis Fluxes)", "lep");
  leg1->AddEntry(fit_func, "Gaussian Fit to Histogram fit", "l");
  leg1->SetBorderSize(0);
  leg1->Draw();

  float x0 = fit_func->GetParameter(1);
  float s0 = fit_func->GetParameter(2);
  float chi2 = fit_func->GetChisquare();
  int ndf = fit_func->GetNDF();

  TLatex *text1 = new TLatex(0.63, 0.68, Form("#splitline{Output from Gaussian Fit:}{#splitline{Mean = %.3f; Sigma = %.3f}{#splitline{Chi2 = %.3f; NDF = %d}{Chi2/NDF = %.3f}}}", x0, s0, chi2, ndf, chi2/ndf));
  text1->SetNDC(kTRUE);
  text1->SetTextFont(43);
  text1->SetTextSize(14);
  text1->Draw();
cout << "From fitting gaussian to resulting histogram:" << endl;
cout << "Mean = " << x0 << " || Std Dev = " << s0 << endl;
cout << " ChiSquare = " << chi2 << " || NDF = " << ndf << " || Chi2/NDF = " << chi2/ndf << endl;
}

//***********************************************************************************************
//PLOTS THE RESULT OF RECO, NON-QE AND QE FITED HISTOGRAMS ON THE SAME CANVAS
void recoFit(double mean, double sigma, TH1D* reco_fit, TH1D* qe_fit, TH1D* nonqe_fit) {

  int x0 = mean*10;
  int s0 = sigma*10;

  TFile *f1 = new TFile(Form("result/reco/nonqe/nonqe_mean_%02d_sigma_%02d.root", x0, s0), "RECREATE");
  TCanvas *c_reco = new TCanvas(Form("c_mean_%02d", x0),"");

  reco_fit->SetTitle(";Reco #nu Energy (GeV); Counts / 50 MeV");
  reco_fit->SetLineColor(kBlack);
  nonqe_fit->SetLineColor(kRed);
  qe_fit->SetLineColor(kBlue);
  reco_fit->GetXaxis()->SetRangeUser(0,3);
  reco_fit->SetStats(kFALSE);
  reco_fit->Draw("hist");
  nonqe_fit->Draw("hist same");
  qe_fit->Draw("hist same");

  TLegend *leg3 = new TLegend(0.55, 0.70, 0.90, 0.85);
  leg3->SetHeader("all non-QE", "C");
  leg3->AddEntry(reco_fit, "reco-energy histograms", "lep");
  leg3->AddEntry(nonqe_fit, "nonQE (reco) histograms", "lep");
  leg3->AddEntry(qe_fit, "QE (reco) histograms", "lep");
  leg3->SetBorderSize(0);
  leg3->Draw();

  TLatex *text = new TLatex(0.4, 0.91, Form("Reco mono-gaussian with input mean = %.2f GeV", mean));
  text->SetNDC(kTRUE);
  text->SetTextFont(43);
  text->SetTextSize(15);
  text->Draw();

  reco_fit->Write(Form("nonqe_reco_%02d", x0));
  qe_fit->Write(Form("nonqe_qe_%02d", x0));
  nonqe_fit->Write(Form("nonqe_nonqe_%02d", x0));

  c_reco->SaveAs(Form("result/reco/nonqe/c_nonqe_mean_%02d.png", x0));

  f1->Close();

}

//***********************************************************************************************
//CREATING A GAUSSIAN HISTOGRAM TO FIT THE TRUE ENERGY HISTOGRAMS TO
TH1D* monoHist(double mean, double sigma) {

  double x0 = mean;
  double s0 = sigma;
  TF1 *f_gaus = new TF1("f_gaus", Form("TMath::Gaus(x, %f, %f)", x0, s0), 0.0, 4.0);
  TH1D *gaus = new TH1D("gaus", "", 80, 0.0, 4.0);
  gaus->FillRandom("f_gaus", 20000);
//  gaus->Draw();

  return gaus; //return the oscillted icarus flux

}

//**********************************************************************
//PLOTING ALL RESULTING OFF-AXIS HISTOGRAMS IN THE SAME CANVAS
void intervalHists(vector<TH1D*> hist, vector<string> offaxis) {

  TCanvas *c_all = new TCanvas();
  c_all->Divide(5, 4, 0.003, 0.003);
  for(Int_t i=0; i<hist.size(); i++) 
  {
    c_all->cd(i+1);

    hist.at(i)->GetXaxis()->SetRangeUser(0.0, 3.0);
    hist.at(i)->SetTitle(";Neutrino Energy (GeV); Counts / 0.5 GeV");
    hist.at(i)->GetXaxis()->CenterTitle();
    hist.at(i)->GetXaxis()->SetTitleSize(0.06);
    hist.at(i)->GetXaxis()->SetLabelSize(0.05);
    hist.at(i)->GetXaxis()->SetTitleOffset(0.80);
    hist.at(i)->GetYaxis()->CenterTitle();
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(0.95);
    hist.at(i)->Draw("hist");

    double rms = hist.at(i)->GetRMS();
    double binmax = hist.at(i)->GetMaximumBin();
    double maxene = hist.at(i)->GetXaxis()->GetBinCenter(binmax);

    TLatex *t = new TLatex(0.45, 0.8, Form("#splitline{OAA = %s}{#splitline{Peak Ene = %g}{RMS = %.4f}}", offaxis.at(i).c_str(), maxene, rms));
    t->SetNDC(kTRUE);
    t->SetTextFont(43);
    t->SetTextSize(12);
    t->Draw();
  }
}


