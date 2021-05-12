//
//  KinematicFitter_code.hpp
//
//
//  Created by chyde on 10/May/21.
//

#ifndef KinematicFitter_code_hpp
#define KinematicFitter_code_hpp

#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TNtuple.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <math.h>
#include "clas12reader.h"
#include <TRandom3.h>
#include <TLegend.h>
#include <TF1.h>
#include <TF2.h>
#include "TMinuit.h"
//#include <TAxis.h>
//#include "histogram.h"
#include "region_particle.h"
//#include "TCanvas_plot.h"
//#include "header.h"
using namespace clas12;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
    p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
               rp->par()->getPz(),p4.M());
    
}
void fcn(int &npar, double *gin, double &f, double *par, int iflag);
//Double_t GaussCosh(Double_t *x, Double_t *par);
bool EC_hit_position_fiducial_cut_homogeneous(region_part_ptr part);
bool EC_hit_position_fiducial_cut(region_part_ptr part);
bool EC_sampling_fraction_cut(region_part_ptr part);
bool DC_hit_position_counts_fiducial_cut(region_part_ptr part, int region);

bool minimal_PCAL_energy_deposition(region_part_ptr part);
bool z_vertex_cut(region_part_ptr part);

//  Physical constants
const double Mp = db->GetParticle(2212)->Mass(); //0.93827;
const double Mp2=Mp*Mp;
const double Me=db->GetParticle(11)->Mass();
const double Mn=db->GetParticle(2112)->Mass();
const double Mphoton=db->GetParticle(22)->Mass();
const double Mpi=db->GetParticle(211)->Mass(); //0.13957;


const double speedoflight=29.9792458;
const Double_t Coupling_const=1./137;
const Double_t rad2deg=180./TMath::Pi();

const Int_t nparMINUIT = 12; //47;

// Define global scope variables for Kinematic Fit
TLorentzVector beam, target, electron, electrontrigger, proton, piplus, piminus;
TLorentzVector elec4vec_fit, prot4vec_fit, pplus4vec_fit, pminus4vec_fit;
////Sigma_values for kinematic fit (not sure these values are correct)need to check
const double sigma_p_electron=0.01002;
const double sigma_p_proton=0.009999;
const double sigma_p_piplus=0.01001;
const double sigma_p_piminus=0.01001;
const double sigma_theta_electron=0.01146;
const double sigma_theta_proton=0.01533;
const double sigma_theta_piplus=0.01266;
const double sigma_theta_piminus=0.009343;
const double sigma_phi_electron=0.004665;
const double sigma_phi_proton=0.002054;
const double sigma_phi_piplus=0.002213;
const double sigma_phi_piminus=0.001622;
const double err_init =0.005;
Double_t gradient[nparMINUIT], val[nparMINUIT],parval[nparMINUIT];
const double wgt_MX2 = Mp*0.2
const double wgt_EX    = Mp*0.1
const double wgt_PY = Mp*0.1



#endif /* KinematicFitter_code_hpp */
