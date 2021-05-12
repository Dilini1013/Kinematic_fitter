//
//  KinematicFitter_hist.hpp
//
//
//  Created by chyde on 10/May/21
//  Global Declaration of Histograms
//

#ifndef KinematicFitter_hist_hpp
#define KinematicFitter_hist_hpp

auto* hmiss=new TH1D("missM","missM",200,0.0,2.0);
auto* hm2g=new TH1D("m2g","m2g",200,0,1);
auto* h_m2k=new TH1D("m2g","M_(K+K-)",200,0,2);
auto* hm2gCut=new TH1D("m2gCut","m2g",200,0,1);
auto* ht=new TH1F("ht","ht;(t0-t)",200,-6.,6.);
auto* hW2 = new TH1D("hW","hW;W(GeV)" ,100, 0.0, 6.0);
auto* hQ2 = new TH1D("hQ2","hQ2;Q^{2}(GeV^{2})",100, 0, 6.);
auto* hQ2_cut = new TH1D("hQ2_cut","hQ2;Q^{2}(GeV^{2})",100, 0, 6.);
auto* h_XBj = new TH1D("hXBj","hXBj;x_{B}",100, 0., 1.);
auto* h_Q2_xBj = new TH2D("h_Q2_xBj","h_Q2_xBj;x_{B};Q2(GeV^{2}",100,0.,1.,100,0.,12.);
auto* h_Q2_xBj_cut = new TH2D("h_Q2_xBj_cut","h_Q2_xBj;x_{B};Q2(GeV^{2}",100,0.,1.,100,0.,12.);
auto* htext1 = new TH2D("htext1","h_Q2_xBj;x_{B};Q2(GeV^{2}",10,0.,1.,100,0.,12.);
auto* htext2 = new TH2D("htext2","h_Q2_xBj;x_{B};Q2(GeV^{2}",10,0.,1.,100,0.,12.);
auto* h_tdiff = new TH1D("tmin-tmax","tmin-tmax",100, 0, 10.);
auto* h_elecEnergy = new TH1D("elecEnergy","electron Energy",100, 0, 12.);
auto* h_MpipiSq_fcut = new TH1D("h_MpipiSq_fcut","Two-pion Invariant Mass^{2};  M_{#pi#pi}^{2} (GeV^{2}) ; Events ",
                                100.,0.0,1.5);

auto* h_pion_invariant_all_cuts = new TH1D("h_pion_invariant_all_cuts","Two-pion Invariant Mass;  M_{#pi#pi} (GeV) ; Events ",
                                           100.,0.2,1.6);
auto* h_pion_invariant_1 = new TH1D("h_pion_invariant","Two-pion Invariant Mass;  M_{#pi#pi} (GeV) ; Events ",
                                    100.,0.0,2.0);
auto* h_proton_piplus_IM_all_cuts_sigmaR = new TH1D("h_proton_piplus_IM_sigmaR","(Proton+piplus) Invariant Mass in mpipi < 0.6;  M_{p#pi^{+}} (GeV) ; Events ",
                                                    200.,1.0,3.0);
auto* h_proton_piminus_IM_all_cuts_sigmaR = new TH1D("h_proton_piminus_IM_sigmaR","(Proton+piminus) Invariant Mass in mpipi < 0.6;  M_{p#pi^{-}} (GeV) ; Events ",
                                                     200.,1.0,3.0);

auto* h_proton_piplus_IM_all_cuts_rhoR = new TH1D("h_proton_piplus_IM_rhoR","(Proton+piplus) Invariant Mass in 0.67 <  mpipi < 0.87;  M_{p#pi^{+}} (GeV) ; Events ",
                                                  200.,1.0,3.0);
auto* h_proton_piminus_IM_all_cuts_rhoR = new TH1D("h_proton_piminus_IM_rhoR","(Proton+piminus) Invariant Mass in 0.67 < mpipi < 0.87  ;  M_{p#pi^{-}} (GeV) ; Events ",
                                                   200.,1.0,3.0);
auto* h_M122_calc = new TH1D("h_MX2","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                             100.,-1.0,3.0);

auto* h_M122_calctrigger_2 = new TH1D("h_MX2trigger_2","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                      100.,-1.0,2.0);
auto* h_M122_calctrigger_3 = new TH1D("h_MX2trigger_3","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                      100.,-1.0,2.0);
auto* h_M122_calctrigger_4 = new TH1D("h_MX2trigger_4","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                      100.,-1.0,2.0);
auto* h_M122_calctrigger_5 = new TH1D("h_MX2trigger_5","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                      100.,-1.0,2.0);
auto* h_M122_calctrigger_6 = new TH1D("h_MX2trigger_6","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                      100.,-1.0,2.0);

auto* h_M122_calctrigger_51 = new TH1D("h_MX2trigger_5","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                       100.,-1.0,2.0);
auto* h_M122_calctrigger_61 = new TH1D("h_MX2trigger_6","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                       100.,-1.0,2.0);

auto* h_M122_calctrigger_71 = new TH1D("h_MX2trigger_6","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                       100.,-1.0,2.0);


auto* h_Miss_proton_fcut = new TH1D("h_Miss_proton_fcut","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                    600.,-1.0,2.0);
auto* h_Miss_proton_all_cuts = new TH1D("h_Miss_proton_all_cuts","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                        600.,-1.0,2.0);

auto* h_Miss_pim_fcut = new TH1D("h_Miss_pim_fcut","Missing pi- Mass^{2};  M_{X=#pi^-}^{2} (GeV^{2}) ; Events ",
                                 600.,-2.0,2.0);

auto* h_Miss_pim_all_cuts = new TH1D("h_Miss_pim_all_cuts","Missing pi- Mass^{2};  M_{X=#pi^-}^{2} (GeV^{2}) ; Events ",   600.,-2.0,2.0);

auto* h_Miss_pip_fcut= new TH1D("h_Miss_pip_fcut","Missing pi+ Mass^{2};  M_{X=#pi^+}^{2} (GeV^{2}) ; Events ",
                                600.,-2.0,2.0);


auto* h_Miss_pip_all_cuts= new TH1D("h_Miss_pip_all_cuts","Missing pi+ Mass^{2};  M_{X=#pi^+}^{2} (GeV^{2}) ; Events ",600.,-2.0,2.0);
auto* hDiff = new TH1D("hDiff","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                       100.,-1.0,3.0);
auto* h_M122_calctriggerwithBremsphoton = new TH1D("h_MX2triggerwithBremsphoton","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                                   100.,-1.0,3.0);


auto* h_M122_calchighE = new TH1D("h_MX2highE","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                                  100.,-1.0,3.0);

auto* h_M124_calc = new TH1D("h_MX2_m4","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                             100.,-1.0,3.0);
auto* h_M122_m2 = new TH1D("h_MX2_m2","Missing proton Mass^{2};  M_{X=p}^{2} (GeV^{2}) ; Events ",
                           100.,0.0,2.0);
auto* h_elecPolar = new TH1D("ele_polar","ele_polar",100, 0,2);
auto* h_gammaPolar = new TH1D("pimmult","pimmult",100, -1, 1);
auto* h_betamom = new TH2D("beta_vs_mom for pion","beta_vs_mom for pion",100, 0.0, 4.0,100,0.2,1.2);
auto* h_betamomkaon = new TH2D("beta_vs_mom for kaon","beta_vs_mom for kaon",100, 0.0, 4.0,100,0.2,1.2);
auto* h_betahadron = new TH2D("beta_vs_mom for hadron","beta_vs_mom for hadron",100, 0.0, 4.0,100,0.2,1.2);
auto* h_betadiff = new TH2D("beta-betapion","beta-betapion",100, 0.0, 4.0,100,-1.0,1.0);
auto* h_ECAL = new TH2D("h_ECALp_vs_mom","ECAL/p_vs_mom ;p(GeV/c);E_{tot}/p",100, 0.0, 8.0,100,0.0,0.5);
auto* h_Sampling_fraction = new TH1D("h_Sampling_fraction","Sampling Fraction; Sampling Fraction",100,0.1, 0.5);
auto* h_projectileY=new TH1D("h_projectileY","ECAL/p projectileY",100,0.0,0.5);
auto* h_vertexZ=new TH1D("VertexZ","Vz disribution for charge particles",200,-30,30);
auto* h_betacalcp = new TH2D(" beta_vs_mom "," betacalc_vs_mom;p(GeV/c);#beta ",200, 0.0, 4.0,200,0.2,1.2);
auto* h_beta_myParticle = new TH2D(" beta_vs_mom "," betamyParticle_vs_mom ",200, 0.0, 4.0,200,0.2,1.2);
auto* h_betaEBp = new TH2D(" beta_vs_mom "," betaEB_vs_mom ",200, 0.0, 4.0,200,0.2,1.2);
auto* h_betacalcn = new TH2D(" beta_vs_mom "," beta_vs_mom ",200, 0.0, 4.0,200,0.2,1.2);
auto* h_betamomCD = new TH2D(" beta_vs_mom CD "," beta_vs_mom CD;p;#beta ",200, 0.0, 4.0,200,0.2,1.2);
auto* h_pid=new TH1D("pid","pid",200,-500,500);
auto* h_pimult=new TH1D("h_pimult","pimult",200,0.0,10.0);
auto* h_electronmult=new TH1D("h_elemult","electronmult",200,0.0,10.0);
auto* h_electronmultFD=new TH1D("h_elemult","electronmultFD",200,0.0,10.0);
auto* h_electronmultiplicity=new TH1D("h_elemult","electron in FD but not trigger particles but ele_energy > 2.0",200,0.0,10.0);
auto* h_protonmult=new TH1D("h_protmult","protonmult",200,0.0,10.0);
auto* h_pippimmult = new TH2D("h_pippimmult","2D histogram for piplus mult and piminus mult",100, 0.0, 10.0,100,0.0,10.0);
auto* h_kdotg=new TH1D("h_kdotg","beam->Dot(Photon)",200,0.0,2);
auto* h_kpdotg=new TH1D("h_kpdotg","eprime->Dot(Photon)",200,0.0,2);
auto* h_kdotg2D=new TH2D("h_kpdotg2D"," beam->Dot(Photon) vs  eprime->Dot(Photon);eprime->Dot(Photon);beam->Dot(Photon)",200,0.0,2.0,200,0.0,2.0);
auto* h_photonenergy=new TH1D("h_photonenergy","photon energy",200,0,10.);


auto* h_Miss_total_fcut = new TH1D("h_Miss_total_fcut","Missing total  Mass^{2};  M_{X}^{2} (GeV^{2}) ; Events ",
                                   100.,-2.0,2.0);

auto* h_Miss_total_Energy_fcut = new TH1D("h_Miss_total_Energy_fcut","Missing total  Energy;  E_{Miss} (GeV) ; Events ",
                                          100.,-2.0,7.0);


auto* h_missing_momentum_fcut = new TH1D("h_Missing_momentum_fcut","Missing Momentum;  P_{X} (GeV) ; Events ",
                                         100.,-2.0,7.0);

//  auto* h_Miss_total_fcut = new TH1F("h_Miss_total_fcut","Missing total  Mass^{2};  M_{X}^{2} (GeV^{2}) ; Events ",
//         100.,-2.0,2.0);

auto* h_Pmiss_ey = new TH1D("h_Pmiss_ey","Pmiss->Dot(ey);  Pmiss.Dot(ey) (GeV) ; Events ",
                            100.,-2.0,7.0);

auto* h_Pz = new TH1D("h_Pz","Pz(along beam direction);  Pz (GeV) ; Events ",
                      100.,-2.0,7.0);

auto* h__a_b = new TH2D("h_a_b","b vs a; a;b ",  100.,-2.0,7.0,100,-2.0,7.0);


auto* h_Miss_pim_p_proton = new TH1D("h_Miss_pim_proton","Missing pi- Mass^{2};  M_{X=#pi^-}^{2} (GeV^{2}) ; Events ",
                                     100.,-2.0,2.0);


auto*h_Emiss_Pmiss_pip_cut= new TH2D("h_Emiss_Pmiss_cut","pi+ Emiss vs Pmiss after cuts; Pmiss_{#pi^{+}};Emiss_{#pi^{+}}",200,0.0,10.0,200,0.0,10.0);
auto*h_Emiss_Pmiss_pip= new TH2D("h_Emiss_Pmiss","pi+ Emiss vs Pmiss; Pmiss_{#pi^{+}};Emiss_{#pi^{+}}",200,0.0,10.0,200,0.0,10.0);

auto*  h_tof_ChiSq = new TH1D(" h_tof_ChiSq","Track to TOF1B match ChiSq; ChiSq; Events  ", 100.,0.0,10.0);

// auto*  h_chi2pid = new TH1F(" h_chi2pid","ChiSqpid; ChiSq pid; Events  ", 100.,-10.0,10.0);

auto*  h_PCAL_DOCA = new TH1D(" h_PCAL_DOCA","PCAL_DOCA; PCAL_DOCA(cm); Events  ", 100.,0.0,10.0);

auto*  h_Rec_Track_DC = new TH1D("h_REC_Track_DC_NDF","NDF; NDF; Events  ", 100.,0.0,50.0);
auto*  h_Rec_Track_DC_integral = new TH1D("h_REC_Track_DC_integral","Chi2_integral; ChiSq; Events  ", 100.,0.0,500.0);

auto*  h_Rec_Track_DC_Chi2 = new TH1D("h_REC_Track_DC_chi2_at_NDF30","Chi2; ChiSq; Events  ", 100.,0.0,300.0);
auto*  h_Rec_Track_DC_Chi2N = new TH1D("h_REC_Track_DC_chi2NDF","Chi2/NDF; ChiSq/NDF; Events  ", 100.,0.0,100.0);
auto*  h_Rec_Track_DC_Chi2N_integral = new TH1D("h_REC_Track_DC_chi2N_integral","Chi2/NDF_integral; ChiSq/NDF_integral; Events  ", 100.,0.0,500.0);

auto* h_st_vt=new TH1D("h_st_vt","startTime-vt",200,-10,10);

auto* h_Egamma=new TH1D("h_Egamma","gamma energy ElectronhighE.Dot(photon)< 0.05",200,0.0,5.0);

auto* h_DC_XY1 = new TH2D("h_DC_XY1","DC_XY1 electron before cut",200,-400, 400,200,-400,400);
auto* h_PCAL_XY = new TH2D("h_PCAL_XY","PCAL_XY only all electrons",200,-400, 400,200,-400,400);
auto* h_PCAL_XY_cut = new TH2D("h_PCAL_XY_cut","PCAL_XY_cut only for electron",200,-400, 400,200,-400,400);
auto* h_DC_XY1_cut = new TH2D("h_DC_XY1_cut","DC_XY1_after fiducial cut)",200,-400, 400,200,-400,400);

auto* h_DC_XY2 = new TH2D("h_DC_XY2","DC_XY2 electron before cut",200,-400, 400,200,-400,400);
auto* h_DC_XY2_cut = new TH2D("h_DC_XY2_cut","DC_XY2 after fiducial cut",200,-400, 400,200,-400,400);
auto* h_DC_XY3 = new TH2D("h_DC_XY3","DC_XY3 electron before cut",200,-400, 400,200,-400,400);
auto* h_DC_XY3_cut = new TH2D("h_DC_XY3_cut","DC_XY3 after fiducial cut",200,-400, 400,200,-400,400);
auto* h_E_PCAL = new TH2D("h_E_PCAL","ECAL vs PCAL; PCAL(GeV);ECAL(GeV)",200,0,0.25,200,0.0,0.35);
auto* h_ECinner_outer = new TH2D("h_ECinner_outer","EC_outer vs EC_inner;EC_{inner}(GeV);EC_{inner}(GeV)",200,0,0.7,200,0.,1.);
auto* h_DC_phi_theta = new TH2D("h_DC_phi_theta","DC phi theta;#phi(deg);#theta(deg)",200,-40, 40,200,0,60);
auto* h_PCAL_Lu=new TH1D("h_PCAL_Lu","PCAL U",200,0,500);
auto* h_PCAL_Lv=new TH1D("h_PCAL_Lv","PCAL V",200,0,500);
auto* h_PCAL_Lw=new TH1D("h_PCAL_Lw","PCAL W",200,0,500);
auto* h_3h_betaDiff = new TH1D("h_3h_betaDiff","#beta-#beta_{p} for protons",100,-0.5,0.5);
auto* h_proton_polar = new TH1D("h_proton_polar","#theta for protons",100,0.0,180.);
auto* h_pip_polar = new TH1D("h_pip_polar","#theta for pim",100,0.0,180.);
auto* h_pim_polar = new TH1D("h_pim_polar","#theta for pip",100,0.0,180.);
auto* h_dvz1 = new TH1D("h_dvz1","(eVz-pVz);e_{Vz}-p_{Vz}(cm);Events",200,-100.,100.);
auto* h_dvz2 = new TH1D("h_dvz2","(eVz-pi+Vz);e_{Vz}-pi+_{Vz}(cm);Events",200,-100.,100.);
auto* h_dvz3 = new TH1D("h_dvz3","(eVz-pi-Vz);e_{Vz}-pi-_{Vz}(cm);Events",200,-100.,100.);
auto* h_Nphe = new TH1D("h_nphe","Nphe; N_{photoelectron} X 10;Events",100,0.0,500.);
auto* h_Nphe_Etot = new TH2D("h_Nphe_Etot","h_Nphe_Etot;E_{tot}/p;N_{photoelectron} X 10",200,0.0, 0.4,200,0.0,500.);
auto* h_GEN_REC_p = new TH2D("h_GEN_REC_p","p_{GEN}/p_{REC} vs p_{REC};p_{REC}(GeV/c);p_{GEN}/p_{REC}",200,0.0,5.0,200,0.8,1.2);
auto* h_GEN_REC_theta = new TH2D("h_GEN_REC_p","p_{GEN}/p_{REC} vs #theta_{REC};#theta^{0};p_{GEN}/p_{REC}",200,0.0,60.0,200,0.8,1.2);
auto* h_GEN_REC_phi = new TH2D("h_GEN_REC_p","p_{GEN}/p_{REC} vs #phi_{REC};#phi^{0};p_{GEN}/p_{REC}",200,0.0,300.,200,0.8,1.2);

auto* h_CsHS_1=new TH1D("h_CsHS","CsHS;cos(#theta_{HS})",200,-1.0,1.0);
auto* h_phiHS=new TH1D("h_phiHS","phiHS;#phi_{HS})",200,-4.0,4.0);
auto* h_phiHS_1=new TH1D("h_phiHS","phiHS;#phi_{HS})",200,-200,200);
auto* h_phi=new TH1D("h_phi","phi;#phi",200,-200.,200.);
auto* h_eprimeDot1=new TH1D("h_eprimeDot1"," ;miss total energy(GeV^{2}) ",200,-2.,7.);
auto* h_eprimeDot2=new TH1D("h_eprimeDot2"," ;miss total energy(GeV^{2})",200,-2.,7.);
auto* h_eprimeDotP=new TH1D("h_eprimeDotP","beam->Dot(miss_total_P);beam->Dot(miss_total_P)",200,-2.,7.);
auto* ggIM_dist=new TH1D("h_ggIM","#gamma #gamma invariant Mass Sq;M_{#gamma #gamma}^{2} ",100,0.0,0.1);
auto* h_collinearity=new TH1D("h_collinearity","proton collinearity theta(deg);Events ",100,0.0,180);
auto* h_chi2pid_electron=new TH1D("h_chi2pid_electorn","chi2pid_elec;Events ",100,-7.0,7.0);
auto* h_chi2pid_proton=new TH1D("h_chi2pid_proton","chi2pid_proton;Events ",100,-7.0,7.0);
auto* h_chi2pid_piplus=new TH1D("h_chi2pid_piplus","chi2pid_piplus;Events ",100,-7.0,7.0);
auto* h_chi2pid_piminus=new TH1D("h_chi2pid_piminus","chi2pid_piminus;Events ",100,-7.0,7.0);
auto* h_missingmass_energy = new TH2D("h_missingmass_energy","missing_energy vs missing_mass;missing_energy;missing_mass",100,-2.0,2.0,100,-2.0,2.0);
auto* h_Pymiss_fcut=new TH1D("h_Pymiss_fcut","yq.Dot(P_{miss};yq.Dot(P_{miss})(GeV);Events ",200,-10,10);
auto* h_Pymiss_cut=new TH1D("h_Pymiss_cut","yq.Dot(P_{miss} after cuts);yq.Dot(P_{miss})(GeV);Events ",200,-10,10);
auto* h_Miss_total_momentum=new TH1D("h_Miss_total_momentum","P_{miss};P_{miss}(GeV);Events ",100,0.0,20.0);
auto* h_miss_doublepion=new TH1D("h_miss_doublepion",";(GeV);Events ",100,0.0,2.0);
auto* h_CsHs_phiHs = new TH2D("h_CsHs_phiHs","phi_{pi+Rest} vs cos(th_pi+Rest); cos(#theta_{#pi^{+}_rest});#phi_{#pi^{+}_rest}",200,-1.0, 1.0,200,-4.0,4.0);
auto* h_CsHS_M12 = new TH2D("h_CsHs_M12"," cos(th_pi+Rest) vs M(#pi,#pi);M_{#pi #pi}; cos(#theta_{#pi^{+}_rest})",200,0.0, 1.4,200,-1.0,1.0);
auto* h_CsHS=new TH1D("h_CsHS","CsHS;cos(#theta_{HS})",200,-1.0,1.0);
auto* h_Dalitz_1=new TH2D("h_Dalitz_1"," M(p #pi^{+}) vs M(#pi^{+},#pi^{-});M_{#pi^{+} #pi^{-}};M(p #pi^{+}))",200,0.0, 1.4,200,1.0,3.0);
auto* h_Dalitz_2=new TH2D("h_Dalitz_2"," M(p #pi^{-}) vs M(#pi^{+},#pi^{-});M_{#pi^{+} #pi^{-}};M(p #pi^{-}))",200,0.0, 1.4,200,1.0,3.0);
auto* h_CsHS_Mppip = new TH2D("h_CsHs_Mppip"," cos(th_pi+Rest) vs M(p,#pi^{+});M_{p #pi^{+}}; cos(#theta_{#pi^{+}_rest})",200,1.0, 3.0,200,-1.0,1.0);
auto* h_CsHS_Mppim = new TH2D("h_CsHs_Mppim"," cos(th_pi+Rest) vs M(p,#pi^{-});M_{p #pi^{-}}; cos(#theta_{#pi^{+}_rest})",200,1.0, 3.0,200,-1.0,1.0);
auto* h_CsHs_pi=new TH2D("h_CsHs_pi","cos(#theta_{HS}) vs p_{#pi^{+}};p_{#pi^{+}};cos(#theta_{HS})",200,0.0,5.0,200,-1.0,1.0);
auto* h_y_vector=new TH2D("h_y_vector","y vs xBj;xBj;y",200,0.0,1.0,200,0.0,2.0);
auto* h_phi_M12 = new TH2D("h_phi_M12"," #phi vs M(#pi,#pi);M_{#pi #pi}; #phi",200,0.0, 1.4,200,-200,200);
auto* h_phiH_M12 = new TH2D("h_phiH_M12"," phiH vs M(#pi,#pi);M_{#pi #pi};phiH",200,0.0, 1.4,200,-200,200);
auto* h_phithetaCD = new TH2D("h_phithetaCD"," phi vs theta for CD (charge>0);#phi^{0};#theta^{0}",200,-200,200,200,0.0,90.);
auto* h_CsHS_t1=new TH1D("h_CsHS_t1","CsHS;cos(#theta_{HS})",200,-1.0,1.0);
auto* h_CsHS_t2=new TH1D("h_CsHS_t2","CsHS;cos(#theta_{HS})",200,-1.0,1.0);
auto* h_CsHS_t3=new TH1D("h_CsHS_t3","CsHS;cos(#theta_{HS})",200,-1.0,1.0);
auto* h_phiHS_t1=new TH1D("h_phiHS_t1","phiHS;#phi_{HS})",200,-200,200);
auto* h_phiHS_t2=new TH1D("h_phiHS_t2","phiHS;#phi_{HS})",200,-200,200);
auto* h_phiHS_t3=new TH1D("h_phiHS_t3","phiHS;#phi_{HS})",200,-200,200);
auto* h_beta_calorimeter = new TH2D("h_beta_calorimeter"," beta vs p;p(GeV/c);#beta",200,0.0,6.0,200,0.6,1.3);
auto* h_pn=new TH1D("h_pn","neutron momentum;p_{n}(GeV/c))",200,1.0,6.0);
auto* h_proton_piplus_no_cuts_1 = new TH1D("h_proton_piplus_no_cuts","(Proton+piplus) Invariant Mass in cos(pi+_rest) < -0.8;  M_{p#pi^{+}} (GeV) ; Events ",
                                           200.,1.0,3.0);
auto* h_proton_piplus_no_cuts_2 = new TH1D("h_proton_piplus_no_cuts","(Proton+piplus) Invariant Mass in cos(pi+_rest) > 0.0;  M_{p#pi^{+}} (GeV) ; Events ",
                                           200.,1.0,3.0);

auto* h_proton_piminus_no_cuts_1 = new TH1D("h_proton_piminus_no_cuts","(Proton+piminus) Invariant Mass in cos(pi+_rest) > 0.8;  M_{p#pi^{+}} (GeV) ; Events ",
                                            200.,1.0,3.0);
auto* h_proton_piminus_no_cuts_2 = new TH1D("h_proton_piminus_no_cuts","(Proton+piminus) Invariant Mass in cos(pi+_rest) < 0.0;  M_{p#pi^{+}} (GeV) ; Events ",
                                            200.,1.0,3.0);

auto* h_pmag_DeltaT_FTOF_positive = new TH2D(" h_pmag_DeltaT_FTOF_positive "," #Delta t vs p for positive in FTOF;p(GeV);#Delta t ",200, 0.0, 3.0,200,-10.0,10.0);
auto* h_pmag_DeltaT_FTOF_negative = new TH2D(" h_pmag_DeltaT_FTOF_negative "," #Delta t vs p for negative in FTOF;p(GeV);#Delta t ",200, 0.0, 3.0,200,-10.0,10.0);
auto* h_pmag_DeltaT_CTOF_positive = new TH2D(" h_pmag_DeltaT_CTOF_positive "," #Delta t vs p for positive in CTOF;p(GeV);#Delta t ",200, 0.0, 3.0,200,-10.0,10.0);
auto* h_pmag_DeltaT_CTOF_negative = new TH2D(" h_pmag_DeltaT_CTOF_negative "," #Delta t vs p for negative in CTOF;p(GeV);#Delta t ",200, 0.0, 3.0,200,-10.0,10.0);

auto* h_pmag_DeltaT_FTOF_proton = new TH2D(" h_pmag_DeltaT_Ftof_proton "," #Delta t vs p for proton in FTOF1B;p(GeV);#Delta t ",200, 0.0, 3.0,200,-6.0,6.0);
auto* h_pmag_DeltaT_FTOF_piplus = new TH2D(" h_pmag_DeltaT_Ftof_piplus "," #Delta t vs p for piplus in FTOF1B;p(GeV) ;#Delta t ",200, 0.0, 3.0,200,-6.0,6.0);
auto* h_pmag_DeltaT_FTOF_piminus = new TH2D(" h_pmag_DeltaT_Ftof_piminus "," #Delta t vs p for piminus in FTOF1B;p(GeV) ;#Delta t ",200, 0.0, 3.0,200,-6.0,6.0);
auto* h_pmag_DeltaT_CTOF_proton = new TH2D(" h_pmag_DeltaT_Ctof_proton "," #Delta t vs p for proton in CTOF;p(GeV) ;#Delta t ",200, 0.0, 3.0,200,-6.0,6.0);
auto* h_pmag_DeltaT_CTOF_piplus = new TH2D(" h_pmag_DeltaT_Ctof_piplus "," #Delta t vs p for piplus in CTOF;p(GeV) ;#Delta t ",200, 0.0, 3.0,200,-6.0,6.0);
auto* h_pmag_DeltaT_CTOF_piminus = new TH2D(" h_pmag_DeltaT_Ctof_piminus "," #Delta t vs p for piminus in CTOF;p(GeV) ;#Delta t ",200, 0.0, 3.0,200,-6.0,6.0);


auto* h_pion_missing = new TH1D("h_pion_missing","Two-pion MM2;  #pi^{+} #pi^{-} (GeV^{2}) ; Events ",
                                200.,0.5,1.5);
auto* h_proton_piplus_missing = new TH1D("h_proton_piplus_missing","(Proton+piplus) MM2;  p#pi^{+} (GeV^{2}) ; Events ",
                                         200.,-0.2,0.4);
auto* h_proton_piminus_missing = new TH1D("h_proton_piminus_missing","(Proton+piminus) MM2;  p#pi^{-} (GeV^{2}) ; Events ",
                                          200.,-0.2,0.4);
auto * h_MX2_fit = new TH1D("h_MX2_fit","M_{X}^2 After Kin. Fit; M_{X}^2 (GeV^{2}); counts  ",
                            200, -0.50,1.5);
auto * h_EX_fit = new TH1D("h_EX_fit","E_{X} After Kin. Fit; E_{X} (GeV); counts  ",
                           200, -0.5,1.5);
auto * h_PX_Y_fit = new TH1D("h_PX_Y_fit","P_{X} out of plane After Kin. Fit; P_{X}^{y} (GeV); counts  ",
                             200, -1.0,1.0);
auto * h_MXP2_fit = new TH1D("h_MXP2_fit","M_{p+X}^2 After Kin. Fit; M_{p+X}^2 (GeV^{2}); counts  ",
                             200, -0.50,1.5);
auto * h_MXpip2_fit = new TH1D("h_MXpip2_fit","M_{pi+,X}^2 After Kin. Fit; M_{pi+,X}^2 (GeV^{2}); counts  ",
                               200, -1.00,1.0);
auto * h_MXpim2_fit = new TH1D("h_MXpim2_fit","M_{pi-,X}^2 After Kin. Fit; M_{pi-,X}^2 (GeV^{2}); counts  ",
                               200, -1.00,1.0);
auto *h_func_before_fit=new TH1D("h_func_before","Func before fit",200,-2.0,2.0);
auto *h_func_after_fit=new TH1D("h_func_after","Func after fit",200,-2.0,2.0);
auto *h_func_before_after_fit=new TH1D("h_func_differrnce","fcn(after) - fcn(initial)",200,-2.0,2.0);

#endif /* KinematicFitter_hist_hpp */
