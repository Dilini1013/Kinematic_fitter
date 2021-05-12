//
//  KinematicFitter_code_11May.cpp
//
//
//  Created by Dilini Bulumulla on 5/11/21.
//

//
//  KinematicFitter_code.cpp
//
//
//  Created by Dilini Bulumulla on 4/25/21.
//


#include "KinematicFitter_code.hpp"
#include "KinematicFitter_hist.hpp"

void KinematicFitter_code_2(){
    //void main(){
    //    double *par; // not relevant is this scope
    
    FILE *outfile;
    outfile = fopen("deeprad.dat", "w");  // using relative path name of file
    if (outfile == NULL) {
        printf("Unable to open file.");
    }
    
    TFile *AnaFile = new TFile("kinematic1.root", "RECREATE", "Test");
    // TTree *t1 = new TTree("t1","var tree");
    Double_t proton_missSq, pip_missSq,pim_missSq;
    
    //  t1->Branch("proton_missSq",&proton_missSq,"proton_missSq/D");
    // t1->Branch("pip_missSq",&pip_missSq,"pip_missSq/D");
    //  t1->Branch("pim_missSq",&pim_missSq,"pim_missSq/D");
    
    //  TTree *AnaTree = (TTree*)AnaFile->Get("tree");
    
    
    
    TNtuple *contNTuple = new TNtuple("contNTuple","contNTuple","hmiss:etime:Ex:Ey:Ez");
    // Record start time
    auto start = std::chrono::high_resolution_clock::now();
    int counter=0;
    int coinc_counter = 0;
    
    /////////////////////////////////////
    // std::string inputFile =
    bool simulation=false;
    bool data=false;
    
    //    cout<<"Analysing hipo file "<<inputFile<<endl;
    
    TChain fake("hipo");
    fake.Add("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_00503*.hipo");
    
    
    
    auto files=fake.GetListOfFiles();
    TRandom3 ran3;
    ran3.SetSeed(1357911);
    //some particles
    auto db=TDatabasePDG::Instance();
    /* // move to global scope
     TLorentzVector beam(0,0,10.6,10.6);
     TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
     */
    beam.SetPxPyPzE(0.0, 0.0, 10.6, 10.6);
    target.SetPxPyPzE(0.0, 0.0, 0.0, db->GetParticle(2212)->Mass());
    TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
    TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass());
    TLorentzVector g1(0,0,0,0);
    TLorentzVector g2(0,0,0,0);
    TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass());
    TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass());
    TLorentzVector kaonp(0,0,0,db->GetParticle(321)->Mass());
    TLorentzVector kaonm(0,0,0,db->GetParticle(-321)->Mass());
    TLorentzVector Delta;
    TLorentzVector P122;
    TLorentzVector W4vec(0,0,0,0);
    TLorentzVector q4vec(0,0,0,0);
    TLorentzVector PpiplusRest(0,0,0,db->GetParticle(211)->Mass());
    TLorentzVector Miss_proton;
    //    TLorentzVector electrontrigger,piplus,piminus,proton;
    //  remove electron from list
    TLorentzVector Lambda,M2K,M122_calctrigger,M122_calchighE,photon,Miss_pim, Miss_pip,electronhighE,photonBrems, M122_calctriggerwithBremsphoton, Miss_total,pionHighE,miss_doublepion,proton_piplus,proton_piminus,pion_invariant;
    TLorentzVector protonFD, protonCD,piplusFD,piplusCD,piminusFD,piminusCD;
    TLorentzVector electron_fit, proton_fit, piplus_fit, piminus_fit;
    TVector3 beta_rest,Xhadron,Yhadron,Zhadron,P3piplusRest,ele3,eprime3,piplus3,piminus3,yq,zq,xq,proton3_miss,proton3_measured;
    TVector3 Xsigma, Ysigma, Zsigma;
    TLorentzVector PX_fit, PXP_fit, PXpip_fit, PXpim_fit;
    TLorentzVector P4X,P4XP,P4Xpip,P4Xpim;
    
    double weight_function;
    
    double M122,phiH,phii,csphi,sinphi;
    double A_par=2.0;
    double B_par=2.0;
    double tau=3.0;
    
    double nu, xB,EpCM,EpprimeCM,q0CM,qCM,PCM,tmin,tmax,P12CM,E12CM,CsHS,y_vector;
    
    
    double px_e,py_e,pz_e,px_p,py_p,pz_p,px_pip,py_pip,pz_pip,px_pim,py_pim,pz_pim;
    double px_e_fit, py_e_fit, pz_e_fit, px_p_fit,py_p_fit,pz_p_fit, px_pip_fit,py_pip_fit, pz_pip_fit, px_pim_fit, py_pim_fit, pz_pim_fit;
    double M1,M2,M3;
    double M_X2_fit,E_X_fit;
    TVector3 eprime3_fit,yq_fit;
    TLorentzVector Miss_total_fit;
    
    
    
    unsigned long int piplusmult=0;
    unsigned long int piminusmult=0;
    int electronmult=0;
    unsigned long int protonmult=0;
    int photonmult=0;
    double eprimevz,pprimevz;
    double X1=0,Y1=0, X2=0,Y2=0, X3=0,Y3=0;
    double sum1=0.0;
    double sum2=0.0;
    double protonBeta=0.0;
    double beta_PCAL;
    
    
    
    gBenchmark->Start("timer");
    
    
    
    
    int nsize=0;
    int nelec=0;
    int pimult=0;
    double beta_calc=0.0;
    double px=0.0,py=0.0, pz=0.0,pmag=0.0 ;
    double theta=0.0;
    double pmagnitude=0.0;
    double phi=0.0;
    double WSq=0.0,QSq=0.0;
    double REC_Track_DC,Rec_Track_DC_Chi2,Rec_Track_DC_Chi2N;
    double PCAL_lu, PCAL_lv, PCAL_lw, PCAL_sector;
    unsigned long int multiplex = 100;
    
    
    const double vertex_cut_neg = -20.0; // cm
    const double vertex_cut_pos =  20.0; // cm
    const double Track_DC_Chi2Max = 20.0;
    double f;
    
    
    
    TMinuit *gMinuit = new TMinuit(nparMINUIT); // only 12 free params
    
    gMinuit->SetFCN(fcn);
    //  Define all free parameters in MINUIT with dummy values;
    //double Value0 = 1.0;
    
    
    
    
    
    
    for(Int_t i=0;i<files->GetEntries();i++){
        //create the event reader
        clas12reader c12(files->At(i)->GetTitle());
        //  clas12reader c12(files->At(i)->GetTitle(),{0});//add tags {tag1,tag2.getg,tag3,...}
        //   auto particles = c12.getDetParticles();
        //Add some event Pid based selections
        //////////c12.AddAtLeastPid(211,1); //at least 1 pi+
        //c12.addExactPid(11,1);    //exactly 1 electron
        //c12.addExactPid(211,1);    //exactly 1 pi+
        //c12.addExactPid(-211,1);    //exactly 1 pi-
        //c12.addExactPid(2212,1);    //exactly 1 proton
        //c12.addExactPid(22,2);    //exactly 2 gamma
        //////c12.addZeroOfRestPid();  //nothing else
        //////c12.useFTBased(); //and use the Pids from RECFT
        double eprimeVz, betapion,piplus1,piplus2,pionPiD;
        
        
        while(c12.next()==true){
            
            
            
            //particle start time
            double startTime= c12.event()->getStartTime();
            auto particles = c12.getDetParticles();
            int runconfig_event=c12.runconfig()->getEvent();
            
            
            
            
            // reset all 4-vectors to zero
            electrontrigger.SetPxPyPzE(0.,0.,0.,0.);
            proton.SetPxPyPzE(0.,0.,0.,0.);
            piplus.SetPxPyPzE(0.,0.,0.,0.);
            piminus.SetPxPyPzE(0.,0.,0.,0.);
            
            int pielect=-1;
            double piMax=0.0;
            int ielect = -1;
            double EeMax=0.0;
            nsize = 0;
            nelec=0;
            piminusmult = 0;
            piplusmult  = 0;
            electronmult=0;
            protonmult=0;
            photonmult=0;
            coinc_counter=0; // CEH
            double xBj=0.0;
            double electronVz=0.0;
            double protonVz=0.0;
            double piplusVz=0.0;
            double piminusVz=0.0;
            bool electvalid=false;
            bool electtrig=false;
            bool protonPid=false;
            bool piplusPid=false;
            bool piminusPid=false;
            bool photonPid=false;
            bool protonCDPid=false;
            bool piplusCDPid=false;
            bool piminusCDPid=false;
            bool photonBremsvalid=false;
            double smear_px = 0.0;
            double smear_py = 0.0;
            double smear_pz = 0.0;
            double  Rec_Track_DC=0.0;
            double Rec_Track_DC_Chi2=0.0;
            double Rec_Track_DC_Chi2N=0.0;
            double PCal_Lu=0.0,  PCal_Lv=0.0, PCal_Lw=0.0;
            double rth=0.0, XS1=0.0,YS1=0.0,phiS1=0.0,r_cut=40, XS2=0.0,YS2=0.0,phiS2=0.0, XS3=0.0,YS3=0.0,phiS3=0.0, PCAL_X=0.0, PCAL_Y=0.0;
            int status;
            unsigned long int i_protons = 0;
            unsigned long int i_piplus = 0;
            unsigned long int i_piminus = 0;
            
            
            
            
            for(int i=0;i<(int)particles.size();i++){
                
                int pid = particles[i]->par()->getPid();
                double px = particles[i]->par()->getPx();
                double py = particles[i]->par()->getPy();
                pz = particles[i]->par()->getPz();
                double pmag  = particles[i]->par()->getP();
                double beta = particles[i]->par()->getBeta();
                short charge= particles[i]->par()->getCharge();
                double vx =particles[i]->par()->getVx();
                double vy = particles[i]->par()->getVy();
                double vz = particles[i]->par()->getVz();
                double vt=particles[i]->par()->getVt();
                status= particles[i]->par()->getStatus();
                double chi2pid=particles[i]->par()->getChi2Pid();
                
                //calculating pmagnitude,phi and theta using px,py and pz
                pmagnitude = sqrt(px*px+py*py+pz*pz);
                theta= acos(pz/pmagnitude);
                phi=TMath::ATan2(py, px);
                
                if(particles[i]->traj(DC,36)->getDetector()){
                    Rec_Track_DC=particles[i]->trk(DC)->getNDF();
                    Rec_Track_DC_Chi2=particles[i]->trk(DC)->getChi2();
                    Rec_Track_DC_Chi2N=particles[i]->trk(DC)->getChi2N();
                    h_Rec_Track_DC_Chi2N->Fill(Rec_Track_DC_Chi2N);
                    // printf(" nnmnmn %f \n",Rec_Track_DC );
                }
                
                
                
                if(!db->GetParticle(pid)) continue;
                nsize++;
                
                
                
                bool Z_vertex_cut=z_vertex_cut(particles[i]);
                bool EC_sampling_fraction= EC_sampling_fraction_cut(particles[i]);
                bool minimal_PCAL= minimal_PCAL_energy_deposition(particles[i]);
                bool PCAL_homogeneous= EC_hit_position_fiducial_cut_homogeneous(particles[i]); //here I use the i particle
                bool PCAL_hit_fiducial= EC_hit_position_fiducial_cut(particles[i]);
                int region[]={1,2,3};
                bool DC_fiducial_region1= DC_hit_position_counts_fiducial_cut(particles[i], region[0]);
                bool DC_fiducial_region2= DC_hit_position_counts_fiducial_cut(particles[i], region[1]);
                bool DC_fiducial_region3= DC_hit_position_counts_fiducial_cut(particles[i], region[2]);
                
                
                if(particles[0]->par()->getPid()==11 &&  status > -3000 && status < -2000 && !electtrig){
                    // only keep first valid trigger electron
                    electronVz=particles[i]->par()->getVz();
                    h_chi2pid_electron->Fill(particles[i]->par()->getChi2Pid());
                    if(PCAL_homogeneous){
                        if(DC_fiducial_region1 && DC_fiducial_region2 && DC_fiducial_region3){
                            if(EC_sampling_fraction){
                                if( minimal_PCAL){
                                    if(Z_vertex_cut){
                                        if((particles[i]->par()->getChi2Pid()) > -2.5 && (particles[i]->par()->getChi2Pid()) < 2.5 ){
                                            px_e=px;
                                            py_e=py;
                                            pz_e=pz;
                                            
                                            electrontrigger.SetXYZM(px_e,py_e,pz_e,Me);
                                            
                                            electtrig=true;
                                            electronmult++;
                                            //  Electron Trigger will be passed to fcn as a global variable
                                            
                                            gMinuit->DefineParameter(0,"px_e_fit", px, px*err_init,0.0,0.0);
                                            gMinuit->DefineParameter(1, "py_e_fit", py, py*err_init,0.0,0.0);
                                            gMinuit->DefineParameter(2, "pz_e_fit", pz, pz*err_init,0.0,0.0);
                                            
                                            
                                        }
                                    }
                                }
                            }
                        }
                    }
                }  else if(pid==2212){
                    protonVz=particles[i]->par()->getVz();
                    h_dvz1->Fill(electronVz-protonVz);
                    h_chi2pid_proton->Fill(particles[i]->par()->getChi2Pid());
                    if((electronVz-protonVz)> vertex_cut_neg && (electronVz-protonVz)<vertex_cut_pos ){
                        if((particles[i]->par()->getChi2Pid()) > -2.5 && (particles[i]->par()->getChi2Pid()) < 2.5 ){
                            if(particles[i]->getRegion()==FD){
                                if(DC_fiducial_region1 && DC_fiducial_region2 && DC_fiducial_region3){
                                    px_p=px;
                                    py_p=py;
                                    pz_p=pz;
                                    proton.SetXYZM(px_p, py_p, pz_p, Mp);
                                    protonPid=true;
                                    i_protons *= multiplex;
                                    i_protons += i;
                                    protonmult++;  //CEH
                                }
                            }
                        }  //  end of Vz cut if statement
                    }//Chi2Pid cut for proton
                    // end of if pid==2212
                } else if(pid==211){
                    piplusVz=particles[i]->par()->getVz();
                    h_dvz2->Fill(electronVz-piplusVz);
                    h_chi2pid_piplus->Fill(particles[i]->par()->getChi2Pid());
                    if((electronVz-piplusVz)> vertex_cut_neg && (electronVz-piplusVz)<vertex_cut_pos ){
                        if((particles[i]->par()->getChi2Pid()) > -2.64 && (particles[i]->par()->getChi2Pid()) < 2.64 ){
                            if(particles[i]->getRegion()==FD){
                                if(DC_fiducial_region1 && DC_fiducial_region2 && DC_fiducial_region3){
                                    
                                    px_pip=px;
                                    py_pip=py;
                                    pz_pip=pz;
                                    piplus.SetXYZM(px_pip,py_pip,pz_pip,Mpi);
                                    piplusPid=true;
                                    i_piplus *= multiplex;
                                    i_piplus += i;
                                    
                                    piplusmult++;
                                }
                            }
                        };
                    }  //  end of Vz cut if statement
                    // end of if pid==211
                } else  if(pid==-211){
                    piminusVz=particles[i]->par()->getVz();
                    h_dvz3->Fill(electronVz-piminusVz);
                    h_chi2pid_piminus->Fill(particles[i]->par()->getChi2Pid());
                    if((electronVz-piminusVz)> vertex_cut_neg && (electronVz-piminusVz)<vertex_cut_pos ){
                        if(particles[i]->getRegion()==FD){
                            if(DC_fiducial_region1 && DC_fiducial_region2 && DC_fiducial_region3){
                                if((particles[i]->par()->getChi2Pid()) > -2.79 && (particles[i]->par()->getChi2Pid()) < 2.79 ){
                                    //DC Fiducial cuts valid
                                    px_pim=px;
                                    py_pim=py;
                                    pz_pim=pz;
                                    piminus.SetXYZM(px_pim,py_pim,pz_pim,Mpi);
                                    piminusPid=true;
                                    i_piminus *= multiplex;
                                    i_piminus += i;
                                    
                                    piminusmult++;
                                }
                            }
                        }
                    }  //  end of Vz cut if statement
                    // end of if pid==-211
                } else if(pid==22){
                    if( PCAL_hit_fiducial){ //PCAL Fiducial cuts
                        if(beta_PCAL < 1.1){
                            h_beta_calorimeter->Fill(pmag,beta);
                            photon.SetPxPyPzE(px,py,pz,sqrt(pow(pmag,2) + pow(Mphoton,2)));
                            photonPid=true;
                        }
                    }
                }
                
            } // end of particle loop in event
            
            if(electronmult*protonmult*piplusmult*piminusmult==0) continue; //  CEH
            
            
            // CEH
            // use i_proton, i_piplus, i_piminus to pick out candidate particles
            // and find all H(e,e'p pi+ pi-)X combinatios in the event
            //  Move all the gMinuit->DefineParameter statements here,
            //  because we are building multiple p(e,e' p pi pi) candidates for each event.
            //  Apply mMINUIT(Migrad) ONLY to candidate exclusive events
            //  We only kept first electron, so electron variables defined in previous loop.
            piMax=0.0;
            int jndex;
            int ii;
            for (int iPart=0; iPart<(int)protonmult; iPart++)
            {
                jndex =(int)(i_protons%multiplex);
                if(jndex==0) {
                    printf("piminus error particle# = %i \n", jndex);
                    continue;
                }
                //    printf("jndex %i, i_protons %lu, mprotons %i \n",jndex, i_protons, protonmult);
                
                px = particles[jndex]->par()->getPx();
                py = particles[jndex]->par()->getPy();
                pz = particles[jndex]->par()->getPz();
                pmag = particles[jndex]->par()->getP();
                proton.SetPxPyPzE(px,py,pz,sqrt(pmag*pmag+Mp*Mp));
                /*
                 //  FCN will now access global 4-vector proton directly
                 //  Values of "..._fit" will be  initialized in FCN
                 //  There was an error here in py_p_fit and pz_p_fit
                 */
                gMinuit->DefineParameter(3,"px_p_fit", px, px*err_init,0.0,0.0);
                gMinuit->DefineParameter(4,"py_p_fit", py, py*err_init,0.0,0.0);
                gMinuit->DefineParameter(5,"pz_p_fit", pz, pz*err_init,0.0,0.0);
                
                for(int jPart=0; jPart<(int)piplusmult; jPart++)
                {
                    jndex =(int) i_piplus%multiplex;
                    if(jndex==0) continue;
                    px = particles[jndex]->par()->getPx();
                    py = particles[jndex]->par()->getPy();
                    pz = particles[jndex]->par()->getPz();
                    pmag = particles[jndex]->par()->getP();
                    piplus.SetPxPyPzE(px,py,pz,sqrt(pmag*pmag+Mpi*Mpi));
                    // PiPlus parameters
                    
                    gMinuit->DefineParameter(6,"px_pip_fit", px, px*err_init,0.0,0.0);
                    gMinuit->DefineParameter(7, "py_pip_fit", py, py*err_init,0.0,0.0);
                    gMinuit->DefineParameter(8, "pz_pip_fit", pz, pz*err_init,0.0,0.0);
                    for(int kPart=0; kPart<(int)piminusmult; kPart++)
                    {
                        jndex =(int) i_piminus%multiplex;
                        if(jndex==0) continue;
                        px = particles[jndex]->par()->getPx();
                        py = particles[jndex]->par()->getPy();
                        pz = particles[jndex]->par()->getPz();
                        pmag = particles[jndex]->par()->getP();
                        piminus.SetPxPyPzE(px,py,pz,sqrt(pmag*pmag+Mpi*Mpi));
                        // Now global initialization done in fcn
                        
                        gMinuit->DefineParameter(9,"px_pim_fit", px, px*err_init,0.0,0.0);
                        gMinuit->DefineParameter(10,"py_pim_fit",py, py*err_init,0.0,0.0);
                        gMinuit->DefineParameter(11,"pz_pim_fit",pz, pz*err_init,0.0,0.0);
                        
                        
                        q4vec=beam-electrontrigger;
                        // Delta=Q2-P122;
                        // Delta=proton-target; // Move inside loop
                        W4vec=beam+target-electrontrigger;   ///  I Hate this 4-vector name  Should be PTotal or W4vec or some such
                        WSq=W4vec.M2();
                        QSq=-q4vec.M2();
                        hQ2->Fill(QSq);
                        xBj=QSq/(2*target.Dot(q4vec));
                        y_vector=(q4vec).Dot(target)/beam.Dot(target);
                        
                        
                        h_Q2_xBj_cut->Fill(xBj,QSq);
                        
                        //  Now we have one unique combination of e, p, pi+, pi- that pass fiducial cuts
                        Delta=proton-target;
                        Miss_total=beam+target-electrontrigger-proton-piminus-piplus;
                        Miss_proton =beam+target-electrontrigger-piminus-piplus;
                        h_Miss_proton_fcut->Fill(Miss_proton.M2());
                        
                        Miss_pim=beam+target-electrontrigger-proton-piplus;
                        h_Miss_pim_fcut-> Fill(Miss_pim.M2());
                        
                        Miss_pip=beam+target-electrontrigger-proton-piminus;
                        h_Miss_pip_fcut->Fill(Miss_pip.M2());
                        
                        h_missing_momentum_fcut->Fill(Miss_total.P());
                        
                        h_proton_polar->Fill(proton.Theta()*180./TMath::Pi());
                        h_pim_polar->Fill(piminus.Theta()*180./TMath::Pi());
                        h_pip_polar->Fill(piplus.Theta()*180./TMath::Pi());
                        
                        h_Miss_total_fcut->Fill(Miss_total.M2());
                        h_Miss_total_Energy_fcut->Fill(Miss_total.E());
                        ele3=beam.Vect();
                        //  Not valid in collider kinematics
                        eprime3=electrontrigger.Vect();
                        yq=ele3.Cross(eprime3);
                        yq.Unit();
                        
                        h_Pymiss_fcut->Fill(yq.Dot(Miss_total.Vect()));
                        pion_invariant=piplus+piminus;
                        proton_piplus= proton+piplus;
                        proton_piminus=proton+piminus;
                        double protonpipMassSq= 2.* proton.Dot(piplus)+ (Mp*Mp)+ (Mpi*Mpi);
                        double protonpimMassSq= 2.* proton.Dot(piminus)+ (Mp*Mp)+ (Mpi*Mpi);
                        
                        
                        
                        
                        M122=P122.M2();
                        h_MpipiSq_fcut->Fill(M122);
                        //  Now apply all exclusivity cuts
                        //  Replace cut limits with constants!!!
                        if(QSq > 1.0 && WSq > 4.0){
                            h_pion_missing->Fill(pion_invariant.M());
                            h_proton_piplus_missing->Fill(proton_piplus.M2());
                            h_proton_piminus_missing->Fill(proton_piminus.M2());
                            h_Pymiss_cut->Fill(yq.Dot(Miss_total.Vect()));
                            if(Miss_total.E() > -0.2 && Miss_total.E() < 0.2){
                                if(Miss_total.M2() > -0.1 && Miss_total.M2() < 0.1){
                                    if(yq.Dot(Miss_total.Vect())> -1 && yq.Dot(Miss_total.Vect())< 1){
                                        if(Rec_Track_DC_Chi2N < Track_DC_Chi2Max){
                                            //  Now we have a candidate exclusive event
                                            // Make histograms, Apply Kinematic Fit, then make revised histograms
                                            // Make initial histograms
                                            proton_missSq=Miss_proton.M2();
                                            pip_missSq=Miss_pip.M2();
                                            pim_missSq=Miss_pim.M2();
                                            //   t1->Fill();
                                            h_Miss_proton_all_cuts->Fill(Miss_proton.M2());
                                            h_Miss_pip_all_cuts->Fill(Miss_pip.M2());
                                            h_Miss_pim_all_cuts->Fill(Miss_pim.M2());
                                            h_pion_invariant_all_cuts->Fill(pion_invariant.M());
                                            
                                            //Mass(proton,pi+) and Mass(proton,pi-) in sigma region m_pi pi < 0.6
                                            if(pion_invariant.M() < 0.6){
                                                h_proton_piplus_IM_all_cuts_sigmaR->Fill(sqrt(protonpipMassSq));
                                                h_proton_piminus_IM_all_cuts_sigmaR->Fill(sqrt(protonpimMassSq));
                                            }
                                            //Mass(proton,pi+) and Mass(proton,pi-) in rho region 0.67 < m_pi pi <  0.87
                                            if(pion_invariant.M() > 0.67 && pion_invariant.M() < 0.87){
                                                h_proton_piplus_IM_all_cuts_rhoR->Fill(sqrt(protonpipMassSq));
                                                h_proton_piminus_IM_all_cuts_rhoR->Fill(sqrt(protonpimMassSq));
                                            }
                                            
                                            //  Apply Kinematic Fit
                                            //  const int npar = gMinuit->GetNumPars();
                                            // Double_t grad[npar], par[npar];
                                            // Int_t iflag = 1;
                                            // Double_t fval;
                                            // gMinuit->Eval(npar,grad, fval,par, iflag);
                                            /*
                                             for (int ipar=0; ipar<nparMINUIT; ipar++) {
                                             printf(" par[%02d] = %8.6f \n", gMINUIT->GetParameter(ipar));
                                             }
                                             
                                             printf (" px_e %8.6f \n", px_e);
                                             printf (" py_e %8.6f \n", py_e);
                                             printf (" pz_e %8.6f \n", pz_e);
                                             printf (" Me %8.6f \n", Me);
                                             printf ("  sigma_p_electron %8.6f \n", sigma_p_electron);
                                             printf (" sigma_costheta_electron %8.6f \n", sigma_theta_electron);
                                             printf (" sigma_phi_electron %8.6f \n", sigma_phi_electron);
                                             printf (" px_e_fit  %8.6f \n", px_e_fit);
                                             printf (" py_e_fit %8.6f \n", py_e_fit);
                                             printf (" pz_e_fit %8.6f \n", pz_e_fit);
                                             
                                             printf (" px_p %8.6f \n", px_p);
                                             printf (" py_p %8.6f \n", py_p);
                                             printf (" pz_p %8.6f \n", pz_p);
                                             printf (" Mp %8.6f \n", Mp);
                                             printf ("  sigma_p_proton %8.6f \n", sigma_p_proton);
                                             printf (" sigma_costheta_proton %8.6f \n", sigma_theta_proton);
                                             printf (" sigma_phi_proton %8.6f \n", sigma_phi_proton);
                                             printf (" px_p_fit  %8.6f \n", px_p_fit);
                                             printf (" py_p_fit %8.6f \n", py_p_fit);
                                             printf (" pz_p_fit %8.6f \n", pz_p_fit);
                                             
                                             printf (" px_pip %8.6f \n", px_pip);
                                             printf (" py_pip %8.6f \n", py_pip);
                                             printf (" pz_pip %8.6f \n", pz_pip);
                                             printf (" Mpi %8.3f \n", Mpi);
                                             printf ("  sigma_p_piplus %8.6f \n", sigma_p_piplus);
                                             printf (" sigma_costheta_piplus %8.6f \n", sigma_theta_piplus);
                                             printf (" sigma_phi_piplus %8.6f \n", sigma_phi_piplus);
                                             printf (" px_pip_fit  %8.6f \n", px_pip_fit);
                                             printf (" py_pip_fit %8.6f \n", py_pip_fit);
                                             printf (" pz_pip_fit %8.6f \n", pz_pip_fit);
                                             
                                             printf (" px_pim %8.6f \n", px_pim);
                                             printf (" py_pim %8.6f \n", py_pim);
                                             printf (" pz_pim %8.6f \n", pz_pim);
                                             printf (" Mpi %8.6f \n", Mpi);
                                             printf ("  sigma_p_piminus %8.6f \n", sigma_p_piminus);
                                             printf (" sigma_costheta_piminus %8.6f \n", sigma_theta_piminus);
                                             printf (" sigma_phi_piminus %8.6f \n", sigma_phi_piminus);
                                             printf (" px_pim_fit  %8.6f \n", px_pim_fit);
                                             printf (" py_pim_fit %8.6f \n", py_pim_fit);
                                             printf (" pz_pim_fit %8.6f \n", pz_pim_fit);
                                             */
                                            // Initialize the fit in FCN
                                            gMinuit->Eval(gMinuit->GetNumFreePars (), gradient,*val,parval,1);
                                            gMinuit->Command("MIGRAD");
                                            gMinuit->Command("HESSE");
                                            //  h_func_after_fit->Fill(f);
                                            //        h_func_before_after_fit->Fill(fcn);
                                            // Get fitted 4-vectors
                                            
                                            double px_efit_error,py_efit_error,pz_efit_error;
                                            ii = 0;
                                            px = gMinuit->GetParameter(ii,px_e_fit, px_efit_error);
                                            py = gMinuit->GetParameter(ii+1,py_e_fit,py_efit_error );
                                            pz = gMinuit->GetParameter(ii+2,pz_e_fit,pz_efit_error);
                                            pmag = px*px + py+py + pz*pz; // pmagSq
                                            electron_fit.SetPxPyPzE(px,py,pz,sqrt(pmag+Me*Me));
                                            double px_pfit_error,py_pfit_error,pz_pfit_error;
                                            ii += 3;
                                            px = gMinuit->GetParameter(ii,px_p_fit,px_pfit_error);
                                            py = gMinuit->GetParameter(ii+1,py_p_fit,py_pfit_error);
                                            pz = gMinuit->GetParameter(ii+2,pz_pfit_error,pz_pfit_error);
                                            pmag = px*px + py+py + pz*pz; // pmagSq
                                            proton_fit.SetPxPyPzE(px,py,pz,sqrt(pmag+Mp*Mp));
                                            double px_pipfit_error,py_pipfit_error,pz_pipfit_error;
                                            ii+=3;
                                            px = gMinuit->GetParameter(ii,px_pip_fit,px_pipfit_error);
                                            py = gMinuit->GetParameter(ii+1,py_pip_fit,py_pipfit_error);
                                            pz = gMinuit->GetParameter(ii+2,pz_pip_fit,pz_pipfit_error);
                                            pmag = px*px + py+py + pz*pz; // pmagSq
                                            piplus_fit.SetPxPyPzE(px,py,pz,sqrt(pmag+Mpi*Mpi));
                                            double px_pimfit_error,py_pimfit_error,pz_pimfit_error;
                                            ii+=3;
                                            px = gMinuit->GetParameter(ii+0,px_pim_fit,px_pimfit_error);
                                            py = gMinuit->GetParameter(ii+1,py_pim_fit,py_pimfit_error);
                                            pz = gMinuit->GetParameter(ii+2,pz_pim_fit,pz_pimfit_error);
                                            pmag = px*px + py+py + pz*pz; // pmagSq
                                            piminus_fit.SetPxPyPzE(px,py,pz,sqrt(pmag+Mpi*Mpi));
                                            // Compute invariants from fitted 4-vectors
                                            P4X = beam + target - electron_fit - proton_fit - piplus_fit - piminus_fit;
                                            P4XP = P4X+proton_fit;
                                            P4Xpip = P4X+piplus_fit;
                                            P4Xpim = P4X+piminus_fit;
                                            printf("MX2 = %8.4f \n",P4X.M2());
                                            return;
                                            
                                            // Make revised histograms
                                            h_MX2_fit->Fill(P4X.M2());
                                            h_EX_fit->Fill(P4X.E());
                                            //  recompute perpendicular unit vector yq_fit
                                            eprime3_fit = electron_fit.Vect();
                                            yq_fit=ele3.Cross(eprime3_fit);
                                            yq_fit.Unit();
                                            h_PX_Y_fit->Fill(yq.Dot(P4X.Vect()));
                                            
                                            h_MXP2_fit->Fill(P4XP.M2());
                                            h_MXpip2_fit->Fill(P4Xpip.M2());
                                            h_MXpim2_fit->Fill(P4Xpim.M2());
                                        }
                                    }
                                }
                            }
                        } //  End of Exclusivity section
                        //cos(theta_rest) calculation
                        Miss_pim.SetE(sqrt(Mpi*Mpi + Miss_pim.Vect()*Miss_pim.Vect()));
                        TLorentzVector P4sigma=piplus+piminus;
                        //    hMpipiSq_Exc->Fill(P4sigma.M2());
                        TVector3 beta_rest=P4sigma.BoostVector();
                        Zsigma = beta_rest;
                        beta_rest *= -1.0;
                        TLorentzVector PpiplusRest(piplus);
                        PpiplusRest.Boost(beta_rest);
                        P3piplusRest=PpiplusRest.Vect();
                        Zsigma.Unit();
                        CsHS=Zsigma.Dot(P3piplusRest)/P3piplusRest.Mag();
                        h_CsHs_pi->Fill(piplus.P(),CsHS);
                        //phiH calculation
                        ele3=beam.Vect();
                        eprime3=electrontrigger.Vect();
                        yq=ele3.Cross(eprime3);
                        yq.Unit();
                        zq=-q4vec.Vect();
                        zq.Unit();
                        Ysigma=zq.Cross(Zsigma);
                        Ysigma.Unit();
                        Xsigma=Ysigma.Cross(Zsigma);
                        Xsigma.Unit();
                        phiH=TMath::ATan2(Ysigma.Dot(P3piplusRest), Xsigma.Dot(P3piplusRest))*rad2deg;
                        
                        //phi Calculation
                        double    csphi=(Zsigma.Cross(zq)).Dot(yq);
                        double    sinphi=yq.Dot(Zsigma);
                        phi=TMath::ATan2(sinphi,csphi)*rad2deg;
                        
                        if(CsHS < -0.8)  h_proton_piplus_no_cuts_1->Fill(proton_piplus.M());
                        if(CsHS > 0.0)  h_proton_piplus_no_cuts_2->Fill(proton_piplus.M());
                        
                        if(CsHS> 0.8) h_proton_piminus_no_cuts_1->Fill(proton_piminus.M());
                        if(CsHS < 0.0)h_proton_piminus_no_cuts_2->Fill(proton_piminus.M());
                        if(QSq > 1.0 && WSq > 4.0){
                            
                            if(electrontrigger.Theta()*180./TMath::Pi() < 35 && electrontrigger.Theta()*180./TMath::Pi()> 5){
                                h_Q2_xBj->Fill(xBj,QSq);
                                
                            }
                        }
                        
                        
                        //Computing flux factor
                        double epsilon=1./1.+2.*(q4vec.M2()+(beam.E()-electrontrigger.E())*(beam.E()-electrontrigger.E()))/(4*(beam.E())*(electrontrigger.E())-(q4vec.M2()));
                        double tau_v=Coupling_const/(8.* TMath::Pi()) * (q4vec.M2())/(Mp2*(beam.E())*(beam.E())) * (1-xBj)/(xBj*xBj*xBj)* 1./(1-epsilon);
                        //       printf("flux factor %f \n",tau_v);
                        if(QSq > 1.0 && WSq > 4.0){
                            ht->Fill(-Delta.M2());
                            hW2->Fill(W4vec.M());
                            hQ2_cut->Fill(QSq);
                            h_XBj->Fill(xBj);
                            
                            h_phiHS_1->Fill(phiH);
                            h_CsHS_1->Fill(CsHS);
                            h_pion_invariant_1->Fill(pion_invariant.M());
                            h_phi->Fill(phi);
                            h_y_vector->Fill(xBj,y_vector);
                            
                            
                            if(Miss_total.E() > -0.1 && Miss_total.E() < 0.1){
                                if(Miss_total.M2() > -0.2 && Miss_total.M2() < 0.2){
                                    if(yq.Dot(Miss_total.Vect())> -1 && yq.Dot(Miss_total.Vect())< 1){
                                        if(Rec_Track_DC_Chi2N < Track_DC_Chi2Max){
                                            
                                            h_phiHS->Fill(phiH);
                                            if(pion_invariant.M() > 0.67 && pion_invariant.M() < 0.87){
                                                if(-Delta.M2() > 0.0 && -Delta.M2() < 0.004){
                                                    h_CsHS_t1->Fill(CsHS);
                                                    h_phiHS_t1->Fill(phiH);
                                                }
                                                if(-Delta.M2() > 0.04 && -Delta.M2() < 0.05){
                                                    h_CsHS_t2->Fill(CsHS);
                                                    h_phiHS_t2->Fill(phiH);
                                                }
                                                if(-Delta.M2() > 0.15 && -Delta.M2() < 0.2){
                                                    h_CsHS_t3->Fill(CsHS);
                                                    h_phiHS_t3->Fill(phiH);
                                                }
                                            }
                                            
                                            if(Miss_pip.M2() > -0.5 && Miss_pip.M2() < 0.5){
                                                if(Miss_pim.M2() > -0.5 && Miss_pim.M2() < 0.5){
                                                    if((Miss_proton.M2()-0.9383*0.9383) > -0.5 && (Miss_proton.M2()-0.9383*0.9383) < 0.5){
                                                        h_CsHS_M12->Fill(pion_invariant.M(),CsHS);
                                                        h_CsHS->Fill(CsHS);
                                                        h_Dalitz_1->Fill(pion_invariant.M(),proton_piplus.M());
                                                        h_Dalitz_2->Fill(pion_invariant.M(),proton_piminus.M());
                                                        h_CsHS_Mppip->Fill(proton_piplus.M(), CsHS);
                                                        h_CsHS_Mppim->Fill(proton_piminus.M(),CsHS);
                                                    }
                                                }
                                            }
                                            h_CsHs_phiHs->Fill(CsHS,phiH);
                                            h_phi_M12->Fill(pion_invariant.M(),phi);
                                            h_phiH_M12->Fill(pion_invariant.M(),phiH);
                                            
                                            P122=piplus+piminus;
                                            if((P122.M())> 2*Mpi) {
                                                double E_pipiCM= (W4vec.M2()-Mp2+P122.M2())/(2.*(W4vec.M()));
                                                double P_pipiCM=sqrt(E_pipiCM*E_pipiCM- P122.M2());
                                                double EqCM = (W4vec.M2()-Mp2-QSq)/(2.*(W4vec.M()));
                                                double q3CM = sqrt(EqCM*EqCM+QSq);
                                                double    t0 = q4vec.M2()+P122.M2() - 2.*E_pipiCM*EqCM + 2.*P_pipiCM*q3CM;
                                                double    t00=q4vec.M2() + P122.M2() -2.*q4vec.Dot(P122);
                                                
                                                ht->Fill(t0-t00);
                                            }
                                        }
                                        //      }
                                    }
                                }
                            }
                        }
                        
                        
                        //Collinearity Angle
                        proton3_miss= M122_calctrigger.Vect();
                        proton3_measured=proton.Vect();
                        double collinearity_protoncsthe= proton3_miss.Dot(proton3_measured)/(proton3_miss.Mag()*proton3_measured.Mag());
                        double collinearity_prototheta=TMath:: ACos(collinearity_protoncsthe);
                        h_collinearity->Fill(collinearity_prototheta*180./TMath::Pi());
                        
                        
                        /*     phiH=TMath::ATan2(piplus3.Dot(Yhadron), piplus3.Dot(Xhadron));
                         h_phiHS->Fill(phiH);
                         csphi=(Zhadron.Cross(zq)).Dot(yq);
                         sinphi= yq.Dot(Zhadron);
                         phii=TMath::ATan2(sinphi,csphi)*180./TMath::Pi();
                         h_phi->Fill(phii);*/
                        TVector3 Pmiss=Miss_total.Vect();
                        h_Pmiss_ey->Fill(Pmiss.Dot(yq));
                        h_Pz->Fill(pz);
                        double a_const=px/sin(theta);
                        double b_const=pz-px/tan(theta);
                        h__a_b->Fill(a_const,b_const);
                        
                        
                        i_piminus = i_piminus/multiplex;
                    } //  end of kPart (piminus)
                    i_piplus = i_piplus/multiplex;
                } // end of jPart (piplus)
                i_protons = i_protons/multiplex;
            }  // end of iPart (proton)
            
            
            //  counter++;
            
            
            
        } // While?
    } //File loop
    int nBin=100;
    for (int jBin = 0; jBin<nBin+2; jBin++){
        sum1 += h_Rec_Track_DC->GetBinContent(jBin);
        h_Rec_Track_DC_integral->SetBinContent(jBin,sum1);
    }
    //interal of the chi2/NDF histogram
    for (int iBin = 0; iBin<nBin+2; iBin++){
        sum2 += h_Rec_Track_DC_Chi2->GetBinContent(iBin);
        h_Rec_Track_DC_Chi2N_integral->SetBinContent(iBin,sum2);
    }
    
    fclose(outfile);
    gBenchmark->Stop("timer");
    gBenchmark->Print("timer");
    
    /*  TF1 *fMX2 = new TF1("fMX2",GaussCosh,-5.0,5.0,6);
     fMX2->SetParameter(0,3200);
     fMX2->FixParameter(1,0.0);
     fMX2->SetParameter(2,0.5);
     fMX2->SetParameter(3,3500.);
     fMX2->SetParameter(4,2.9);
     fMX2->SetParameter(5,1.0);*/
    
    TF2 *f2 = new TF2("f2","xygaus",2,8,0.2,0.3);
    TCanvas* can111= new TCanvas();
    can111->cd(1);
    ggIM_dist->DrawCopy();
    
    
    TCanvas* can1 = new TCanvas();
    can1->Divide(2,2);
    can1->cd(1);
    h_proton_piplus_IM_all_cuts_sigmaR->DrawCopy();
    can1->cd(2);
    h_proton_piminus_IM_all_cuts_sigmaR->DrawCopy();
    can1->cd(3);
    h_proton_piplus_IM_all_cuts_rhoR->DrawCopy();
    can1->cd(4);
    h_proton_piminus_IM_all_cuts_rhoR->DrawCopy();
    
    TCanvas* can=new TCanvas();
    can->Divide(3,1);
    can->cd(1);
    hmiss->DrawCopy();
    can->cd(2);
    ht->DrawCopy();
    can->cd(3);
    //  h_tdiff->DrawCopy();
    h_MpipiSq_fcut->DrawCopy();
    
    TCanvas* cEl = new TCanvas();
    cEl->Divide(2,4);
    cEl->cd(1);
    hW2->DrawCopy();
    cEl->cd(2);
    hQ2_cut->DrawCopy();
    double Q2x1=1.5;
    for(int iline=0 ; iline < 4; iline++ ){
        TLine * line = new TLine(Q2x1+0.325*iline,0.0,Q2x1+0.325*iline,20000);
        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        line->Draw();
    }
    for(int kline=0 ; kline < 3; kline++ ){
        double Q2x2=2.8;
        TLine * line1 = new TLine(Q2x2+0.76*kline,0.0,Q2x2+0.76*kline,20000);
        line1->SetLineColor(kRed);
        line1->SetLineWidth(2);
        line1->Draw();
    }
    cEl->cd(3);
    ht->DrawCopy();
    for(int bline=0 ; bline < 6; bline++ ){
        double xt1=0.10;
        TLine * line3 = new TLine(xt1+0.30*bline,0.0,xt1+0.30*bline,9000);
        line3->SetLineColor(kRed);
        line3->SetLineWidth(2);
        line3->Draw();
    }
    for(int cline=0 ; cline < 3; cline++ ){
        double xt2=1.9;
        TLine * line4 = new TLine(xt2+0.80*cline,0.0,xt2+0.80*cline,9000);
        line4->SetLineColor(kRed);
        line4->SetLineWidth(2);
        line4->Draw();
    }
    cEl->cd(4);
    h_XBj->DrawCopy();
    for(int aline=0 ; aline < 6; aline++ ){
        double xBj=0.15;
        TLine * line5 = new TLine(xBj+0.06*aline,0.0,xBj+0.06*aline,20000);
        line5->SetLineColor(kRed);
        line5->SetLineWidth(2);
        line5->Draw();
    }
    cEl->cd(5);
    h_CsHS_1->DrawCopy();
    for(int eline=0 ; eline < 8; eline++ ){
        double xCsHS=-1.0;
        TLine * line7 = new TLine(xCsHS+0.28*eline,0.0,xCsHS+0.28*eline,14000);
        line7->SetLineColor(kRed);
        line7->SetLineWidth(2);
        line7->Draw();
    }
    cEl->cd(6);
    h_phiHS_1->DrawCopy();
    for(int gline=0 ; gline < 8; gline++ ){
        double xphiH=-180.;
        TLine * line9 = new TLine(xphiH+51.42*gline,0.0,xphiH+51.42*gline,3800);
        line9->SetLineColor(kRed);
        line9->SetLineWidth(2);
        line9->Draw();
    }
    cEl->cd(7);
    h_phi->DrawCopy();
    for(int dline=0 ; dline < 8; dline++ ){
        double xphi=-180;
        TLine * line6 = new TLine(xphi+51.42*dline,0.0,xphi+51.42*dline,20000);
        line6->SetLineColor(kRed);
        line6->SetLineWidth(2);
        line6->Draw();
    }
    cEl->cd(8);
    h_pion_invariant_1->DrawCopy();
    for(int fline=0 ; fline < 45; fline++ ){
        double xpion=0.26;
        TLine * line8 = new TLine(xpion+0.04*fline,0.0,xpion+0.04*fline,15000);
        line8->SetLineColor(kRed);
        line8->SetLineWidth(2);
        line8->Draw();
    }
    
    TCanvas* CQ2_xBj = new TCanvas();
    CQ2_xBj->cd(1);
    //   htext1->SetMarkerSize(0.5);
    // htext2->SetMarkerSize(0.5);
    // htext2->SetMarkerColor(kYellow);
    h_Q2_xBj->DrawCopy("colz");
    //   h_Q2_xBj_cut->DrawCopy("same");
    
    TLine * line9 = new TLine(0.05,2.0,0.05,9.0);
    line9->SetLineColor(kRed);
    line9->SetLineWidth(2);
    line9->Draw();
    
    TLine * line10 = new TLine(0.16,2.0,0.16,9.0);
    line10->SetLineColor(kRed);
    line10->SetLineWidth(2);
    line10->Draw();
    
    TLine * line11 = new TLine(0.22,2.0,0.22,9.0);
    line11->SetLineColor(kRed);
    line11->SetLineWidth(2);
    line11->Draw();
    
    TLine * line12 = new TLine(0.28,2.0,0.28,9.0);
    line12->SetLineColor(kRed);
    line12->SetLineWidth(2);
    line12->Draw();
    
    TLine * line13 = new TLine(0.36,2.0,0.36,9.0);
    line13->SetLineColor(kRed);
    line13->SetLineWidth(2);
    line13->Draw();
    
    TLine * line14 = new TLine(0.46,2.0,0.46,9.0);
    line14->SetLineColor(kRed);
    line14->SetLineWidth(2);
    line14->Draw();
    
    TLine * line21 = new TLine(0.59,2.0,0.59,9.0);
    line21->SetLineColor(kRed);
    line21->SetLineWidth(2);
    line21->Draw();
    
    TLine * line15 = new TLine(0.7,2.0,0.7,9.0);
    line15->SetLineColor(kRed);
    line15->SetLineWidth(2);
    line15->Draw();
    
    
    TLine * line16 = new TLine(0.05,2.0,0.7,2.0);
    line16->SetLineColor(kRed);
    line16->SetLineWidth(2);
    line16->Draw();
    
    TLine * line17 = new TLine(0.05,3.0,0.7,3.0);
    line17->SetLineColor(kRed);
    line17->SetLineWidth(2);
    line17->Draw();
    
    
    TLine * line18 = new TLine(0.05,3.6,0.7,3.6);
    line18->SetLineColor(kRed);
    line18->SetLineWidth(2);
    line18->Draw();
    
    TLine * line20 = new TLine(0.05,6.0,0.7,6.0);
    line20->SetLineColor(kRed);
    line20->SetLineWidth(2);
    line20->Draw();
    
    TLine * line19 = new TLine(0.05,9.0,0.7,9.0);
    line19->SetLineColor(kRed);
    line19->SetLineWidth(2);
    line19->Draw();
    
    
    
    
    
    
    
    TCanvas* Cy_xBj = new TCanvas();
    Cy_xBj->cd(1);
    h_y_vector->DrawCopy("colz");
    
    
    
    TCanvas* Cmpipi_phi = new TCanvas();
    Cmpipi_phi->cd(1);
    h_phi_M12->DrawCopy("colz");
    
    TCanvas* Cmpipi_phiH = new TCanvas();
    Cmpipi_phiH->cd(1);
    h_phiH_M12->DrawCopy("colz");
    
    
    TCanvas* eleEnergy = new TCanvas();
    eleEnergy->cd(1);
    h_elecEnergy->DrawCopy();
    
    
    TCanvas* polar = new TCanvas();
    polar->cd(1);
    h_elecPolar->DrawCopy();
    
    
    TCanvas* Tof=new TCanvas();
    Tof->Divide(4,1);
    Tof->cd(1);
    h_betamom->DrawCopy("colz");
    Tof->cd(2);
    h_betamomkaon->DrawCopy("colz");
    Tof->cd(3);
    h_betahadron->DrawCopy("colz");
    Tof->cd(4);
    h_betadiff->DrawCopy("colz");
    
    TCanvas* ECAL=new TCanvas();
    ECAL->Divide(3,1);
    ECAL->cd(1);
    h_Sampling_fraction->DrawCopy();
    ECAL->cd(2);
    h_ECAL->DrawCopy("colz");
    ECAL->cd(3);
    h_projectileY->DrawCopy();
    f2->SetLineColor(2);
    f2->DrawCopy("same");
    
    TCanvas* vertexZ=new TCanvas();
    vertexZ->cd(1);
    h_vertexZ->DrawCopy("colz");
    
    TCanvas* Tofcalc=new TCanvas();
    Tofcalc->Divide(2,2);
    Tofcalc->cd(1);
    h_betacalcp->DrawCopy("colz");
    TF1 *f_proton = new TF1("f_proton","(1+x/(sqrt([mass]*[mass]+x*x)))/2", 2.6, 4.0);
    TF1 *f_x = new TF1("f_x","1.05", 2.6, 4.0);
    //TF1 *f_proton = new TF1("f_proton","(1+x/(sqrt([mass]*[mass]+x*x)))/2", 2.6, 4.0);
    f_proton->SetParameter(0, 0.9383);
    f_proton->SetLineWidth(2);
    f_x->SetLineWidth(2);
    f_x->SetLineColor(2);
    f_proton->SetLineColor(2);
    f_proton->Draw("same");
    f_x->Draw("same");
    Tofcalc->cd(2);
    h_betaEBp->DrawCopy("colz");
    Tofcalc->cd(3);
    h_beta_myParticle->DrawCopy("colz");
    
    gStyle->SetOptStat(1111);
    gPad->Update();
    TPaveStats *statsbox = (TPaveStats*)gPad->GetPrimitive("stats");
    double y1 = statsbox->GetY1NDC(); // (lower) y start position of stats box
    double y2 = statsbox->GetY2NDC(); // (upper) y start position of stats box
    double newy1 = 2 * y1 - y2;   // new (lower) y start position of stats box
    double newy2 = y1;            // new (upper) y start position of stats box
    statsbox->SetY1NDC(newy1);    //set new y start position
    statsbox->SetY2NDC(newy2);    //set new y end position
    
    TCanvas* M12_calc=new TCanvas();
    M12_calc->Divide(2,1);
    M12_calc->cd(1);
    h_Miss_proton_fcut->DrawCopy();
    M12_calc->cd(2);
    h_Miss_proton_all_cuts->SetLineColor(2);
    h_Miss_proton_all_cuts->DrawCopy();
    
    TAxis *axis10 =  h_Miss_proton_all_cuts->GetXaxis();
    int bmin10 = axis10->FindBin(-1.0); //in  case xmin=-1.0
    int bmax10 = axis10->FindBin(1.15); //in  case xmax=1.15
    double integral10=  h_Miss_proton_all_cuts->Integral(bmin10,bmax10);
    cout<<"Integral from -1.0 to 1.15 for proton peak "<<integral10 <<endl;
    
    TAxis *axis11 =  h_Miss_proton_all_cuts->GetXaxis();
    int bmin11 = axis11->FindBin(-0.5); //in  case xmin=-0.5
    int bmax11 = axis11->FindBin(1.15); //in  case xmax=1.15
    double integral11=  h_Miss_proton_all_cuts->Integral(bmin11,bmax11);
    cout<<"Integral from -0.5 to 1.15 for proton peak "<<integral11 <<endl;
    
    TAxis *axis12 =  h_Miss_proton_all_cuts->GetXaxis();
    int bmin12 = axis12->FindBin(0.0); //in  case xmin=0.0
    int bmax12 = axis12->FindBin(1.15); //in  case xmax=1.15
    double integral12=  h_Miss_proton_all_cuts->Integral(bmin12,bmax12);
    cout<<"Integral from 0.0 to 1.15 for proton peak "<<integral12 <<endl;
    
    TAxis *axis13 =  h_Miss_proton_all_cuts->GetXaxis();
    int bmin13 = axis13->FindBin(0.5); //in  case xmin=0.5
    int bmax13 = axis13->FindBin(1.15); //in  case xmax=1.15
    double integral13=  h_Miss_proton_all_cuts->Integral(bmin13,bmax13);
    cout<<"Integral from 0.5 to 1.15 for proton peak "<<integral13 <<endl;
    
    TAxis *axis14 =  h_Miss_proton_all_cuts->GetXaxis();
    int bmin14 = axis14->FindBin(0.7); //in  case xmin=0.7
    int bmax14 = axis14->FindBin(1.15); //in  case xmax=1.15
    double integral14=  h_Miss_proton_all_cuts->Integral(bmin14,bmax14);
    cout<<"Integral from 0.7 to 1.15 for proton peak "<<integral14 <<endl;
    
    
    TCanvas* Cmisspim=new TCanvas();
    Cmisspim->Divide(2,1);
    Cmisspim->cd(1);
    h_Miss_pim_fcut->DrawCopy();
    Cmisspim->cd(2);
    h_Miss_pim_all_cuts->SetLineColor(2);
    h_Miss_pim_all_cuts->DrawCopy();
    
    
    TCanvas* Cmisspip=new TCanvas();
    Cmisspip->Divide(2,1);
    Cmisspip->cd(1);
    h_Miss_pip_fcut->DrawCopy();
    Cmisspip->cd(2);
    h_Miss_pip_all_cuts->SetLineColor(2);
    h_Miss_pip_all_cuts->DrawCopy();
    
    TCanvas* Cprotonpiplus=new TCanvas();
    Cprotonpiplus->Divide(2,1);
    Cprotonpiplus->cd(1);
    h_proton_piplus_no_cuts_1->DrawCopy();
    Cprotonpiplus->cd(2);
    h_proton_piplus_no_cuts_2->DrawCopy();
    
    TCanvas* Cprotonpiminus=new TCanvas();
    Cprotonpiminus->Divide(2,1);
    Cprotonpiminus->cd(1);
    h_proton_piminus_no_cuts_1->DrawCopy();
    Cprotonpiminus->cd(2);
    h_proton_piminus_no_cuts_2->DrawCopy();
    
    
    
    TAxis *axis1 =  h_M122_calctrigger_2->GetXaxis();
    int bmin1 = axis1->FindBin(0.0); //in  case xmin=0.0
    int bmax1 = axis1->FindBin(2.0); //in  case xmax=2.0
    double integral1= h_M122_calctrigger_2->Integral(bmin1,bmax1);
    
    integral1 -= h_M122_calctrigger_2->GetBinContent(bmin1)*(0.0-axis1->GetBinLowEdge(bmin1))/axis1->GetBinWidth(bmin1);
    integral1 -=  h_M122_calctrigger_2->GetBinContent(bmax1)*(axis1->GetBinUpEdge(bmax1)-2.0)/axis1->GetBinWidth(bmax1);
    cout<<"Integral from 0 to 2 for electron trigger histogram(blue) "<<integral1 <<endl;
    /*
     h_M122_calchighE->SetLineColor(2);
     h_M122_calchighE->DrawCopy("sames");
     TAxis *axis2 = h_M122_calchighE->GetXaxis();
     int bmin2 = axis2->FindBin(0.0); //in  case xmin=0.0
     int bmax2 = axis2->FindBin(2.0); //in  case xmax=2.0
     double integral2=  h_M122_calchighE->Integral(bmin2,bmax2);
     integral2 -= h_M122_calchighE->GetBinContent(bmin2)*(0.0-axis2->GetBinLowEdge(bmin2))/axis2->GetBinWidth(bmin2);
     integral2 -=  h_M122_calchighE->GetBinContent(bmax2)*(axis2->GetBinUpEdge(bmax2)-2.0)/axis2->GetBinWidth(bmax2);
     
     cout<<"Integral from 0 to 2 for  high energy electron(red) "<<integral2 <<endl;
     h_M122_calctriggerwithBremsphoton->SetLineColor(8);
     h_M122_calctriggerwithBremsphoton->DrawCopy("sames");*/
    
    
    
    TCanvas* m2k=new TCanvas();
    m2k->cd(1);
    h_m2k->DrawCopy();
    
    TCanvas* Cpid=new TCanvas();
    Cpid->cd(1);
    h_pid->DrawCopy();
    
    TCanvas* Cpimult=new TCanvas();
    Cpimult->cd(1);
    h_pimult->DrawCopy();
    
    TCanvas* Cpimult2D=new TCanvas();
    Cpimult2D->cd(1);
    h_pippimmult->DrawCopy("colz");
    
    
    TCanvas* Ckdotg=new TCanvas();
    Ckdotg->Divide(2,2);
    Ckdotg->cd(1);
    h_kdotg->DrawCopy();
    Ckdotg->cd(2);
    h_kpdotg->DrawCopy();
    Ckdotg->cd(3);
    h_kdotg2D->Draw("colz");
    Ckdotg->cd(4);
    h_photonenergy->Draw();
    
    TCanvas* Cphotonenergy=new TCanvas();
    Cphotonenergy->cd(1);
    h_photonenergy->DrawCopy();
    
    
    
    
    TCanvas* CEmissmisspip=new TCanvas();
    CEmissmisspip->Divide(2,1);
    CEmissmisspip->cd(1);
    h_Emiss_Pmiss_pip->DrawCopy("colz");
    CEmissmisspip->cd(2);
    h_Emiss_Pmiss_pip_cut->DrawCopy("colz");
    
    TCanvas* Cmisstotal=new TCanvas();
    Cmisstotal->Divide(2,1);
    Cmisstotal->cd(1);
    h_Miss_total_fcut->DrawCopy();
    //  h_Miss_total_cut->SetLineColor(2);
    //  h_Miss_total_cut->DrawCopy("sames");
    Cmisstotal->cd(2);
    h_Miss_total_Energy_fcut->DrawCopy();
    //  h_Miss_total_Energy_cut->SetLineColor(2);
    //  h_Miss_total_Energy_cut->DrawCopy("same");
    
    
    
    TCanvas* Celectronmult=new TCanvas();
    Celectronmult->Divide(3,1);
    Celectronmult->cd(1);
    h_electronmult->DrawCopy();
    Celectronmult->cd(2);
    h_electronmultFD->DrawCopy();
    Celectronmult->cd(3);
    h_electronmultiplicity->DrawCopy();
    
    
    TCanvas* Cprotonmult=new TCanvas();
    Cprotonmult->cd(1);
    h_protonmult->DrawCopy();
    
    TCanvas* Cst_vt=new TCanvas();
    Cst_vt->cd(1);
    h_st_vt->DrawCopy();
    
    
    TCanvas* ChiSq= new TCanvas();
    ChiSq->cd(1);
    h_tof_ChiSq->DrawCopy();
    
    
    TCanvas* PCAL_DOCA= new TCanvas();
    PCAL_DOCA->cd(1);
    h_PCAL_DOCA->DrawCopy();
    
    TCanvas* REC_Track= new TCanvas();
    REC_Track->Divide(3,1);
    REC_Track->cd(1);
    h_Rec_Track_DC->DrawCopy();
    REC_Track->cd(2);
    h_Rec_Track_DC_Chi2->DrawCopy();
    REC_Track->cd(3);
    h_Rec_Track_DC_Chi2N->DrawCopy();
    
    
    
    
    TCanvas* CEgamma= new TCanvas();
    CEgamma->cd(1);
    h_Egamma->DrawCopy();
    //    cout << "chi-squre = " << fval << endl;
    
    
    TCanvas* CPCALUVW= new TCanvas();
    CPCALUVW->Divide(3,1);
    CPCALUVW->cd(1);
    h_PCAL_Lu->DrawCopy();
    CPCALUVW->cd(2);
    h_PCAL_Lv->DrawCopy();
    CPCALUVW->cd(3);
    h_PCAL_Lw->DrawCopy();
    
    TCanvas* DC_XY= new TCanvas();
    DC_XY->Divide(2,4);
    DC_XY->cd(1);
    h_PCAL_XY->DrawCopy("colz");
    DC_XY->cd(2);
    h_PCAL_XY_cut->DrawCopy("colz");
    DC_XY->cd(3);
    h_DC_XY1->DrawCopy("colz");
    DC_XY->cd(4);
    h_DC_XY1_cut->DrawCopy("colz");
    DC_XY->cd(5);
    h_DC_XY2->DrawCopy("colz");
    DC_XY->cd(6);
    h_DC_XY2_cut->DrawCopy("colz");
    DC_XY->cd(7);
    h_DC_XY3->DrawCopy("colz");
    DC_XY->cd(8);
    h_DC_XY3_cut->DrawCopy("colz");
    
    TCanvas* CEgamma1= new TCanvas();
    CEgamma1->cd(1);
    h_E_PCAL->DrawCopy("colz");
    
    
    TCanvas* DC_phi_theta= new TCanvas();
    DC_phi_theta->cd(1);
    h_DC_phi_theta->DrawCopy("colz");
    
    TCanvas* Beta_Diff= new TCanvas();
    Beta_Diff->cd(1);
    h_3h_betaDiff->DrawCopy();
    
    TCanvas* Protonmultipli= new TCanvas();
    Protonmultipli->cd(1);
    h_M122_calctrigger_5->DrawCopy();
    h_M122_calctrigger_6->SetLineColor(2);
    h_M122_calctrigger_6->DrawCopy("sames");
    h_M122_calctrigger_6->Clone("hDiff");
    hDiff->Add( h_M122_calctrigger_5, -1);
    hDiff->SetLineColor(7);
    hDiff->DrawCopy("sames");
    
    TCanvas* Protonpolar= new TCanvas();
    Protonpolar->cd(1);
    h_proton_polar->DrawCopy();
    h_pip_polar->SetLineColor(2);
    h_pip_polar->DrawCopy("sames");
    h_pim_polar->SetLineColor(7);
    h_pim_polar->DrawCopy("sames");
    TLegend *legend = new TLegend(0.55,0.65,0.76,0.82);
    legend->AddEntry(h_proton_polar,"Proton Polar angle","l");
    legend->AddEntry(h_pip_polar,"Pip polar angle ","l");
    legend->AddEntry(h_pim_polar,"Pim polar angle ","l");
    legend->Draw();
    
    
    /*  TCanvas* ZvertexDiff= new TCanvas();
     ZvertexDiff->cd(1);
     h_dvz1->DrawCopy();
     h_dvz2->SetLineColor(2);
     h_dvz2->DrawCopy("sames");
     TAxis *axis =  h_dvz1->GetXaxis();
     int bmin = axis->FindBin(-10.0); //in  case xmin=-10.0
     int bmax = axis->FindBin(10.0); //in  case xmax=10.0
     double integral= h_dvz1->Integral(bmin,bmax);
     
     //  integral1 -= h_M122_calctrigger_2->GetBinContent(bmin1)*(0.0-axis1->GetBinLowEdge(bmin1))/axis1->GetBinWidth(bmin1);
     //    integral1 -=  h_M122_calctrigger_2->GetBinContent(bmax1)*(axis1->GetBinUpEdge(bmax1)-2.0)/axis1->GetBinWidth(bmax1);
     cout<<"Integral from -10 to 10 for vertex difference "<<integral <<endl;*/
    
    TCanvas* ECinnerouter= new TCanvas();
    ECinnerouter->cd(1);
    h_ECinner_outer->DrawCopy("colz");
    
    
    TCanvas* Nphe= new TCanvas();
    Nphe->cd(1);
    h_Nphe->DrawCopy();
    
    TCanvas* NpheEtot= new TCanvas();
    NpheEtot->cd(1);
    h_Nphe_Etot->DrawCopy("colz");
    
    
    TCanvas* GENREC= new TCanvas();
    GENREC->Divide(3,1);
    GENREC->cd(1);
    h_GEN_REC_p->DrawCopy("colz");
    GENREC->cd(2);
    h_GEN_REC_theta->DrawCopy("colz");
    GENREC->cd(3);
    h_GEN_REC_phi->DrawCopy("colz");
    
    TCanvas* Q2xBj= new TCanvas();
    auto hBox1  = new TH2F("hbox1","Q2 vs xBj bins",6,0.15,0.55,7,1.5,5.1);
    //  Q2xBj->Divide(5,5);
    Q2xBj->cd(1);
    
    hBox1->DrawCopy();
    Q2xBj->cd(5*1+2);
    hmiss->DrawCopy();
    
    
    TCanvas* eprimeDot= new TCanvas();
    eprimeDot->cd(1);
    // h_Miss_total_Energy->SetLineColor(1);
    h_eprimeDot1->SetLineColor(7);
    // h_eprimeDot1->DrawCopy("sames");
    h_eprimeDot2->SetLineColor(2);
    //  h_eprimeDot2->DrawCopy("sames");
    
    TLegend *leg = new TLegend(0.5, 0.7, 0.4, 0.89);
    leg->AddEntry(h_Miss_total_Energy_fcut->DrawCopy(),"no any selection","l");
    leg->AddEntry(h_eprimeDot1->DrawCopy("sames"), "eprime.Dot(Photon)> 0.05", "l");
    leg->AddEntry(h_eprimeDot2->DrawCopy("sames"), "eprime.Dot(Photon)< 0.05", "l");
    leg->Draw();
    
    
    TCanvas* eprimeDotP= new TCanvas();
    eprimeDotP->cd(1);
    h_eprimeDotP->DrawCopy();
    
    
    TCanvas* ggIM= new TCanvas();
    ggIM->cd(1);
    ggIM_dist->DrawCopy();
    
    TCanvas* Missingmomentum= new TCanvas();
    Missingmomentum->cd(1);
    h_missing_momentum_fcut->DrawCopy();
    
    TCanvas* Pmissey= new TCanvas();
    Pmissey->cd(1);
    h_Pmiss_ey->DrawCopy();
    
    TCanvas* Pz= new TCanvas();
    Pz->cd(1);
    h_Pz->DrawCopy();
    
    TCanvas* ab= new TCanvas();
    ab->cd(1);
    h__a_b->DrawCopy("colz");
    
    
    TCanvas* CQ2= new TCanvas();
    CQ2->cd(1);
    hQ2->DrawCopy();
    hQ2_cut->SetLineColor(2);
    hQ2_cut->DrawCopy("sames");
    
    TCanvas* collinearityProton= new TCanvas();
    collinearityProton->cd(1);
    h_collinearity->DrawCopy();
    
    TCanvas* CChi2pid= new TCanvas();
    CChi2pid->Divide(4,1);
    CChi2pid->cd(1);
    h_chi2pid_electron->DrawCopy();
    CChi2pid->cd(2);
    h_chi2pid_proton->DrawCopy();
    CChi2pid->cd(3);
    h_chi2pid_piplus->DrawCopy();
    CChi2pid->cd(4);
    h_chi2pid_piminus->DrawCopy();
    
    
    TCanvas* C_betamomCD= new TCanvas();
    C_betamomCD->cd(1);
    h_betamomCD->DrawCopy("colz");
    
    TCanvas* C_missingmass_energy= new TCanvas();
    C_missingmass_energy->cd(1);
    h_missingmass_energy->DrawCopy("colz");
    
    TCanvas* C_Pymiss= new TCanvas();
    C_Pymiss->Divide(2,1);
    C_Pymiss->cd(1);
    h_Pymiss_fcut->DrawCopy();
    //  fMX2->SetLineColor(2);
    // fMX2->DrawCopy("sames");
    C_Pymiss->cd(2);
    h_Pymiss_cut->DrawCopy();
    //    fMX2->SetLineColor(2);
    // fMX2->DrawCopy("sames");
    
    TCanvas* C_P_miss= new TCanvas();
    C_P_miss->cd(1);
    h_Miss_total_momentum->DrawCopy();
    
    TCanvas* C_double_pion= new TCanvas();
    C_double_pion->cd(1);
    h_miss_doublepion->DrawCopy();
    
    TCanvas* C_CsHs_phiHs= new TCanvas();
    C_CsHs_phiHs->cd(1);
    h_CsHs_phiHs->DrawCopy("colz");
    
    TCanvas* C_CsHs_M12= new TCanvas();
    C_CsHs_M12->cd(1);
    h_CsHS_M12->DrawCopy("LEGO2Z");
    
    TCanvas* C_CsHs= new TCanvas();
    C_CsHs->cd(1);
    h_CsHS->DrawCopy();
    
    TCanvas* C_Dalitz= new TCanvas();
    C_Dalitz->Divide(2,2);
    C_Dalitz->cd(1);
    h_Dalitz_1->DrawCopy("colz");
    C_Dalitz->cd(2);
    h_Dalitz_1->DrawCopy("lego2z");
    C_Dalitz->cd(3);
    h_Dalitz_2->DrawCopy("colz");
    C_Dalitz->cd(4);
    h_Dalitz_2->DrawCopy("lego2z");
    
    TCanvas* C_CsHs_ppip= new TCanvas();
    C_CsHs_ppip->Divide(2,1);
    C_CsHs_ppip->cd(1);
    h_CsHS_Mppip->DrawCopy("colz");
    C_CsHs_ppip->cd(2);
    h_CsHS_Mppip->DrawCopy("lego2z");
    
    TCanvas* C_CsHs_ppim= new TCanvas();
    C_CsHs_ppim->Divide(2,1);
    C_CsHs_ppim->cd(1);
    h_CsHS_Mppim->DrawCopy("colz");
    C_CsHs_ppim->cd(2);
    h_CsHS_Mppim->DrawCopy("lego2z");
    
    
    TCanvas* C_CsHs_pi= new TCanvas();
    C_CsHs_pi->cd(1);
    h_CsHs_pi->DrawCopy("colz");
    
    TCanvas* CphithetaCD= new TCanvas();
    CphithetaCD->cd(1);
    h_phithetaCD->DrawCopy("colz");
    
    
    TCanvas* C_CsHS_phiH_t= new TCanvas();
    C_CsHS_phiH_t->Divide(2,3);
    C_CsHS_phiH_t->cd(1);
    h_CsHS_t1->DrawCopy();
    C_CsHS_phiH_t->cd(2);
    h_phiHS_t1->DrawCopy();
    C_CsHS_phiH_t->cd(3);
    h_CsHS_t2->DrawCopy();
    C_CsHS_phiH_t->cd(4);
    h_phiHS_t2->DrawCopy();
    C_CsHS_phiH_t->cd(5);
    h_CsHS_t3->DrawCopy();
    C_CsHS_phiH_t->cd(6);
    h_phiHS_t3->DrawCopy();
    
    
    TCanvas* CVertexdiff= new TCanvas();
    CVertexdiff->Divide(3,1);
    CVertexdiff->cd(1);
    h_dvz1->DrawCopy();
    CVertexdiff->cd(2);
    h_dvz2->DrawCopy();
    CVertexdiff->cd(3);
    h_dvz3->DrawCopy();
    
    
    TCanvas* Cbetacalorimeter= new TCanvas();
    Cbetacalorimeter->cd(1);
    h_beta_calorimeter->DrawCopy("colz");
    
    
    TCanvas* Cpn= new TCanvas();
    Cpn->cd(1);
    h_pn->DrawCopy();
    
    TCanvas* CpmagDelta_FTOF= new TCanvas();
    CpmagDelta_FTOF->Divide(2,1);
    CpmagDelta_FTOF->cd(1);
    h_pmag_DeltaT_FTOF_positive->DrawCopy("colz");
    CpmagDelta_FTOF->cd(2);
    h_pmag_DeltaT_FTOF_negative->DrawCopy("colz");
    
    
    TCanvas* CpmagDelta_CTOF= new TCanvas();
    CpmagDelta_CTOF->Divide(2,1);
    CpmagDelta_CTOF->cd(1);
    h_pmag_DeltaT_CTOF_positive->DrawCopy("colz");
    CpmagDelta_CTOF->cd(2);
    h_pmag_DeltaT_CTOF_negative->DrawCopy("colz");
    
    
    TCanvas* CpmagDelta_FTOF_PID= new TCanvas();
    CpmagDelta_FTOF_PID->Divide(3,1);
    CpmagDelta_FTOF_PID->cd(1);
    h_pmag_DeltaT_FTOF_proton->DrawCopy("colz");
    CpmagDelta_FTOF_PID->cd(2);
    h_pmag_DeltaT_FTOF_piplus->DrawCopy("colz");
    CpmagDelta_FTOF_PID->cd(3);
    h_pmag_DeltaT_FTOF_piminus->DrawCopy("colz");
    
    TCanvas* CpmagDelta_CTOF_PID= new TCanvas();
    CpmagDelta_CTOF_PID->Divide(3,1);
    CpmagDelta_CTOF_PID->cd(1);
    h_pmag_DeltaT_CTOF_proton->DrawCopy("colz");
    CpmagDelta_CTOF_PID->cd(2);
    h_pmag_DeltaT_CTOF_piplus->DrawCopy("colz");
    CpmagDelta_CTOF_PID->cd(3);
    h_pmag_DeltaT_CTOF_piminus->DrawCopy("colz");
    
    
    TCanvas* Cmissingtwoparticle= new TCanvas();
    Cmissingtwoparticle->Divide(3,1);
    Cmissingtwoparticle->cd(1);
    h_pion_missing->DrawCopy();
    Cmissingtwoparticle->cd(2);
    h_proton_piplus_missing->DrawCopy();
    Cmissingtwoparticle->cd(3);
    h_proton_piminus_missing->DrawCopy();
    
    
    /*  TCanvas* Kin_fit_1= new TCanvas();
     Kin_fit_1->cd(1);
     h_MX2_fit->DrawCopy();
     
     TCanvas* Kin_fit_2= new TCanvas();
     Kin_fit_2->cd(1);
     h_MX2_fit->DrawCopy();
     
     TCanvas* Kin_fit_3= new TCanvas();
     Kin_fit_3->cd(1);
     h_MX2_fit->DrawCopy();
     
     TCanvas* Kin_fit_4= new TCanvas();
     Kin_fit_4->cd(1);
     h_MX2_fit->DrawCopy();*/
    
    
    
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<counter<< " s\n";
    
    //  t1->Write();
    
    // AnaFile->Close();
    AnaFile->Write();
    
    //  delete AnaFile;
    //
    
}

/*Double_t GaussCosh(Double_t *x, Double_t *par){
 Double_t result = 0.0;
 result  = par[0]*exp(-pow((x[0]-par[1])/par[2],2)/2.);
 result += par[3]/pow(cosh(x[0]/par[4]),par[5]);
 return result;
 }*/

void fcn(int &npar, double *gin, double &f, double *par, int iflag){
    // function to calculate ChisSq for Kinematic fit of ep--> e p pi+ pi-
    
    double chiSq=0;
    double MXSq=0;
    TLorentzVector epTot4vec;
    TLorentzVector P4X;
    TLorentzVector elec4vec_fit, prot4vec_fit, pplus4vec_fit, pminus4vec_fit;
    TLorentzVector elec4vec_0,prot4vec_0,pplus4vec_0,pminus4vec_0;
    TVector3 elec3vec_0, prot3vec_0, pplus3vec_0, pminus3vec_0, yq;
    TVector3 elec3vec_fit,prot3vec_fit,pplus3vec_fit,pminus3vec_fit;
    int ii;
    
    switch (iflag) {
        case 1:
            // Intitialize 4-vectors?
            epTot4vec    = beam+target;
            elec3vec_0   = electrontrigger.Vect();
            prot3vec_0   = proton.Vect();
            pplus3vec_0  = piplus.Vect();
            pminus3vec_0 = piminus.Vect();
            //  keep going and evaluate the initial fit values
            // break;
            
        default:
            ii=0;
            elec3vec_fit.SetXYZ(par[ii+0],par[ii+1],par[ii+2]);
            ii+=3;
            prot3vec_fit.SetXYZ(par[ii+0],par[ii+1],par[ii+2]);
            ii+=3;
            pplus3vec_fit.SetXYZ(par[ii+0],par[ii+1],par[ii+2]);
            ii+=3;
            pminus3vec_fit.SetXYZ(par[ii+0],par[ii+1],par[ii+2]);
            elec4vec_fit.SetVectM(elec3vec_fit,Me);
            prot4vec_fit.SetVectM(prot3vec_fit,Mp);
            pplus4vec_fit.SetVectM(pplus3vec_fit,Mpi);
            pminus4vec_fit.SetVectM(pminus3vec_fit,Mpi);
    }
    // Electron
    chiSq += pow((elec4vec_fit.P()-elec4vec_0.P())/sigma_p_electron,2);
    printf("chiSq = %8.4f \n",chiSq);
    chiSq += pow((elec4vec_fit.Theta()-elec4vec_0.Theta())/sigma_theta_electron,2);
    printf("chiSq = %8.4f\n ",chiSq);
    chiSq += pow((elec4vec_fit.Phi()-elec4vec_0.Phi())/sigma_phi_electron,2);
    printf("chiSq = %8.4f \n",chiSq);
    // Proton
    chiSq += pow((prot4vec_fit.P()-prot4vec_0.P())/sigma_p_proton,2);
    printf("chiSq = %8.4f \n ",chiSq);
    chiSq += pow((prot4vec_fit.Theta()-prot4vec_0.Theta())/sigma_theta_proton,2);
    printf("chiSq = %8.4f \n ",chiSq);
    chiSq += pow((prot4vec_fit.Phi()-prot4vec_0.Phi())/sigma_phi_proton,2);
    printf("chiSq = %8.4f \n ",chiSq);
    // Pi_plus
    chiSq += pow((pplus4vec_fit.P()-pplus4vec_0.P())/sigma_p_piplus,2);
    printf("chiSq = %8.4f \n ",chiSq);
    chiSq += pow((pplus4vec_fit.Theta()-pplus4vec_0.Theta())/sigma_theta_piplus,2);
    printf("chiSq = %8.4f \n ",chiSq);
    chiSq += pow((pplus4vec_fit.Phi()-pplus4vec_0.Phi())/sigma_phi_piplus,2);
    printf("chiSq = %8.4f \n ",chiSq);
    //  Pi_minus
    chiSq += pow((pminus4vec_fit.P()-pminus4vec_0.P())/sigma_p_piminus,2);
    printf("chiSq = %8.4f \n ",chiSq);
    chiSq += pow((pminus4vec_fit.Theta()-pminus4vec_0.Theta())/sigma_theta_piminus,2);
    printf("chiSq = %8.4f \n ",chiSq);
    chiSq += pow((pminus4vec_fit.Phi()-pminus4vec_0.Phi())/sigma_phi_piminus,2);
    printf("chiSq = %8.4f \n ",chiSq);
    //
    P4X = epTot4vec - elec4vec_fit-prot4vec_fit-pplus4vec_fit-pminus4vec_fit;
    chiSq += P4X.M2()*P4X.M2()/(wgt_MX2*wgt_MX2);
    printf("chiSq = %8.4f \n ",chiSq);
    chiSq += P4X.E()*P4X.E()/(wgt_EX*wgt_EX);
    printf("chiSq = %8.4f \n ",chiSq);
    //  Calculate PY
    //  This is the shortcut for fixed target
    yq=beam.Vect();
    yq=yq.Cross(elec3vec_fit);
    yq.Unit();
    double PY = yq.Dot(P4X.Vect());
    chiSq += PY*PY/(wgt_PY*wgt_PY);
    printf("chiSq = %8.4f \n",chiSq);
    f=chiSq;
    
    /*
     for(int j=0; j<nparMINUIT ;j++){
     printf("par[%d] = %8.6f\n", j,par[j]);
     }
     */
}

/// //////////////////////////////////////////////////////////////////////////////////////////////
/// 1. Fiducial cuts for the PCAL (required for electrons and photons)
///

/// 1.1 A homgenous cut for all sectors which does not consider variations of the sampling fraction
///     This cut assumes, that for the electron ID only a proper cluster formation is required, which
///     is given, if the center of the cluster has enough distance to the edges.
///     loose: cuts 1 bar from the outer side (4.5 cm)
///     medium: cuts 2 bars from the outer side (9.0 cm)
///     tight: cuts 3 bars from the outer side (13.5 cm)


bool EC_hit_position_fiducial_cut_homogeneous(region_part_ptr part){
    
    ///////////////////////////
    bool tight = true;
    bool medium = false;
    bool loose = false;
    bool inbending=true;
    bool outbending=false;
    //////////////////////////
    
    // Cut using the natural directions of the scintillator bars/ fibers:
    int PCAL_sector=part->cal(PCAL)->getSector();
    double u = part->cal(PCAL)->getLu();
    double v = part->cal(PCAL)->getLv();
    double w = part->cal(PCAL)->getLw();
    
    /// v + w is going from the side to the back end of the PCAL, u is going from side to side
    /// 1 scintillator bar is 4.5 cm wide. In the outer regions (back) double bars are used.
    
    ///////////////////////////////////////////////////////////////////
    /// inbending:
    //
    double min_u_loose_inb[]  = {4.5,  4.5,  4.5,  4.5,  4.5,  4.5 };
    double min_u_med_inb[]   = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
    double min_u_tight_inb[]   = {13.5, 13.5, 13.5, 13.5, 13.5, 13.5};
    //
    double max_u_loose_inb[] = {398, 398, 398, 398, 398, 398};
    double max_u_med_inb[]   = {408, 408, 408, 408, 408, 408};
    double max_u_tight_inb[] = {420, 420, 420, 420, 420, 420};
    //
    double min_v_loose_inb[] = {4.5,  4.5,  4.5,  4.5,  4.5,  4.5 };
    double min_v_med_inb[]   = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
    double min_v_tight_inb[] = {13.5, 13.5, 13.5, 13.5, 13.5, 13.5};
    //
    double max_v_loose_inb[] = {400, 400, 400, 400, 400, 400};
    double max_v_med_inb[]   = {400, 400, 400, 400, 400, 400};
    double max_v_tight_inb[]   = {400, 400, 400, 400, 400, 400};
    //
    double min_w_loose_inb[] = {4.5,  4.5,  4.5,  4.5,  4.5,  4.5 };
    double min_w_med_inb[]   = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
    double min_w_tight_inb[] = {13.5, 13.5, 13.5, 13.5, 13.5, 13.5};
    //
    double max_w_loose_inb[] = {400, 400, 400, 400, 400, 400};
    double max_w_med_inb[]   = {400, 400, 400, 400, 400, 400};
    double max_w_tight_inb[]  = {400, 400, 400, 400, 400, 400};
    
    ///////////////////////////////////////////////////////////////////////
    /// outbending (not adjusted up to now, same as inbending!):
    //
    double min_u_loose_out[] = {4.5,  4.5,  4.5,  4.5,  4.5,  4.5 };
    double min_u_med_out[]   = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
    double min_u_tight_out[]   = {13.5, 13.5, 13.5, 13.5, 13.5, 13.5};
    //
    double max_u_loose_out[] = {398, 398, 398, 398, 398, 398};
    double max_u_med_out[]   = {408, 408, 408, 408, 408, 408};
    double max_u_tight_out[]  = {420, 420, 420, 420, 420, 420};
    //
    double min_v_loose_out[] = {4.5,  4.5,  4.5,  4.5,  4.5,  4.5 };
    double min_v_med_out[]   = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
    double min_v_tight_out[]  = {13.5, 13.5, 13.5, 13.5, 13.5, 13.5};
    //
    double max_v_loose_out[] = {400, 400, 400, 400, 400, 400};
    double max_v_med_out[]   = {400, 400, 400, 400, 400, 400};
    double max_v_tight_out[]  = {400, 400, 400, 400, 400, 400};
    //
    double min_w_loose_out[] = {4.5,  4.5,  4.5,  4.5,  4.5,  4.5 };
    double min_w_med_out[]   = {9.0,  9.0,  9.0,  9.0,  9.0,  9.0 };
    double min_w_tight_out[]  = {13.5, 13.5, 13.5, 13.5, 13.5, 13.5};
    //
    double max_w_loose_out[] = {400, 400, 400, 400, 400, 400};
    double max_w_med_out[]   = {400, 400, 400, 400, 400, 400};
    double max_w_tight_out[]  = {400, 400, 400, 400, 400, 400};
    
    //////////////////////////////////////////////////////////////
    
    double min_u = 0; double max_u = 0; double min_v = 0; double max_v = 0; double min_w = 0; double max_w = 0;
    
    for(Int_t k = 0; k < 6; k++){
        if(PCAL_sector-1 == k && inbending == true){
            if(tight == true){
                min_u = min_u_tight_inb[k]; max_u = max_u_tight_inb[k];
                min_v = min_v_tight_inb[k]; max_v = max_v_tight_inb[k];
                min_w = min_w_tight_inb[k]; max_w = max_w_tight_inb[k];
            }
            if(medium == true){
                min_u = min_u_med_inb[k]; max_u = max_u_med_inb[k];
                min_v = min_v_med_inb[k]; max_v = max_v_med_inb[k];
                min_w = min_w_med_inb[k]; max_w = max_w_med_inb[k];
            }
            if(loose == true){
                min_u = min_u_loose_inb[k]; max_u = max_u_loose_inb[k];
                min_v = min_v_loose_inb[k]; max_v = max_v_loose_inb[k];
                min_w = min_w_loose_inb[k]; max_w = max_w_loose_inb[k];
            }
        }
        if(PCAL_sector-1 == k && outbending == true){
            if(tight == true){
                min_u = min_u_tight_out[k]; max_u = max_u_tight_out[k];
                min_v = min_v_tight_out[k]; max_v = max_v_tight_out[k];
                min_w = min_w_tight_out[k]; max_w = max_w_tight_out[k];
            }
            if(medium == true){
                min_u = min_u_med_out[k]; max_u = max_u_med_out[k];
                min_v = min_v_med_out[k]; max_v = max_v_med_out[k];
                min_w = min_w_med_out[k]; max_w = max_w_med_out[k];
            }
            if(loose == true){
                min_u = min_u_loose_out[k]; max_u = max_u_loose_out[k];
                min_v = min_v_loose_out[k]; max_v = max_v_loose_out[k];
                min_w = min_w_loose_out[k]; max_w = max_w_loose_out[k];
            }
        }
    }
    
    if(v > min_v && v < max_v && w > min_w && w < max_w) return true;
    else return false;
}

/// 1.2 PCAL fiducial cut based on fitted sampling fraction (I woudl recommend this cut for photons)
///     For this cut, the cut criterium is teh drop of the sampling fraction
///     loose: The mean value of E/p is allowed to drop to 2 RMS of the distribution
///     medium: The mean value of E/p is allowed to drop to 1.5 RMS of the distribution
///     tight: The mean value of E/p is allowed to drop to 1 RMS of the distribution

bool EC_hit_position_fiducial_cut(region_part_ptr part){
    double Pival=TMath::Pi();
    //////////////////////////////////////////////
    bool tight   = false;   // MEAN - 1.0 RMS
    bool medium  = true;    // MEAN - 1.5 RMS
    bool loose   = false;   // MEAN - 2.0 RMS
    bool inbending=true;
    bool outbending=false;
    //////////////////////////////////////////////
    int PCAL_sector=part->cal(PCAL)->getSector();
    double PCAL_x=part->cal(PCAL)->getX();
    double PCAL_y=part->cal(PCAL)->getY();
    double PCAL_z=part->cal(PCAL)->getZ();
    
    double theta_PCAL = 180/Pival * acos(PCAL_z/sqrt(pow(PCAL_x, 2) + pow(PCAL_y, 2) + pow(PCAL_z,2)));
    double phi_PCAL_raw = 180/Pival * atan2(PCAL_y/sqrt(pow(PCAL_x, 2) + pow(PCAL_y, 2) + pow(PCAL_z,2)),
                                            PCAL_x/sqrt(pow(PCAL_x, 2) + pow(PCAL_y, 2) + pow(PCAL_z,2)));
    
    double phi_PCAL = 0;
    if(PCAL_sector == 1) phi_PCAL = phi_PCAL_raw;
    if(PCAL_sector == 2) phi_PCAL = phi_PCAL_raw - 60;
    if(PCAL_sector == 3) phi_PCAL = phi_PCAL_raw - 120;
    if(PCAL_sector == 4 && phi_PCAL_raw > 0) phi_PCAL = phi_PCAL_raw - 180;
    if(PCAL_sector == 4 && phi_PCAL_raw < 0) phi_PCAL = phi_PCAL_raw + 180;
    if(PCAL_sector == 5) phi_PCAL = phi_PCAL_raw + 120;
    if(PCAL_sector == 6) phi_PCAL = phi_PCAL_raw + 60;
    
    // 2 sigma (inb adjusted):
    
    double par_0_min_2sigma[] = {13.771, 25.639, 28.4616, 34.2333, 41.777, 18.041};
    double par_1_min_2sigma[] = {-14.3952, -27.4437, -29.5074, -40.9268, -33.056, -19.1539};
    double par_2_min_2sigma[] = {-0.0579567, 2.12184, 2.7033, 4.72287, 1.80765, 0.648548};
    double par_3_min_2sigma[] = {0.0190146, -0.0315046, -0.0523381, -0.0919215, -0.00542624, 0.00836905};
    double par_4_min_2sigma[] = {-0.000222315, 0.000243574, 0.000453275, 0.000771317, -0.000150837, -0.000208545};
    double par_0_max_2sigma[] = {-19.2009, -16.4848, -47.8295, -24.0029, -25.096, -19.2967};
    double par_1_max_2sigma[] = {19.3148, 11.0556, 47.4188, 23.6525, 19.3032, 19.9627};
    double par_2_max_2sigma[] = {-0.66582, 1.35067, -5.54184, -1.20742, 0.0744728, -0.714339};
    double par_3_max_2sigma[] = {-0.00592537, -0.0623736, 0.123022, 0.00276304, -0.0334298, -0.0077081};
    double par_4_max_2sigma[] = {0.000187643, 0.000759023, -0.00120291, 0.000128345, 0.000486216, 0.000224336};
    
    double par_0_min_2sigma_out[] = {64.895, 69.4726, 37.4087, 57.2897, 36.5138, 66.2482};
    double par_1_min_2sigma_out[] = {-49.7813, -52.7634, -28.7868, -44.398, -27.2841, -50.6646};
    double par_2_min_2sigma_out[] = {3.79889, 4.06934, 1.52759, 3.26992, 1.25399, 3.8817};
    double par_3_min_2sigma_out[] = {-0.0389919, -0.0418169, -0.0120817, -0.0329207, -0.00772183, -0.0398076};
    double par_4_min_2sigma_out[] = {0,0,0,0,0,0};
    double par_0_max_2sigma_out[] = {-58.8252, -42.022, -35.6843, -62.0889, -25.7336, -62.4078};
    double par_1_max_2sigma_out[] = {46.3788, 32.7105, 29.2649, 48.8274, 21.4091, 49.5489};
    double par_2_max_2sigma_out[] = {-3.47241, -2.00905, -1.83462, -3.7453, -0.830711, -3.86809};
    double par_3_max_2sigma_out[] = {0.0350037, 0.018892, 0.0188074, 0.0381712, 0.00422423, 0.040186};
    double par_4_max_2sigma_out[] = {0,0,0,0,0,0};
    
    // 1.5 sigma (inb adjusted, 4 min replaced by 2 min!):
    
    double par_0_min_15sigma[] = {25.2996, 19.3705, 59.5003, 19.3705, 26.9823, 21.8217};
    double par_1_min_15sigma[] = {-26.1158, -19.5271, -55.9639, -19.5271, -22.4489, -23.0262};
    double par_2_min_15sigma[] = {2.09145, 0.72118, 6.4372, 0.72118, 0.890449, 1.48469};
    double par_3_min_15sigma[] = {-0.041483, 0.00293465, -0.13059, 0.00293465, -0.000774943, -0.0190315};
    double par_4_min_15sigma[] = {0.000456261, -0.000109323, 0.00115246, -0.000109323, -7.12074e-05, 0.000125463};
    double par_0_max_15sigma[] = {-26.7425, -15.9009, -47.556, -21.5038, -33.9197, -24.0325};
    double par_1_max_15sigma[] = {26.004, 10.5989, 41.6295, 21.1734, 31.1811, 23.3122};
    double par_2_max_15sigma[] = {-1.76638, 1.10168, -3.66934, -0.969572, -2.51229, -1.23308};
    double par_3_max_15sigma[] = {0.0227414, -0.0455311, 0.0493171, 0.00373945, 0.0443308, 0.00762314};
    double par_4_max_15sigma[] = {-0.000111102, 0.000503536, -0.000151053, 5.22425e-05, -0.000403627, 4.24553e-05};
    
    double par_0_min_15sigma_out[] = {57.0314, 70.411, 74.9683, 53.9146, 44.7614, 64.012};
    double par_1_min_15sigma_out[] = {-43.0803, -53.3163, -59.0842, -41.6436, -35.2193, -48.8726};
    double par_2_min_15sigma_out[] = {2.99452, 4.11397, 5.27234, 2.96164, 2.44681, 3.67536};
    double par_3_min_15sigma_out[] = {-0.0287862, -0.042257, -0.0638554, -0.0293148, -0.0264986, -0.0372214};
    double par_4_min_15sigma_out[] = {0,0,0,0,0,0};
    double par_0_max_15sigma_out[] = {-52.0537, -48.3703, -94.6197, -54.6123, -79.9164, -55.3222};
    double par_1_max_15sigma_out[] = {40.7573, 38.333, 73.5425, 42.9251, 64.9277, 42.8186};
    double par_2_max_15sigma_out[] = {-2.8105, -2.83403, -6.88649, -3.08431, -6.18973, -2.99337};
    double par_3_max_15sigma_out[] = {0.0266371, 0.0318955, 0.0851474, 0.0301575, 0.0787016, 0.0286132};
    double par_4_max_15sigma_out[] = {0,0,0,0,0,0};
    
    // 1 sigma (inb adjusted):
    
    double par_0_min_1sigma[] = {34.1128, 26.6064, 65.3241, 44.0344, 54.5539, 25.7356};
    double par_1_min_1sigma[] = {-30.7179, -26.3373, -58.7761, -35.918, -47.3194, -25.3968};
    double par_2_min_1sigma[] = {2.31272, 1.85141, 6.48495, 2.34733, 4.58872, 1.76128};
    double par_3_min_1sigma[] = {-0.0351384, -0.023687, -0.121042, -0.0170119, -0.0778135, -0.0243688};
    double par_4_min_1sigma[] = {0.000262438, 0.000120765, 0.000936822, -7.66933e-05, 0.000539922, 0.000156061};
    double par_0_max_1sigma[] = {-31.6364, -28.7094, -35.6017, -30.1334, -61.5491, -30.9496};
    double par_1_max_1sigma[] = {27.253, 20.8471, 26.4139, 27.7852, 55.5266, 28.8408};
    double par_2_max_1sigma[] = {-1.53381, -0.254236, -0.77312, -1.85849, -6.03473, -2.07467};
    double par_3_max_1sigma[] = {0.00817262, -0.0259995, -0.0300504, 0.0211806, 0.113112, 0.0285026};
    double par_4_max_1sigma[] = {0.000129402, 0.000495665, 0.000740296, -6.30543e-05, -0.000871381, -0.000162669};
    
    double par_0_min_1sigma_out[] = {58.8694, 75.494, 119.951, 52.9731, 111.593, 53.6272};
    double par_1_min_1sigma_out[] = {-43.6605, -55.9065, -89.9595, -41.2703, -85.093, -40.5409};
    double par_2_min_1sigma_out[] = {3.03906, 4.28761, 8.38977, 3.00994, 8.06472, 2.75155};
    double par_3_min_1sigma_out[] = {-0.0301444, -0.0440019, -0.0999269, -0.0313006, -0.0984858, -0.0264477};
    double par_4_min_1sigma_out[] = {0,0,0,0,0,0};
    double par_0_max_1sigma_out[] = {-40.256, -58.3938, -60.3614, -57.7244, -102.98, -51.0424};
    double par_1_max_1sigma_out[] = {31.4367, 43.3923, 45.8203, 42.9619, 79.613, 38.5613};
    double par_2_max_1sigma_out[] = {-1.78797, -3.14225, -3.49161, -2.9105, -7.51346, -2.46405};
    double par_3_max_1sigma_out[] = {0.0146775, 0.0334355, 0.0366689, 0.0277125, 0.0916122, 0.0223185};
    double par_4_max_1sigma_out[] = {0,0,0,0,0,0};
    
    
    double par0_min = 0;
    double par1_min = 0;
    double par2_min = 0;
    double par3_min = 0;
    double par4_min = 0;
    double par0_max = 0;
    double par1_max = 0;
    double par2_max = 0;
    double par3_max = 0;
    double par4_max = 0;
    
    
    if(tight == true){
        for(int d=0; d<6; d++){
            if(PCAL_sector == d+1){
                par0_min = par_0_min_1sigma[d]; par1_min = par_1_min_1sigma[d]; par2_min = par_2_min_1sigma[d]; par3_min = par_3_min_1sigma[d]; par4_min = par_4_min_1sigma[d];
                par0_max = par_0_max_1sigma[d]; par1_max = par_1_max_1sigma[d]; par2_max = par_2_max_1sigma[d]; par3_max = par_3_max_1sigma[d]; par4_max = par_4_max_1sigma[d];
                if(outbending == true){
                    par0_min = par_0_min_1sigma_out[d]; par1_min = par_1_min_1sigma_out[d]; par2_min = par_2_min_1sigma_out[d]; par3_min = par_3_min_1sigma_out[d]; par4_min = par_4_min_1sigma_out[d];
                    par0_max = par_0_max_1sigma_out[d]; par1_max = par_1_max_1sigma_out[d]; par2_max = par_2_max_1sigma_out[d]; par3_max = par_3_max_1sigma_out[d]; par4_max = par_4_max_1sigma_out[d];
                }
            }
        }
    }
    
    if(medium == true){
        for(int d=0; d<6; d++){
            if(PCAL_sector == d+1){
                par0_min = par_0_min_15sigma[d]; par1_min = par_1_min_15sigma[d]; par2_min = par_2_min_15sigma[d]; par3_min = par_3_min_15sigma[d]; par4_min = par_4_min_15sigma[d];
                par0_max = par_0_max_15sigma[d]; par1_max = par_1_max_15sigma[d]; par2_max = par_2_max_15sigma[d]; par3_max = par_3_max_15sigma[d]; par4_max = par_4_max_15sigma[d];
                if(outbending == true){
                    par0_min = par_0_min_15sigma_out[d]; par1_min = par_1_min_15sigma_out[d]; par2_min = par_2_min_15sigma_out[d]; par3_min = par_3_min_15sigma_out[d]; par4_min = par_4_min_15sigma_out[d];
                    par0_max = par_0_max_15sigma_out[d]; par1_max = par_1_max_15sigma_out[d]; par2_max = par_2_max_15sigma_out[d]; par3_max = par_3_max_15sigma_out[d]; par4_max = par_4_max_15sigma_out[d];
                }
            }
        }
    }
    
    if(loose == true){
        for(int d=0; d<6; d++){
            if(PCAL_sector == d+1){
                par0_min = par_0_min_2sigma[d]; par1_min = par_1_min_2sigma[d]; par2_min = par_2_min_2sigma[d]; par3_min = par_3_min_2sigma[d]; par4_min = par_4_min_2sigma[d];
                par0_max = par_0_max_2sigma[d]; par1_max = par_1_max_2sigma[d]; par2_max = par_2_max_2sigma[d]; par3_max = par_3_max_2sigma[d]; par4_max = par_4_max_2sigma[d];
                if(outbending == true){
                    par0_min = par_0_min_2sigma_out[d]; par1_min = par_1_min_2sigma_out[d]; par2_min = par_2_min_2sigma_out[d]; par3_min = par_3_min_2sigma_out[d]; par4_min = par_4_min_2sigma_out[d];
                    par0_max = par_0_max_2sigma_out[d]; par1_max = par_1_max_2sigma_out[d]; par2_max = par_2_max_2sigma_out[d]; par3_max = par_3_max_2sigma_out[d]; par4_max = par_4_max_2sigma_out[d];
                }
            }
        }
    }
    
    double min_phi = par0_min + par1_min * log(theta_PCAL) + par2_min * theta_PCAL + par3_min *theta_PCAL*theta_PCAL + par4_min *theta_PCAL*theta_PCAL*theta_PCAL;
    double max_phi = par0_max + par1_max * log(theta_PCAL) + par2_max * theta_PCAL + par3_max *theta_PCAL*theta_PCAL + par4_max *theta_PCAL*theta_PCAL*theta_PCAL;
    
    if(phi_PCAL >= min_phi && phi_PCAL <= max_phi) return true;
    else return false;
}



bool EC_sampling_fraction_cut(region_part_ptr part){
    
    /////////////////////////////////////
    bool tight = true;
    bool medium = false;
    bool loose = false;
    bool inbending=true;
    bool outbending=false;
    ////////////////////////////////////
    
    int PCAL_sector=part->cal(PCAL)->getSector();
    double p_mag=part->par()->getP();
    double ECAL_total= (part->cal(PCAL)->getEnergy())+(part->cal(ECIN)->getEnergy())+(part->cal(ECOUT)->getEnergy());
    double p_min =2.0;
    double sigma_range = 4;
    
    if(tight == true){  sigma_range = 3.5;}
    if(medium == true){ sigma_range = 4;}
    if(loose == true){  sigma_range = 5;}
    
    ///
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    /// a) cut based on sampling fraction versus drift chamber momentum
    
    // 10.6 GeV inbending (adjusted to pass 1 data)
    
    double p0mean[] = {0.116405, 0.118783, 0.115558, 0.112439, 0.120496,
        0.117107};
    double p1mean[] = {-0.0155433, 0.109417, -0.0088641, -0.0704795,
        0.100674, -0.0164561};
    double p2mean[] = {0.00722919, 0.0064115, 0.00773618, 0.00769267,
        0.00393716, 0.00702487};
    double p3mean[] = {-0.000859442, -0.000704544, -0.000853684,
        -0.00073044, -0.000403603, -0.000825457};
    double p0sigma[] = {0.0164086, 0.00448347, 0.00770184, 0.00412631,
        0.00687372, 0.00747944};
    double p1sigma[] = {0.00771924, 0.020591, 0.0164028, 0.0217636,
        0.0183004, 0.0184456};
    double p2sigma[] = {-0.00125094, 0.000362296, -0.000381947, 0.000232047,
        -0.000439471, -0.000560982};
    double p3sigma[] = {8.43022e-05, -2.94663e-05, 2.47408e-05,
        -5.70412e-06, 5.14683e-05, 5.80948e-05};
    
    // 10.6 GeV outbending (not adjusted)
    
    double p0mean_out[] = {0.11946, 0.120319, 0.121442, 0.117984, 0.126656,
        0.120343};
    double p1mean_out[] = {-0.00388941, 0.566384, 0.346713, -0.100043,
        0.766969, 0.394986};
    double p2mean_out[] = {0.00693105, 0.00702178, 0.0040979, 0.00381667,
        0.000243406, 0.00563105};
    double p3mean_out[] = {-0.000721964, -0.000688366, -0.000459883,
        -0.000359391, -0.000113102, -0.000620295};
    double p0sigma_out[] = {0.00831965, 0.0151432, 0.0122673, 0.0114592,
        0.014943, 0.00711617};
    double p1sigma_out[] = {0.0156242, 0.00691968, 0.0106717, 0.0116251,
        0.00837599, 0.0158003};
    double p2sigma_out[] = {-0.000538767, -0.000810607, -0.00120787,
        -0.000734793, -0.00160916, -2.0166e-05};
    double p3sigma_out[] = {3.91481e-05, 2.00367e-05, 7.92897e-05,
        2.99526e-05, 0.00011959, -2.08114e-05};
    
    
    if(outbending == true){    // overwrite variables with outbending parameters
        for(Int_t i = 0; i < 6; i++){
            p0mean[i] = p0mean_out[i];
            p1mean[i] = p1mean_out[i];
            p2mean[i] = p2mean_out[i];
            p3mean[i] = p3mean_out[i];
            p0sigma[i] = p0sigma_out[i];
            p1sigma[i] = p1sigma_out[i];
            p2sigma[i] = p2sigma_out[i];
            p3sigma[i] = p3sigma_out[i];
        }
    }
    
    double mean = 0;
    double sigma = 0;
    double upper_lim_total = 0;
    double lower_lim_total = 0;
    
    for(Int_t k = 0; k < 6; k++){
        if(PCAL_sector-1 == k){
            mean = p0mean[k] *( 1 + p_mag/sqrt(pow(p_mag,2) +
                                               p1mean[k])) + p2mean[k] * p_mag + p3mean[k] * pow(p_mag,2);
            sigma = p0sigma[k] + p1sigma[k] / sqrt(p_mag) + p2sigma[k] *
            p_mag + p3sigma[k] * pow(p_mag,2);
            upper_lim_total = mean + sigma_range * sigma;
            lower_lim_total = mean - sigma_range * sigma;
        }
    }
    
    ///
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    if(ECAL_total/p_mag <= upper_lim_total &&
       ECAL_total/p_mag >= lower_lim_total && p_mag > p_min) return true;
    else return false;
}


/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 2.1 DC fiducial cuts for the 3 regions based on the drop of the count rate at the edeges of the DC  (only for electrons!)
///     loose: count rate can drop to 35 %
///     medium: count rate can drop to 50 %
///     tight: count rate can drop to 65 %


bool DC_hit_position_counts_fiducial_cut(region_part_ptr part, int region){
    double Pival=TMath::Pi();
    ///////////////////////////
    bool tight = false;
    bool medium = true;
    bool loose = false;
    bool inbending=true;
    bool outbending=false;
    //////////////////////////
    
    double  part_DC_c1x=part->traj(DC,6)->getX();
    double  part_DC_c1y=part->traj(DC,6)->getY();
    double  part_DC_c1z=part->traj(DC,6)->getZ();
    double  part_DC_c2x=part->traj(DC,18)->getX();
    double  part_DC_c2y=part->traj(DC,18)->getY();
    double  part_DC_c2z=part->traj(DC,18)->getZ();
    double  part_DC_c3x=part->traj(DC,36)->getX();
    double  part_DC_c3y=part->traj(DC,36)->getY();
    double  part_DC_c3z=part->traj(DC,36)->getY();
    
    
    // calculate theta and phi local:
    
    double theta_DC = 0;
    double phi_DC_raw = 0;
    
    if(region == 1){
        theta_DC = 180/Pival * acos(part_DC_c1z/sqrt(pow(part_DC_c1x, 2) + pow(part_DC_c1y, 2) + pow(part_DC_c1z,2)));
        phi_DC_raw = 180/Pival * atan2(part_DC_c1y/sqrt(pow(part_DC_c1x, 2) + pow(part_DC_c1y, 2) + pow(part_DC_c1z,2)),
                                       part_DC_c1x/sqrt(pow(part_DC_c1x, 2) + pow(part_DC_c1y, 2) + pow(part_DC_c1z,2)));
    }
    if(region == 2){
        theta_DC = 180/Pival * acos(part_DC_c2z/sqrt(pow(part_DC_c2x, 2) + pow(part_DC_c2y, 2) + pow(part_DC_c2z,2)));
        phi_DC_raw = 180/Pival * atan2(part_DC_c2y/sqrt(pow(part_DC_c2x, 2) + pow(part_DC_c2y, 2) + pow(part_DC_c2z,2)),
                                       part_DC_c2x/sqrt(pow(part_DC_c2x, 2) + pow(part_DC_c2y, 2) + pow(part_DC_c2z,2)));
    }
    if(region == 3){
        theta_DC = 180/Pival * acos(part_DC_c3z/sqrt(pow(part_DC_c3x, 2) + pow(part_DC_c3y, 2) + pow(part_DC_c3z,2)));
        phi_DC_raw = 180/Pival * atan2(part_DC_c3y/sqrt(pow(part_DC_c3x, 2) + pow(part_DC_c3y, 2) + pow(part_DC_c3z,2)),
                                       part_DC_c3x/sqrt(pow(part_DC_c3x, 2) + pow(part_DC_c3y, 2) + pow(part_DC_c3z,2)));
    }
    
    
    double phi_DC = 0;
    if(part->trk(DC)->getSector() == 1) phi_DC = phi_DC_raw;
    if(part->trk(DC)->getSector() == 2) phi_DC = phi_DC_raw - 60;
    if(part->trk(DC)->getSector() == 3) phi_DC = phi_DC_raw - 120;
    if(part->trk(DC)->getSector() == 4 && phi_DC_raw > 0) phi_DC = phi_DC_raw - 180;
    if(part->trk(DC)->getSector() == 4 && phi_DC_raw < 0) phi_DC = phi_DC_raw + 180;
    if(part->trk(DC)->getSector() == 5) phi_DC = phi_DC_raw + 120;
    if(part->trk(DC)->getSector() == 6) phi_DC = phi_DC_raw + 60;
    
    // inbending cut parameters (adjusted):
    
    // medium (50 %):
    
    double reg1_min_sec1_inb[] = {25.6464, -26.0938, 1.36205, 0.0083525, -0.000433074};
    double reg1_max_sec1_inb[] = {-21.002, 19.3141, 0.183817, -0.061424, 0.00117881};
    double reg2_min_sec1_inb[] = {20.8156, -19.1708, -0.244432, 0.0652386, -0.00120129};
    double reg2_max_sec1_inb[] = {-20.2319, 19.2371, 0.0491417, -0.0527963, 0.000995253};
    double reg3_min_sec1_inb[] = {22.9342, -24.6621, 1.25766, 0.0133507, -0.00049441};
    double reg3_max_sec1_inb[] = {-8.72949, 4.50811, 3.13299, -0.15452, 0.00234593};
    double reg1_min_sec2_inb[] = {16.0752, -15.8399, -0.65495, 0.0740595, -0.00133068};
    double reg1_max_sec2_inb[] = {-40.921, 40.2558, -3.67186, 0.0498104, -0.000151358};
    double reg2_min_sec2_inb[] = {20.5524, -23.7002, 1.20081, 0.0126066, -0.000501927};
    double reg2_max_sec2_inb[] = {-22.7933, 18.397, 0.707315, -0.0848617, 0.00150965};
    double reg3_min_sec2_inb[] = {7.84857, -9.62769, -1.39224, 0.0876156, -0.00137369};
    double reg3_max_sec2_inb[] = {-24.5099, 24.8373, -1.22162, -0.0135951, 0.00051455};
    double reg1_min_sec3_inb[] = {18.2067, -12.322, -1.79675, 0.111661, -0.0017831};
    double reg1_max_sec3_inb[] = {-20.4639, 11.3285, 2.56498, -0.149755, 0.00243423};
    double reg2_min_sec3_inb[] = {17.5093, -13.0897, -1.44801, 0.0995944, -0.00163973};
    double reg2_max_sec3_inb[] = {-29.4241, 25.1723, -0.643032, -0.0433457, 0.0010157};
    double reg3_min_sec3_inb[] = {19.6412, -18.3363, -0.139482, 0.0623644, -0.00122319};
    double reg3_max_sec3_inb[] = {-23.4069, 20.1087, 0.0111006, -0.0595069, 0.00123094};
    double reg1_min_sec4_inb[] = {27.7329, -24.8668, 0.855257, 0.0229504, -0.000588103};
    double reg1_max_sec4_inb[] = {-25.1595, 24.3626, -0.826633, -0.0280833, 0.000694613};
    double reg2_min_sec4_inb[] = {35.3061, -35.9815, 3.34709, -0.0587318, 0.000517383};
    double reg2_max_sec4_inb[] = {-23.5005, 23.5232, -0.872864, -0.0214562, 0.000570724};
    double reg3_min_sec4_inb[] = {13.1517, -8.73139, -2.18634, 0.115133, -0.00170354};
    double reg3_max_sec4_inb[] = {-14.1635, 13.6887, 0.949846, -0.0787818, 0.00128683};
    double reg1_min_sec5_inb[] = {23.5653, -17.285, -1.21053, 0.108096, -0.00189499};
    double reg1_max_sec5_inb[] = {-23.3017, 16.1351, 1.45328, -0.112999, 0.00193088};
    double reg2_min_sec5_inb[] = {30.5631, -30.3887, 2.11072, -0.00879513, -0.000281135};
    double reg2_max_sec5_inb[] = {-18.7762, 11.5513, 2.23634, -0.132, 0.00208225};
    double reg3_min_sec5_inb[] = {14.2375, -10.0129, -2.03418, 0.1236, -0.00199283};
    double reg3_max_sec5_inb[] = {-23.8159, 21.7783, -0.373554, -0.0463973, 0.0010055};
    double reg1_min_sec6_inb[] = {34.9766, -38.2253, 3.80455, -0.0624294, 0.000385789};
    double reg1_max_sec6_inb[] = {-23.6287, 21.0793, 0.0916929, -0.0641014, 0.00125089};
    double reg2_min_sec6_inb[] = {21.6462, -23.9106, 1.1985, 0.0108284, -0.000435214};
    double reg2_max_sec6_inb[] = {-22.1881, 21.3921, -0.325758, -0.0406509, 0.000821165};
    double reg3_min_sec6_inb[] = {12.3393, -12.5734, -1.07645, 0.0836258, -0.00137295};
    double reg3_max_sec6_inb[] = {-20.6041, 23.7403, -1.43195, 0.00507544, 0.000160285};
    
    // loose (35 %):
    double reg1_min_sec1_inb_loose[] = {24.9198, -26.6682, 1.22328, 0.0326184, 0.0326184};
    double reg1_max_sec1_inb_loose[] = {-38.0903, 44.8413, -5.13684, 0.0913599, -0.000566229};
    double reg2_min_sec1_inb_loose[] = {34.3872, -41.6892, 4.78076, -0.0904016, 0.000691157};
    double reg2_max_sec1_inb_loose[] = {-36.582, 45.4679, -5.70527, 0.121105, -0.00107292};
    double reg3_min_sec1_inb_loose[] = {38.5882, -51.7707, 7.65247, -0.198183, 0.00222671};
    double reg3_max_sec1_inb_loose[] = {-35.1098, 46.9821, -6.62713, 0.16512, -0.00177347};
    double reg1_min_sec2_inb_loose[] = {26.3964, -32.5314, 2.87032, -0.0261837, -0.0261837};
    double reg1_max_sec2_inb_loose[] = {-45.4528, 53.4467, -6.83895, 0.144609, -0.00127358};
    double reg2_min_sec2_inb_loose[] = {29.5885, -38.6103, 4.45494, -0.084114, 0.000626459};
    double reg2_max_sec2_inb_loose[] = {-42.6042, 50.783, -6.44381, 0.135305, -0.00117332};
    double reg3_min_sec2_inb_loose[] = {21.3286, -33.2316, 4.16718, -0.0947928, 0.000950229};
    double reg3_max_sec2_inb_loose[] = {-26.2955, 32.832, -3.21406, 0.0443669, -0.000110018};
    double reg1_min_sec3_inb_loose[] = {38.3683, -39.555, 3.48786, -0.0311045, -0.0311045};
    double reg1_max_sec3_inb_loose[] = {-44.7002, 48.5238, -5.33933, 0.0857334, -0.000389759};
    double reg2_min_sec3_inb_loose[] = {26.2411, -26.4956, 1.17454, 0.0314168, -0.000923764};
    double reg2_max_sec3_inb_loose[] = {-49.1042, 53.7936, -6.50797, 0.125435, -0.000940508};
    double reg3_min_sec3_inb_loose[] = {30.2087, -36.2303, 3.94272, -0.070659, 0.000484111};
    double reg3_max_sec3_inb_loose[] = {-25.1075, 27.5492, -1.81589, -0.00372993, 0.000534551};
    double reg1_min_sec4_inb_loose[] = {45.5355, -50.2561, 5.91736, -0.116203, -0.116203};
    double reg1_max_sec4_inb_loose[] = {-45.3956, 52.8208, -6.57583, 0.13245, -0.00106037};
    double reg2_min_sec4_inb_loose[] = {42.9007, -48.8609, 6.02101, -0.13073, 0.00125121};
    double reg2_max_sec4_inb_loose[] = {-32.1523, 38.6933, -4.15879, 0.0722104, -0.000474606};
    double reg3_min_sec4_inb_loose[] = {34.0991, -40.7811, 4.8997, -0.107033, 0.00105086};
    double reg3_max_sec4_inb_loose[] = {-39.189, 52.9641, -8.01756, 0.214872, -0.00250226};
    double reg1_min_sec5_inb_loose[] = {52.0273, -61.3244, 8.45187, -0.193011, -0.193011};
    double reg1_max_sec5_inb_loose[] = {-50.6284, 54.3673, -6.38031, 0.116905, -0.000778886};
    double reg2_min_sec5_inb_loose[] = {32.1495, -35.3084, 3.11177, -0.0314703, -8.39833e-05};
    double reg2_max_sec5_inb_loose[] = {-47.0576, 50.676, -5.80697, 0.10389, -0.000677858};
    double reg3_min_sec5_inb_loose[] = {27.5659, -34.8935, 3.81989, -0.0655272, 0.000389246};
    double reg3_max_sec5_inb_loose[] = {-25.7138, 27.0011, -1.53335, -0.0152295, 0.000694336};
    double reg1_min_sec6_inb_loose[] = {35.6586, -43.3351, 5.05374, -0.0961372, -0.0961372};
    double reg1_max_sec6_inb_loose[] = {-43.6986, 52.1413, -6.61602, 0.137247, -0.00115149};
    double reg2_min_sec6_inb_loose[] = {33.3586, -42.3939, 5.17582, -0.106612, 0.000910642};
    double reg2_max_sec6_inb_loose[] = {-38.3305, 46.3099, -5.57615, 0.109409, -0.000869664};
    double reg3_min_sec6_inb_loose[] = {20.9702, -29.641, 3.08715, -0.0542694, 0.00036852};
    double reg3_max_sec6_inb_loose[] = {-36.3081, 47.6882, -6.51571, 0.152934, -0.00153014};
    
    // tight (65 %):
    double reg1_min_sec1_inb_tight[] = {10.8127, 0.867417, -5.23315, 0.236924, 0.236924};
    double reg1_max_sec1_inb_tight[] = {-1.05905, -11.11, 6.99604, -0.280718, 0.00394115};
    double reg2_min_sec1_inb_tight[] = {15.6429, -8.51554, -2.8828, 0.158688, -0.00249913};
    double reg2_max_sec1_inb_tight[] = {4.94591, -19.5039, 8.81259, -0.33865, 0.00465444};
    double reg3_min_sec1_inb_tight[] = {1.34092, 9.28247, -6.60044, 0.281318, -0.00409396};
    double reg3_max_sec1_inb_tight[] = {7.95213, -22.3926, 9.38108, -0.362289, 0.00504018};
    double reg1_min_sec2_inb_tight[] = {6.83952, 1.15885, -4.75388, 0.212145, 0.212145};
    double reg1_max_sec2_inb_tight[] = {-13.5338, 2.90293, 4.04606, -0.18406, 0.00270977};
    double reg2_min_sec2_inb_tight[] = {2.23926, 6.12481, -5.72273, 0.242722, -0.00348261};
    double reg2_max_sec2_inb_tight[] = {-19.2228, 12.2586, 1.80921, -0.109461, 0.00172974};
    double reg3_min_sec2_inb_tight[] = {-7.44705, 15.9478, -7.51342, 0.2998, -0.00422782};
    double reg3_max_sec2_inb_tight[] = {-8.05779, -1.88087, 4.79989, -0.208209, 0.00301469};
    double reg1_min_sec3_inb_tight[] = {5.36591, 8.06724, -6.39578, 0.262518, 0.262518};
    double reg1_max_sec3_inb_tight[] = {-8.20889, -5.7771, 5.98507, -0.245252, 0.00348177};
    double reg2_min_sec3_inb_tight[] = {-2.61308, 17.4864, -8.28425, 0.321424, -0.00449555};
    double reg2_max_sec3_inb_tight[] = {-15.7028, 5.06188, 3.58882, -0.169598, 0.00252206};
    double reg3_min_sec3_inb_tight[] = {-1.47028, 12.199, -6.67515, 0.261044, -0.00361062};
    double reg3_max_sec3_inb_tight[] = {-2.49195, -10.9866, 6.77373, -0.267217, 0.0037212};
    double reg1_min_sec4_inb_tight[] = {7.28085, 4.43634, -5.4954, 0.22469, 0.22469};
    double reg1_max_sec4_inb_tight[] = {-2.88726, -7.75256, 6.11348, -0.246024, 0.00342};
    double reg2_min_sec4_inb_tight[] = {11.1628, -0.875717, -4.40181, 0.191682, -0.00269305};
    double reg2_max_sec4_inb_tight[] = {4.98009, -20.3121, 9.08305, -0.347362, 0.00476636};
    double reg3_min_sec4_inb_tight[] = {-1.53387, 16.9129, -8.4497, 0.334303, -0.00465195};
    double reg3_max_sec4_inb_tight[] = {5.76932, -17.8998, 8.19437, -0.317309, 0.00436422};
    double reg1_min_sec5_inb_tight[] = {16.2619, -7.2257, -3.02427, 0.151797, 0.151797};
    double reg1_max_sec5_inb_tight[] = {-8.7963, -3.03534, 5.17438, -0.214586, 0.0030239};
    double reg2_min_sec5_inb_tight[] = {14.656, -7.43444, -2.74998, 0.141668, -0.00216617};
    double reg2_max_sec5_inb_tight[] = {2.24964, -18.1672, 8.48444, -0.321566, 0.00438927};
    double reg3_min_sec5_inb_tight[] = {3.38717, 7.1616, -5.73249, 0.231652, -0.00320412};
    double reg3_max_sec5_inb_tight[] = {4.63583, -20.0558, 8.78226, -0.33313, 0.0045258};
    double reg1_min_sec6_inb_tight[] = {8.74274, 0.225052, -4.62671, 0.205157, 0.205157};
    double reg1_max_sec6_inb_tight[] = {-8.82084, 1.73358, 3.85433, -0.168173, 0.0024055};
    double reg2_min_sec6_inb_tight[] = {3.86614, 7.51372, -6.36287, 0.268641, -0.00385432};
    double reg2_max_sec6_inb_tight[] = {-7.0834, -0.809866, 4.46815, -0.18982, 0.00266762};
    double reg3_min_sec6_inb_tight[] = {3.12244, 3.84008, -5.05955, 0.223398, -0.00326745};
    double reg3_max_sec6_inb_tight[] = {5.70817, -16.9736, 7.8074, -0.295122, 0.00398425};
    
    // outbending cut parameters (not adjusted, same as inbending):
    
    //_medium (50 %):
    double reg1_min_sec1_outb[] = {25.6464, -26.0938, 1.36205, 0.0083525, 0.0083525};
    double reg1_max_sec1_outb[] = {-21.002, 19.3141, 0.183817, -0.061424, 0.00117881};
    double reg2_min_sec1_outb[] = {20.8156, -19.1708, -0.244432, 0.0652386, -0.00120129};
    double reg2_max_sec1_outb[] = {-20.2319, 19.2371, 0.0491417, -0.0527963, 0.000995253};
    double reg3_min_sec1_outb[] = {22.9342, -24.6621, 1.25766, 0.0133507, -0.00049441};
    double reg3_max_sec1_outb[] = {-8.72949, 4.50811, 3.13299, -0.15452, 0.00234593};
    double reg1_min_sec2_outb[] = {16.0752, -15.8399, -0.65495, 0.0740595, 0.0740595};
    double reg1_max_sec2_outb[] = {-40.921, 40.2558, -3.67186, 0.0498104, -0.000151358};
    double reg2_min_sec2_outb[] = {20.5524, -23.7002, 1.20081, 0.0126066, -0.000501927};
    double reg2_max_sec2_outb[] = {-22.7933, 18.397, 0.707315, -0.0848617, 0.00150965};
    double reg3_min_sec2_outb[] = {7.84857, -9.62769, -1.39224, 0.0876156, -0.00137369};
    double reg3_max_sec2_outb[] = {-24.5099, 24.8373, -1.22162, -0.0135951, 0.00051455};
    double reg1_min_sec3_outb[] = {18.2067, -12.322, -1.79675, 0.111661, 0.111661};
    double reg1_max_sec3_outb[] = {-20.4639, 11.3285, 2.56498, -0.149755, 0.00243423};
    double reg2_min_sec3_outb[] = {17.5093, -13.0897, -1.44801, 0.0995944, -0.00163973};
    double reg2_max_sec3_outb[] = {-29.4241, 25.1723, -0.643032, -0.0433457, 0.0010157};
    double reg3_min_sec3_outb[] = {19.6412, -18.3363, -0.139482, 0.0623644, -0.00122319};
    double reg3_max_sec3_outb[] = {-23.4069, 20.1087, 0.0111006, -0.0595069, 0.00123094};
    double reg1_min_sec4_outb[] = {27.7329, -24.8668, 0.855257, 0.0229504, 0.0229504};
    double reg1_max_sec4_outb[] = {-25.1595, 24.3626, -0.826633, -0.0280833, 0.000694613};
    double reg2_min_sec4_outb[] = {35.3061, -35.9815, 3.34709, -0.0587318, 0.000517383};
    double reg2_max_sec4_outb[] = {-23.5005, 23.5232, -0.872864, -0.0214562, 0.000570724};
    double reg3_min_sec4_outb[] = {13.1517, -8.73139, -2.18634, 0.115133, -0.00170354};
    double reg3_max_sec4_outb[] = {-14.1635, 13.6887, 0.949846, -0.0787818, 0.00128683};
    double reg1_min_sec5_outb[] = {23.5653, -17.285, -1.21053, 0.108096, 0.108096};
    double reg1_max_sec5_outb[] = {-23.3017, 16.1351, 1.45328, -0.112999, 0.00193088};
    double reg2_min_sec5_outb[] = {30.5631, -30.3887, 2.11072, -0.00879513, -0.000281135};
    double reg2_max_sec5_outb[] = {-18.7762, 11.5513, 2.23634, -0.132, 0.00208225};
    double reg3_min_sec5_outb[] = {14.2375, -10.0129, -2.03418, 0.1236, -0.00199283};
    double reg3_max_sec5_outb[] = {-23.8159, 21.7783, -0.373554, -0.0463973, 0.0010055};
    double reg1_min_sec6_outb[] = {34.9766, -38.2253, 3.80455, -0.0624294, -0.0624294};
    double reg1_max_sec6_outb[] = {-23.6287, 21.0793, 0.0916929, -0.0641014, 0.00125089};
    double reg2_min_sec6_outb[] = {21.6462, -23.9106, 1.1985, 0.0108284, -0.000435214};
    double reg2_max_sec6_outb[] = {-22.1881, 21.3921, -0.325758, -0.0406509, 0.000821165};
    double reg3_min_sec6_outb[] = {12.3393, -12.5734, -1.07645, 0.0836258, -0.00137295};
    double reg3_max_sec6_outb[] = {-20.6041, 23.7403, -1.43195, 0.00507544, 0.000160285};
    
    // loose (35 %):
    double reg1_min_sec1_outb_loose[] = {24.9198, -26.6682, 1.22328, 0.0326184, 0.0326184};
    double reg1_max_sec1_outb_loose[] = {-38.0903, 44.8413, -5.13684, 0.0913599, -0.000566229};
    double reg2_min_sec1_outb_loose[] = {34.3872, -41.6892, 4.78076, -0.0904016, 0.000691157};
    double reg2_max_sec1_outb_loose[] = {-36.582, 45.4679, -5.70527, 0.121105, -0.00107292};
    double reg3_min_sec1_outb_loose[] = {38.5882, -51.7707, 7.65247, -0.198183, 0.00222671};
    double reg3_max_sec1_outb_loose[] = {-35.1098, 46.9821, -6.62713, 0.16512, -0.00177347};
    double reg1_min_sec2_outb_loose[] = {26.3964, -32.5314, 2.87032, -0.0261837, -0.0261837};
    double reg1_max_sec2_outb_loose[] = {-45.4528, 53.4467, -6.83895, 0.144609, -0.00127358};
    double reg2_min_sec2_outb_loose[] = {29.5885, -38.6103, 4.45494, -0.084114, 0.000626459};
    double reg2_max_sec2_outb_loose[] = {-42.6042, 50.783, -6.44381, 0.135305, -0.00117332};
    double reg3_min_sec2_outb_loose[] = {21.3286, -33.2316, 4.16718, -0.0947928, 0.000950229};
    double reg3_max_sec2_outb_loose[] = {-26.2955, 32.832, -3.21406, 0.0443669, -0.000110018};
    double reg1_min_sec3_outb_loose[] = {38.3683, -39.555, 3.48786, -0.0311045, -0.0311045};
    double reg1_max_sec3_outb_loose[] = {-44.7002, 48.5238, -5.33933, 0.0857334, -0.000389759};
    double reg2_min_sec3_outb_loose[] = {26.2411, -26.4956, 1.17454, 0.0314168, -0.000923764};
    double reg2_max_sec3_outb_loose[] = {-49.1042, 53.7936, -6.50797, 0.125435, -0.000940508};
    double reg3_min_sec3_outb_loose[] = {30.2087, -36.2303, 3.94272, -0.070659, 0.000484111};
    double reg3_max_sec3_outb_loose[] = {-25.1075, 27.5492, -1.81589, -0.00372993, 0.000534551};
    double reg1_min_sec4_outb_loose[] = {45.5355, -50.2561, 5.91736, -0.116203, -0.116203};
    double reg1_max_sec4_outb_loose[] = {-45.3956, 52.8208, -6.57583, 0.13245, -0.00106037};
    double reg2_min_sec4_outb_loose[] = {42.9007, -48.8609, 6.02101, -0.13073, 0.00125121};
    double reg2_max_sec4_outb_loose[] = {-32.1523, 38.6933, -4.15879, 0.0722104, -0.000474606};
    double reg3_min_sec4_outb_loose[] = {34.0991, -40.7811, 4.8997, -0.107033, 0.00105086};
    double reg3_max_sec4_outb_loose[] = {-39.189, 52.9641, -8.01756, 0.214872, -0.00250226};
    double reg1_min_sec5_outb_loose[] = {52.0273, -61.3244, 8.45187, -0.193011, -0.193011};
    double reg1_max_sec5_outb_loose[] = {-50.6284, 54.3673, -6.38031, 0.116905, -0.000778886};
    double reg2_min_sec5_outb_loose[] = {32.1495, -35.3084, 3.11177, -0.0314703, -8.39833e-05};
    double reg2_max_sec5_outb_loose[] = {-47.0576, 50.676, -5.80697, 0.10389, -0.000677858};
    double reg3_min_sec5_outb_loose[] = {27.5659, -34.8935, 3.81989, -0.0655272, 0.000389246};
    double reg3_max_sec5_outb_loose[] = {-25.7138, 27.0011, -1.53335, -0.0152295, 0.000694336};
    double reg1_min_sec6_outb_loose[] = {35.6586, -43.3351, 5.05374, -0.0961372, -0.0961372};
    double reg1_max_sec6_outb_loose[] = {-43.6986, 52.1413, -6.61602, 0.137247, -0.00115149};
    double reg2_min_sec6_outb_loose[] = {33.3586, -42.3939, 5.17582, -0.106612, 0.000910642};
    double reg2_max_sec6_outb_loose[] = {-38.3305, 46.3099, -5.57615, 0.109409, -0.000869664};
    double reg3_min_sec6_outb_loose[] = {20.9702, -29.641, 3.08715, -0.0542694, 0.00036852};
    double reg3_max_sec6_outb_loose[] = {-36.3081, 47.6882, -6.51571, 0.152934, -0.00153014};
    
    // tight (65 %):
    double reg1_min_sec1_outb_tight[] = {10.8127, 0.867417, -5.23315, 0.236924, 0.236924};
    double reg1_max_sec1_outb_tight[] = {-1.05905, -11.11, 6.99604, -0.280718, 0.00394115};
    double reg2_min_sec1_outb_tight[] = {15.6429, -8.51554, -2.8828, 0.158688, -0.00249913};
    double reg2_max_sec1_outb_tight[] = {4.94591, -19.5039, 8.81259, -0.33865, 0.00465444};
    double reg3_min_sec1_outb_tight[] = {1.34092, 9.28247, -6.60044, 0.281318, -0.00409396};
    double reg3_max_sec1_outb_tight[] = {7.95213, -22.3926, 9.38108, -0.362289, 0.00504018};
    double reg1_min_sec2_outb_tight[] = {6.83952, 1.15885, -4.75388, 0.212145, 0.212145};
    double reg1_max_sec2_outb_tight[] = {-13.5338, 2.90293, 4.04606, -0.18406, 0.00270977};
    double reg2_min_sec2_outb_tight[] = {2.23926, 6.12481, -5.72273, 0.242722, -0.00348261};
    double reg2_max_sec2_outb_tight[] = {-19.2228, 12.2586, 1.80921, -0.109461, 0.00172974};
    double reg3_min_sec2_outb_tight[] = {-7.44705, 15.9478, -7.51342, 0.2998, -0.00422782};
    double reg3_max_sec2_outb_tight[] = {-8.05779, -1.88087, 4.79989, -0.208209, 0.00301469};
    double reg1_min_sec3_outb_tight[] = {5.36591, 8.06724, -6.39578, 0.262518, 0.262518};
    double reg1_max_sec3_outb_tight[] = {-8.20889, -5.7771, 5.98507, -0.245252, 0.00348177};
    double reg2_min_sec3_outb_tight[] = {-2.61308, 17.4864, -8.28425, 0.321424, -0.00449555};
    double reg2_max_sec3_outb_tight[] = {-15.7028, 5.06188, 3.58882, -0.169598, 0.00252206};
    double reg3_min_sec3_outb_tight[] = {-1.47028, 12.199, -6.67515, 0.261044, -0.00361062};
    double reg3_max_sec3_outb_tight[] = {-2.49195, -10.9866, 6.77373, -0.267217, 0.0037212};
    double reg1_min_sec4_outb_tight[] = {7.28085, 4.43634, -5.4954, 0.22469, 0.22469};
    double reg1_max_sec4_outb_tight[] = {-2.88726, -7.75256, 6.11348, -0.246024, 0.00342};
    double reg2_min_sec4_outb_tight[] = {11.1628, -0.875717, -4.40181, 0.191682, -0.00269305};
    double reg2_max_sec4_outb_tight[] = {4.98009, -20.3121, 9.08305, -0.347362, 0.00476636};
    double reg3_min_sec4_outb_tight[] = {-1.53387, 16.9129, -8.4497, 0.334303, -0.00465195};
    double reg3_max_sec4_outb_tight[] = {5.76932, -17.8998, 8.19437, -0.317309, 0.00436422};
    double reg1_min_sec5_outb_tight[] = {16.2619, -7.2257, -3.02427, 0.151797, 0.151797};
    double reg1_max_sec5_outb_tight[] = {-8.7963, -3.03534, 5.17438, -0.214586, 0.0030239};
    double reg2_min_sec5_outb_tight[] = {14.656, -7.43444, -2.74998, 0.141668, -0.00216617};
    double reg2_max_sec5_outb_tight[] = {2.24964, -18.1672, 8.48444, -0.321566, 0.00438927};
    double reg3_min_sec5_outb_tight[] = {3.38717, 7.1616, -5.73249, 0.231652, -0.00320412};
    double reg3_max_sec5_outb_tight[] = {4.63583, -20.0558, 8.78226, -0.33313, 0.0045258};
    double reg1_min_sec6_outb_tight[] = {8.74274, 0.225052, -4.62671, 0.205157, 0.205157};
    double reg1_max_sec6_outb_tight[] = {-8.82084, 1.73358, 3.85433, -0.168173, 0.0024055};
    double reg2_min_sec6_outb_tight[] = {3.86614, 7.51372, -6.36287, 0.268641, -0.00385432};
    double reg2_max_sec6_outb_tight[] = {-7.0834, -0.809866, 4.46815, -0.18982, 0.00266762};
    double reg3_min_sec6_outb_tight[] = {3.12244, 3.84008, -5.05955, 0.223398, -0.00326745};
    double reg3_max_sec6_outb_tight[] = {5.70817, -16.9736, 7.8074, -0.295122, 0.00398425};
    
    
    double p0_min = 0; double p1_min = 0; double p2_min = 0; double p3_min = 0; double p4_min = 0;
    double p0_max = 0; double p1_max = 0; double p2_max = 0; double p3_max = 0; double p4_max = 0;
    
    
    if(region == 1){
        if(part->trk(DC)->getSector() == 1 && inbending == true){
            p0_min = reg1_min_sec1_inb[0]; p1_min = reg1_min_sec1_inb[1]; p2_min = reg1_min_sec1_inb[2]; p3_min = reg1_min_sec1_inb[3]; p4_min = reg1_min_sec1_inb[4];
            p0_max = reg1_max_sec1_inb[0]; p1_max = reg1_max_sec1_inb[1]; p2_max = reg1_max_sec1_inb[2]; p3_max = reg1_max_sec1_inb[3]; p4_max = reg1_max_sec1_inb[4];
        }
        if(part->trk(DC)->getSector() == 2 && inbending == true){
            p0_min = reg1_min_sec2_inb[0]; p1_min = reg1_min_sec2_inb[1]; p2_min = reg1_min_sec2_inb[2]; p3_min = reg1_min_sec2_inb[3]; p4_min = reg1_min_sec2_inb[4];
            p0_max = reg1_max_sec2_inb[0]; p1_max = reg1_max_sec2_inb[1]; p2_max = reg1_max_sec2_inb[2]; p3_max = reg1_max_sec2_inb[3]; p4_max = reg1_max_sec2_inb[4];
        }
        if(part->trk(DC)->getSector() == 3 && inbending == true){
            p0_min = reg1_min_sec3_inb[0]; p1_min = reg1_min_sec3_inb[1]; p2_min = reg1_min_sec3_inb[2]; p3_min = reg1_min_sec3_inb[3]; p4_min = reg1_min_sec3_inb[4];
            p0_max = reg1_max_sec3_inb[0]; p1_max = reg1_max_sec3_inb[1]; p2_max = reg1_max_sec3_inb[2]; p3_max = reg1_max_sec3_inb[3]; p4_max = reg1_max_sec3_inb[4];
        }
        if(part->trk(DC)->getSector() == 4 && inbending == true){
            p0_min = reg1_min_sec4_inb[0]; p1_min = reg1_min_sec4_inb[1]; p2_min = reg1_min_sec4_inb[2]; p3_min = reg1_min_sec4_inb[3]; p4_min = reg1_min_sec4_inb[4];
            p0_max = reg1_max_sec4_inb[0]; p1_max = reg1_max_sec4_inb[1]; p2_max = reg1_max_sec4_inb[2]; p3_max = reg1_max_sec4_inb[3]; p4_max = reg1_max_sec4_inb[4];
        }
        if(part->trk(DC)->getSector() == 5 && inbending == true){
            p0_min = reg1_min_sec5_inb[0]; p1_min = reg1_min_sec5_inb[1]; p2_min = reg1_min_sec5_inb[2]; p3_min = reg1_min_sec5_inb[3]; p4_min = reg1_min_sec5_inb[4];
            p0_max = reg1_max_sec5_inb[0]; p1_max = reg1_max_sec5_inb[1]; p2_max = reg1_max_sec5_inb[2]; p3_max = reg1_max_sec5_inb[3]; p4_max = reg1_max_sec5_inb[4];
        }
        if(part->trk(DC)->getSector() == 6 && inbending == true){
            p0_min = reg1_min_sec6_inb[0]; p1_min = reg1_min_sec6_inb[1]; p2_min = reg1_min_sec6_inb[2]; p3_min = reg1_min_sec6_inb[3]; p4_min = reg1_min_sec6_inb[4];
            p0_max = reg1_max_sec6_inb[0]; p1_max = reg1_max_sec6_inb[1]; p2_max = reg1_max_sec6_inb[2]; p3_max = reg1_max_sec6_inb[3]; p4_max = reg1_max_sec6_inb[4];
        }
        if(part->trk(DC)->getSector() == 1 && outbending == true){
            p0_min = reg1_min_sec1_outb[0]; p1_min = reg1_min_sec1_outb[1]; p2_min = reg1_min_sec1_outb[2]; p3_min = reg1_min_sec1_outb[3]; p4_min = reg1_min_sec1_outb[4];
            p0_max = reg1_max_sec1_outb[0]; p1_max = reg1_max_sec1_outb[1]; p2_max = reg1_max_sec1_outb[2]; p3_max = reg1_max_sec1_outb[3]; p4_max = reg1_max_sec1_outb[4];
        }
        if(part->trk(DC)->getSector() == 2 && outbending == true){
            p0_min = reg1_min_sec2_outb[0]; p1_min = reg1_min_sec2_outb[1]; p2_min = reg1_min_sec2_outb[2]; p3_min = reg1_min_sec2_outb[3]; p4_min = reg1_min_sec2_outb[4];
            p0_max = reg1_max_sec2_outb[0]; p1_max = reg1_max_sec2_outb[1]; p2_max = reg1_max_sec2_outb[2]; p3_max = reg1_max_sec2_outb[3]; p4_max = reg1_max_sec2_outb[4];
        }
        if(part->trk(DC)->getSector() == 3 && outbending == true){
            p0_min = reg1_min_sec3_outb[0]; p1_min = reg1_min_sec3_outb[1]; p2_min = reg1_min_sec3_outb[2]; p3_min = reg1_min_sec3_outb[3]; p4_min = reg1_min_sec3_outb[4];
            p0_max = reg1_max_sec3_outb[0]; p1_max = reg1_max_sec3_outb[1]; p2_max = reg1_max_sec3_outb[2]; p3_max = reg1_max_sec3_outb[3]; p4_max = reg1_max_sec3_outb[4];
        }
        if(part->trk(DC)->getSector() == 4 && outbending == true){
            p0_min = reg1_min_sec4_outb[0]; p1_min = reg1_min_sec4_outb[1]; p2_min = reg1_min_sec4_outb[2]; p3_min = reg1_min_sec4_outb[3]; p4_min = reg1_min_sec4_outb[4];
            p0_max = reg1_max_sec4_outb[0]; p1_max = reg1_max_sec4_outb[1]; p2_max = reg1_max_sec4_outb[2]; p3_max = reg1_max_sec4_outb[3]; p4_max = reg1_max_sec4_outb[4];
        }
        if(part->trk(DC)->getSector() == 5 && outbending == true){
            p0_min = reg1_min_sec5_outb[0]; p1_min = reg1_min_sec5_outb[1]; p2_min = reg1_min_sec5_outb[2]; p3_min = reg1_min_sec5_outb[3]; p4_min = reg1_min_sec5_outb[4];
            p0_max = reg1_max_sec5_outb[0]; p1_max = reg1_max_sec5_outb[1]; p2_max = reg1_max_sec5_outb[2]; p3_max = reg1_max_sec5_outb[3]; p4_max = reg1_max_sec5_outb[4];
        }
        if(part->trk(DC)->getSector() == 6 && outbending == true){
            p0_min = reg1_min_sec6_outb[0]; p1_min = reg1_min_sec6_outb[1]; p2_min = reg1_min_sec6_outb[2]; p3_min = reg1_min_sec6_outb[3]; p4_min = reg1_min_sec6_outb[4];
            p0_max = reg1_max_sec6_outb[0]; p1_max = reg1_max_sec6_outb[1]; p2_max = reg1_max_sec6_outb[2]; p3_max = reg1_max_sec6_outb[3]; p4_max = reg1_max_sec6_outb[4];
        }
    }
    
    if(region == 2){
        if(part->trk(DC)->getSector()  == 1 && inbending == true){
            p0_min = reg2_min_sec1_inb[0]; p1_min = reg2_min_sec1_inb[1]; p2_min = reg2_min_sec1_inb[2]; p3_min = reg2_min_sec1_inb[3]; p4_min = reg2_min_sec1_inb[4];
            p0_max = reg2_max_sec1_inb[0]; p1_max = reg2_max_sec1_inb[1]; p2_max = reg2_max_sec1_inb[2]; p3_max = reg2_max_sec1_inb[3]; p4_max = reg2_max_sec1_inb[4];
        }
        if(part->trk(DC)->getSector() ==  2 && inbending == true){
            p0_min = reg2_min_sec2_inb[0]; p1_min = reg2_min_sec2_inb[1]; p2_min = reg2_min_sec2_inb[2]; p3_min = reg2_min_sec2_inb[3]; p4_min = reg2_min_sec2_inb[4];
            p0_max = reg2_max_sec2_inb[0]; p1_max = reg2_max_sec2_inb[1]; p2_max = reg2_max_sec2_inb[2]; p3_max = reg2_max_sec2_inb[3]; p4_max = reg2_max_sec2_inb[4];
        }
        if(part->trk(DC)->getSector() ==  3 && inbending == true){
            p0_min = reg2_min_sec3_inb[0]; p1_min = reg2_min_sec3_inb[1]; p2_min = reg2_min_sec3_inb[2]; p3_min = reg2_min_sec3_inb[3]; p4_min = reg2_min_sec3_inb[4];
            p0_max = reg2_max_sec3_inb[0]; p1_max = reg2_max_sec3_inb[1]; p2_max = reg2_max_sec3_inb[2]; p3_max = reg2_max_sec3_inb[3]; p4_max = reg2_max_sec3_inb[4];
        }
        if(part->trk(DC)->getSector() ==  4 && inbending == true){
            p0_min = reg2_min_sec4_inb[0]; p1_min = reg2_min_sec4_inb[1]; p2_min = reg2_min_sec4_inb[2]; p3_min = reg2_min_sec4_inb[3]; p4_min = reg2_min_sec4_inb[4];
            p0_max = reg2_max_sec4_inb[0]; p1_max = reg2_max_sec4_inb[1]; p2_max = reg2_max_sec4_inb[2]; p3_max = reg2_max_sec4_inb[3]; p4_max = reg2_max_sec4_inb[4];
        }
        if(part->trk(DC)->getSector() ==  5 && inbending == true){
            p0_min = reg2_min_sec5_inb[0]; p1_min = reg2_min_sec5_inb[1]; p2_min = reg2_min_sec5_inb[2]; p3_min = reg2_min_sec5_inb[3]; p4_min = reg2_min_sec5_inb[4];
            p0_max = reg2_max_sec5_inb[0]; p1_max = reg2_max_sec5_inb[1]; p2_max = reg2_max_sec5_inb[2]; p3_max = reg2_max_sec5_inb[3]; p4_max = reg2_max_sec5_inb[4];
        }
        if(part->trk(DC)->getSector() == 6 && inbending == true){
            p0_min = reg2_min_sec6_inb[0]; p1_min = reg2_min_sec6_inb[1]; p2_min = reg2_min_sec6_inb[2]; p3_min = reg2_min_sec6_inb[3]; p4_min = reg2_min_sec6_inb[4];
            p0_max = reg2_max_sec6_inb[0]; p1_max = reg2_max_sec6_inb[1]; p2_max = reg2_max_sec6_inb[2]; p3_max = reg2_max_sec6_inb[3]; p4_max = reg2_max_sec6_inb[4];
        }
        if(part->trk(DC)->getSector() == 1 && outbending == true){
            p0_min = reg2_min_sec1_outb[0]; p1_min = reg2_min_sec1_outb[1]; p2_min = reg2_min_sec1_outb[2]; p3_min = reg2_min_sec1_outb[3]; p4_min = reg2_min_sec1_outb[4];
            p0_max = reg2_max_sec1_outb[0]; p1_max = reg2_max_sec1_outb[1]; p2_max = reg2_max_sec1_outb[2]; p3_max = reg2_max_sec1_outb[3]; p4_max = reg2_max_sec1_outb[4];
        }
        if(part->trk(DC)->getSector() ==  2 && outbending == true){
            p0_min = reg2_min_sec2_outb[0]; p1_min = reg2_min_sec2_outb[1]; p2_min = reg2_min_sec2_outb[2]; p3_min = reg2_min_sec2_outb[3]; p4_min = reg2_min_sec2_outb[4];
            p0_max = reg2_max_sec2_outb[0]; p1_max = reg2_max_sec2_outb[1]; p2_max = reg2_max_sec2_outb[2]; p3_max = reg2_max_sec2_outb[3]; p4_max = reg2_max_sec2_outb[4];
        }
        if(part->trk(DC)->getSector() == 3 && outbending == true){
            p0_min = reg2_min_sec3_outb[0]; p1_min = reg2_min_sec3_outb[1]; p2_min = reg2_min_sec3_outb[2]; p3_min = reg2_min_sec3_outb[3]; p4_min = reg2_min_sec3_outb[4];
            p0_max = reg2_max_sec3_outb[0]; p1_max = reg2_max_sec3_outb[1]; p2_max = reg2_max_sec3_outb[2]; p3_max = reg2_max_sec3_outb[3]; p4_max = reg2_max_sec3_outb[4];
        }
        if(part->trk(DC)->getSector() == 4 && outbending == true){
            p0_min = reg2_min_sec4_outb[0]; p1_min = reg2_min_sec4_outb[1]; p2_min = reg2_min_sec4_outb[2]; p3_min = reg2_min_sec4_outb[3]; p4_min = reg2_min_sec4_outb[4];
            p0_max = reg2_max_sec4_outb[0]; p1_max = reg2_max_sec4_outb[1]; p2_max = reg2_max_sec4_outb[2]; p3_max = reg2_max_sec4_outb[3]; p4_max = reg2_max_sec4_outb[4];
        }
        if(part->trk(DC)->getSector() == 5 && outbending == true){
            p0_min = reg2_min_sec5_outb[0]; p1_min = reg2_min_sec5_outb[1]; p2_min = reg2_min_sec5_outb[2]; p3_min = reg2_min_sec5_outb[3]; p4_min = reg2_min_sec5_outb[4];
            p0_max = reg2_max_sec5_outb[0]; p1_max = reg2_max_sec5_outb[1]; p2_max = reg2_max_sec5_outb[2]; p3_max = reg2_max_sec5_outb[3]; p4_max = reg2_max_sec5_outb[4];
        }
        if(part->trk(DC)->getSector() == 6 && outbending == true){
            p0_min = reg2_min_sec6_outb[0]; p1_min = reg2_min_sec6_outb[1]; p2_min = reg2_min_sec6_outb[2]; p3_min = reg2_min_sec6_outb[3]; p4_min = reg2_min_sec6_outb[4];
            p0_max = reg2_max_sec6_outb[0]; p1_max = reg2_max_sec6_outb[1]; p2_max = reg2_max_sec6_outb[2]; p3_max = reg2_max_sec6_outb[3]; p4_max = reg2_max_sec6_outb[4];
        }
    }
    
    if(region == 3){
        if(part->trk(DC)->getSector()  == 1 && inbending == true){
            p0_min = reg3_min_sec1_inb[0]; p1_min = reg3_min_sec1_inb[1]; p2_min = reg3_min_sec1_inb[2]; p3_min = reg3_min_sec1_inb[3]; p4_min = reg3_min_sec1_inb[4];
            p0_max = reg3_max_sec1_inb[0]; p1_max = reg3_max_sec1_inb[1]; p2_max = reg3_max_sec1_inb[2]; p3_max = reg3_max_sec1_inb[3]; p4_max = reg3_max_sec1_inb[4];
        }
        if(part->trk(DC)->getSector() ==  2 && inbending == true){
            p0_min = reg3_min_sec2_inb[0]; p1_min = reg3_min_sec2_inb[1]; p2_min = reg3_min_sec2_inb[2]; p3_min = reg3_min_sec2_inb[3]; p4_min = reg3_min_sec2_inb[4];
            p0_max = reg3_max_sec2_inb[0]; p1_max = reg3_max_sec2_inb[1]; p2_max = reg3_max_sec2_inb[2]; p3_max = reg3_max_sec2_inb[3]; p4_max = reg3_max_sec2_inb[4];
        }
        if(part->trk(DC)->getSector() ==  3 && inbending == true){
            p0_min = reg3_min_sec3_inb[0]; p1_min = reg3_min_sec3_inb[1]; p2_min = reg3_min_sec3_inb[2]; p3_min = reg3_min_sec3_inb[3]; p4_min = reg3_min_sec3_inb[4];
            p0_max = reg3_max_sec3_inb[0]; p1_max = reg3_max_sec3_inb[1]; p2_max = reg3_max_sec3_inb[2]; p3_max = reg3_max_sec3_inb[3]; p4_max = reg3_max_sec3_inb[4];
        }
        if(part->trk(DC)->getSector() == 4 && inbending == true){
            p0_min = reg3_min_sec4_inb[0]; p1_min = reg3_min_sec4_inb[1]; p2_min = reg3_min_sec4_inb[2]; p3_min = reg3_min_sec4_inb[3]; p4_min = reg3_min_sec4_inb[4];
            p0_max = reg3_max_sec4_inb[0]; p1_max = reg3_max_sec4_inb[1]; p2_max = reg3_max_sec4_inb[2]; p3_max = reg3_max_sec4_inb[3]; p4_max = reg3_max_sec4_inb[4];
        }
        if(part->trk(DC)->getSector() == 5 && inbending == true){
            p0_min = reg3_min_sec5_inb[0]; p1_min = reg3_min_sec5_inb[1]; p2_min = reg3_min_sec5_inb[2]; p3_min = reg3_min_sec5_inb[3]; p4_min = reg3_min_sec5_inb[4];
            p0_max = reg3_max_sec5_inb[0]; p1_max = reg3_max_sec5_inb[1]; p2_max = reg3_max_sec5_inb[2]; p3_max = reg3_max_sec5_inb[3]; p4_max = reg3_max_sec5_inb[4];
        }
        if(part->trk(DC)->getSector() == 6 && inbending == true){
            p0_min = reg3_min_sec6_inb[0]; p1_min = reg3_min_sec6_inb[1]; p2_min = reg3_min_sec6_inb[2]; p3_min = reg3_min_sec6_inb[3]; p4_min = reg3_min_sec6_inb[4];
            p0_max = reg3_max_sec6_inb[0]; p1_max = reg3_max_sec6_inb[1]; p2_max = reg3_max_sec6_inb[2]; p3_max = reg3_max_sec6_inb[3]; p4_max = reg3_max_sec6_inb[4];
        }
        if(part->trk(DC)->getSector() ==  1 && outbending == true){
            p0_min = reg3_min_sec1_outb[0]; p1_min = reg3_min_sec1_outb[1]; p2_min = reg3_min_sec1_outb[2]; p3_min = reg3_min_sec1_outb[3]; p4_min = reg3_min_sec1_outb[4];
            p0_max = reg3_max_sec1_outb[0]; p1_max = reg3_max_sec1_outb[1]; p2_max = reg3_max_sec1_outb[2]; p3_max = reg3_max_sec1_outb[3]; p4_max = reg3_max_sec1_outb[4];
        }
        if(part->trk(DC)->getSector() == 2 && outbending == true){
            p0_min = reg3_min_sec2_outb[0]; p1_min = reg3_min_sec2_outb[1]; p2_min = reg3_min_sec2_outb[2]; p3_min = reg3_min_sec2_outb[3]; p4_min = reg3_min_sec2_outb[4];
            p0_max = reg3_max_sec2_outb[0]; p1_max = reg3_max_sec2_outb[1]; p2_max = reg3_max_sec2_outb[2]; p3_max = reg3_max_sec2_outb[3]; p4_max = reg3_max_sec2_outb[4];
        }
        if(part->trk(DC)->getSector() == 3 && outbending == true){
            p0_min = reg3_min_sec3_outb[0]; p1_min = reg3_min_sec3_outb[1]; p2_min = reg3_min_sec3_outb[2]; p3_min = reg3_min_sec3_outb[3]; p4_min = reg3_min_sec3_outb[4];
            p0_max = reg3_max_sec3_outb[0]; p1_max = reg3_max_sec3_outb[1]; p2_max = reg3_max_sec3_outb[2]; p3_max = reg3_max_sec3_outb[3]; p4_max = reg3_max_sec3_outb[4];
        }
        if(part->trk(DC)->getSector() == 4 && outbending == true){
            p0_min = reg3_min_sec4_outb[0]; p1_min = reg3_min_sec4_outb[1]; p2_min = reg3_min_sec4_outb[2]; p3_min = reg3_min_sec4_outb[3]; p4_min = reg3_min_sec4_outb[4];
            p0_max = reg3_max_sec4_outb[0]; p1_max = reg3_max_sec4_outb[1]; p2_max = reg3_max_sec4_outb[2]; p3_max = reg3_max_sec4_outb[3]; p4_max = reg3_max_sec4_outb[4];
        }
        if(part->trk(DC)->getSector() == 5 && outbending == true){
            p0_min = reg3_min_sec5_outb[0]; p1_min = reg3_min_sec5_outb[1]; p2_min = reg3_min_sec5_outb[2]; p3_min = reg3_min_sec5_outb[3]; p4_min = reg3_min_sec5_outb[4];
            p0_max = reg3_max_sec5_outb[0]; p1_max = reg3_max_sec5_outb[1]; p2_max = reg3_max_sec5_outb[2]; p3_max = reg3_max_sec5_outb[3]; p4_max = reg3_max_sec5_outb[4];
        }
        if(part->trk(DC)->getSector() ==  6 && outbending == true){
            p0_min = reg3_min_sec6_outb[0]; p1_min = reg3_min_sec6_outb[1]; p2_min = reg3_min_sec6_outb[2]; p3_min = reg3_min_sec6_outb[3]; p4_min = reg3_min_sec6_outb[4];
            p0_max = reg3_max_sec6_outb[0]; p1_max = reg3_max_sec6_outb[1]; p2_max = reg3_max_sec6_outb[2]; p3_max = reg3_max_sec6_outb[3]; p4_max = reg3_max_sec6_outb[4];
        }
    }
    
    if(tight == true){
        if(region == 1){
            if(part->trk(DC)->getSector() == 1 && inbending == true){
                p0_min = reg1_min_sec1_inb_tight[0]; p1_min = reg1_min_sec1_inb_tight[1]; p2_min = reg1_min_sec1_inb_tight[2]; p3_min = reg1_min_sec1_inb_tight[3]; p4_min = reg1_min_sec1_inb_tight[4];
                p0_max = reg1_max_sec1_inb_tight[0]; p1_max = reg1_max_sec1_inb_tight[1]; p2_max = reg1_max_sec1_inb_tight[2]; p3_max = reg1_max_sec1_inb_tight[3]; p4_max = reg1_max_sec1_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 2 && inbending == true){
                p0_min = reg1_min_sec2_inb_tight[0]; p1_min = reg1_min_sec2_inb_tight[1]; p2_min = reg1_min_sec2_inb_tight[2]; p3_min = reg1_min_sec2_inb_tight[3]; p4_min = reg1_min_sec2_inb_tight[4];
                p0_max = reg1_max_sec2_inb_tight[0]; p1_max = reg1_max_sec2_inb_tight[1]; p2_max = reg1_max_sec2_inb_tight[2]; p3_max = reg1_max_sec2_inb_tight[3]; p4_max = reg1_max_sec2_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 3 && inbending == true){
                p0_min = reg1_min_sec3_inb_tight[0]; p1_min = reg1_min_sec3_inb_tight[1]; p2_min = reg1_min_sec3_inb_tight[2]; p3_min = reg1_min_sec3_inb_tight[3]; p4_min = reg1_min_sec3_inb_tight[4];
                p0_max = reg1_max_sec3_inb_tight[0]; p1_max = reg1_max_sec3_inb_tight[1]; p2_max = reg1_max_sec3_inb_tight[2]; p3_max = reg1_max_sec3_inb_tight[3]; p4_max = reg1_max_sec3_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 4 && inbending == true){
                p0_min = reg1_min_sec4_inb_tight[0]; p1_min = reg1_min_sec4_inb_tight[1]; p2_min = reg1_min_sec4_inb_tight[2]; p3_min = reg1_min_sec4_inb_tight[3]; p4_min = reg1_min_sec4_inb_tight[4];
                p0_max = reg1_max_sec4_inb_tight[0]; p1_max = reg1_max_sec4_inb_tight[1]; p2_max = reg1_max_sec4_inb_tight[2]; p3_max = reg1_max_sec4_inb_tight[3]; p4_max = reg1_max_sec4_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 5 && inbending == true){
                p0_min = reg1_min_sec5_inb_tight[0]; p1_min = reg1_min_sec5_inb_tight[1]; p2_min = reg1_min_sec5_inb_tight[2]; p3_min = reg1_min_sec5_inb_tight[3]; p4_min = reg1_min_sec5_inb_tight[4];
                p0_max = reg1_max_sec5_inb_tight[0]; p1_max = reg1_max_sec5_inb_tight[1]; p2_max = reg1_max_sec5_inb_tight[2]; p3_max = reg1_max_sec5_inb_tight[3]; p4_max = reg1_max_sec5_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 6 && inbending == true){
                p0_min = reg1_min_sec6_inb_tight[0]; p1_min = reg1_min_sec6_inb_tight[1]; p2_min = reg1_min_sec6_inb_tight[2]; p3_min = reg1_min_sec6_inb_tight[3]; p4_min = reg1_min_sec6_inb_tight[4];
                p0_max = reg1_max_sec6_inb_tight[0]; p1_max = reg1_max_sec6_inb_tight[1]; p2_max = reg1_max_sec6_inb_tight[2]; p3_max = reg1_max_sec6_inb_tight[3]; p4_max = reg1_max_sec6_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 1 && outbending == true){
                p0_min = reg1_min_sec1_outb_tight[0]; p1_min = reg1_min_sec1_outb_tight[1]; p2_min = reg1_min_sec1_outb_tight[2]; p3_min = reg1_min_sec1_outb_tight[3]; p4_min = reg1_min_sec1_outb_tight[4];
                p0_max = reg1_max_sec1_outb_tight[0]; p1_max = reg1_max_sec1_outb_tight[1]; p2_max = reg1_max_sec1_outb_tight[2]; p3_max = reg1_max_sec1_outb_tight[3]; p4_max = reg1_max_sec1_outb_tight[4];
            }
            if(part->trk(DC)->getSector() ==  2 && outbending == true){
                p0_min = reg1_min_sec2_outb_tight[0]; p1_min = reg1_min_sec2_outb_tight[1]; p2_min = reg1_min_sec2_outb_tight[2]; p3_min = reg1_min_sec2_outb_tight[3]; p4_min = reg1_min_sec2_outb_tight[4];
                p0_max = reg1_max_sec2_outb_tight[0]; p1_max = reg1_max_sec2_outb_tight[1]; p2_max = reg1_max_sec2_outb_tight[2]; p3_max = reg1_max_sec2_outb_tight[3]; p4_max = reg1_max_sec2_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 3 && outbending == true){
                p0_min = reg1_min_sec3_outb_tight[0]; p1_min = reg1_min_sec3_outb_tight[1]; p2_min = reg1_min_sec3_outb_tight[2]; p3_min = reg1_min_sec3_outb_tight[3]; p4_min = reg1_min_sec3_outb_tight[4];
                p0_max = reg1_max_sec3_outb_tight[0]; p1_max = reg1_max_sec3_outb_tight[1]; p2_max = reg1_max_sec3_outb_tight[2]; p3_max = reg1_max_sec3_outb_tight[3]; p4_max = reg1_max_sec3_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 4 && outbending == true){
                p0_min = reg1_min_sec4_outb_tight[0]; p1_min = reg1_min_sec4_outb_tight[1]; p2_min = reg1_min_sec4_outb_tight[2]; p3_min = reg1_min_sec4_outb_tight[3]; p4_min = reg1_min_sec4_outb_tight[4];
                p0_max = reg1_max_sec4_outb_tight[0]; p1_max = reg1_max_sec4_outb_tight[1]; p2_max = reg1_max_sec4_outb_tight[2]; p3_max = reg1_max_sec4_outb_tight[3]; p4_max = reg1_max_sec4_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 5 && outbending == true){
                p0_min = reg1_min_sec5_outb_tight[0]; p1_min = reg1_min_sec5_outb_tight[1]; p2_min = reg1_min_sec5_outb_tight[2]; p3_min = reg1_min_sec5_outb_tight[3]; p4_min = reg1_min_sec5_outb_tight[4];
                p0_max = reg1_max_sec5_outb_tight[0]; p1_max = reg1_max_sec5_outb_tight[1]; p2_max = reg1_max_sec5_outb_tight[2]; p3_max = reg1_max_sec5_outb_tight[3]; p4_max = reg1_max_sec5_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 6 && outbending == true){
                p0_min = reg1_min_sec6_outb_tight[0]; p1_min = reg1_min_sec6_outb_tight[1]; p2_min = reg1_min_sec6_outb_tight[2]; p3_min = reg1_min_sec6_outb_tight[3]; p4_min = reg1_min_sec6_outb_tight[4];
                p0_max = reg1_max_sec6_outb_tight[0]; p1_max = reg1_max_sec6_outb_tight[1]; p2_max = reg1_max_sec6_outb_tight[2]; p3_max = reg1_max_sec6_outb_tight[3]; p4_max = reg1_max_sec6_outb_tight[4];
            }
        }
        
        if(region == 2){
            if(part->trk(DC)->getSector() == 1 && inbending == true){
                p0_min = reg2_min_sec1_inb_tight[0]; p1_min = reg2_min_sec1_inb_tight[1]; p2_min = reg2_min_sec1_inb_tight[2]; p3_min = reg2_min_sec1_inb_tight[3]; p4_min = reg2_min_sec1_inb_tight[4];
                p0_max = reg2_max_sec1_inb_tight[0]; p1_max = reg2_max_sec1_inb_tight[1]; p2_max = reg2_max_sec1_inb_tight[2]; p3_max = reg2_max_sec1_inb_tight[3]; p4_max = reg2_max_sec1_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 2 && inbending == true){
                p0_min = reg2_min_sec2_inb_tight[0]; p1_min = reg2_min_sec2_inb_tight[1]; p2_min = reg2_min_sec2_inb_tight[2]; p3_min = reg2_min_sec2_inb_tight[3]; p4_min = reg2_min_sec2_inb_tight[4];
                p0_max = reg2_max_sec2_inb_tight[0]; p1_max = reg2_max_sec2_inb_tight[1]; p2_max = reg2_max_sec2_inb_tight[2]; p3_max = reg2_max_sec2_inb_tight[3]; p4_max = reg2_max_sec2_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 3 && inbending == true){
                p0_min = reg2_min_sec3_inb_tight[0]; p1_min = reg2_min_sec3_inb_tight[1]; p2_min = reg2_min_sec3_inb_tight[2]; p3_min = reg2_min_sec3_inb_tight[3]; p4_min = reg2_min_sec3_inb_tight[4];
                p0_max = reg2_max_sec3_inb_tight[0]; p1_max = reg2_max_sec3_inb_tight[1]; p2_max = reg2_max_sec3_inb_tight[2]; p3_max = reg2_max_sec3_inb_tight[3]; p4_max = reg2_max_sec3_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 4 && inbending == true){
                p0_min = reg2_min_sec4_inb_tight[0]; p1_min = reg2_min_sec4_inb_tight[1]; p2_min = reg2_min_sec4_inb_tight[2]; p3_min = reg2_min_sec4_inb_tight[3]; p4_min = reg2_min_sec4_inb_tight[4];
                p0_max = reg2_max_sec4_inb_tight[0]; p1_max = reg2_max_sec4_inb_tight[1]; p2_max = reg2_max_sec4_inb_tight[2]; p3_max = reg2_max_sec4_inb_tight[3]; p4_max = reg2_max_sec4_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 5 && inbending == true){
                p0_min = reg2_min_sec5_inb_tight[0]; p1_min = reg2_min_sec5_inb_tight[1]; p2_min = reg2_min_sec5_inb_tight[2]; p3_min = reg2_min_sec5_inb_tight[3]; p4_min = reg2_min_sec5_inb_tight[4];
                p0_max = reg2_max_sec5_inb_tight[0]; p1_max = reg2_max_sec5_inb_tight[1]; p2_max = reg2_max_sec5_inb_tight[2]; p3_max = reg2_max_sec5_inb_tight[3]; p4_max = reg2_max_sec5_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 6 && inbending == true){
                p0_min = reg2_min_sec6_inb_tight[0]; p1_min = reg2_min_sec6_inb_tight[1]; p2_min = reg2_min_sec6_inb_tight[2]; p3_min = reg2_min_sec6_inb_tight[3]; p4_min = reg2_min_sec6_inb_tight[4];
                p0_max = reg2_max_sec6_inb_tight[0]; p1_max = reg2_max_sec6_inb_tight[1]; p2_max = reg2_max_sec6_inb_tight[2]; p3_max = reg2_max_sec6_inb_tight[3]; p4_max = reg2_max_sec6_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 1 && outbending == true){
                p0_min = reg2_min_sec1_outb_tight[0]; p1_min = reg2_min_sec1_outb_tight[1]; p2_min = reg2_min_sec1_outb_tight[2]; p3_min = reg2_min_sec1_outb_tight[3]; p4_min = reg2_min_sec1_outb_tight[4];
                p0_max = reg2_max_sec1_outb_tight[0]; p1_max = reg2_max_sec1_outb_tight[1]; p2_max = reg2_max_sec1_outb_tight[2]; p3_max = reg2_max_sec1_outb_tight[3]; p4_max = reg2_max_sec1_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 2 && outbending == true){
                p0_min = reg2_min_sec2_outb_tight[0]; p1_min = reg2_min_sec2_outb_tight[1]; p2_min = reg2_min_sec2_outb_tight[2]; p3_min = reg2_min_sec2_outb_tight[3]; p4_min = reg2_min_sec2_outb_tight[4];
                p0_max = reg2_max_sec2_outb_tight[0]; p1_max = reg2_max_sec2_outb_tight[1]; p2_max = reg2_max_sec2_outb_tight[2]; p3_max = reg2_max_sec2_outb_tight[3]; p4_max = reg2_max_sec2_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 3 && outbending == true){
                p0_min = reg2_min_sec3_outb_tight[0]; p1_min = reg2_min_sec3_outb_tight[1]; p2_min = reg2_min_sec3_outb_tight[2]; p3_min = reg2_min_sec3_outb_tight[3]; p4_min = reg2_min_sec3_outb_tight[4];
                p0_max = reg2_max_sec3_outb_tight[0]; p1_max = reg2_max_sec3_outb_tight[1]; p2_max = reg2_max_sec3_outb_tight[2]; p3_max = reg2_max_sec3_outb_tight[3]; p4_max = reg2_max_sec3_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 4 && outbending == true){
                p0_min = reg2_min_sec4_outb_tight[0]; p1_min = reg2_min_sec4_outb_tight[1]; p2_min = reg2_min_sec4_outb_tight[2]; p3_min = reg2_min_sec4_outb_tight[3]; p4_min = reg2_min_sec4_outb_tight[4];
                p0_max = reg2_max_sec4_outb_tight[0]; p1_max = reg2_max_sec4_outb_tight[1]; p2_max = reg2_max_sec4_outb_tight[2]; p3_max = reg2_max_sec4_outb_tight[3]; p4_max = reg2_max_sec4_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 5 && outbending == true){
                p0_min = reg2_min_sec5_outb_tight[0]; p1_min = reg2_min_sec5_outb_tight[1]; p2_min = reg2_min_sec5_outb_tight[2]; p3_min = reg2_min_sec5_outb_tight[3]; p4_min = reg2_min_sec5_outb_tight[4];
                p0_max = reg2_max_sec5_outb_tight[0]; p1_max = reg2_max_sec5_outb_tight[1]; p2_max = reg2_max_sec5_outb_tight[2]; p3_max = reg2_max_sec5_outb_tight[3]; p4_max = reg2_max_sec5_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 6 && outbending == true){
                p0_min = reg2_min_sec6_outb_tight[0]; p1_min = reg2_min_sec6_outb_tight[1]; p2_min = reg2_min_sec6_outb_tight[2]; p3_min = reg2_min_sec6_outb_tight[3]; p4_min = reg2_min_sec6_outb_tight[4];
                p0_max = reg2_max_sec6_outb_tight[0]; p1_max = reg2_max_sec6_outb_tight[1]; p2_max = reg2_max_sec6_outb_tight[2]; p3_max = reg2_max_sec6_outb_tight[3]; p4_max = reg2_max_sec6_outb_tight[4];
            }
        }
        
        if(region == 3){
            if(part->trk(DC)->getSector() == 1 && inbending == true){
                p0_min = reg3_min_sec1_inb_tight[0]; p1_min = reg3_min_sec1_inb_tight[1]; p2_min = reg3_min_sec1_inb_tight[2]; p3_min = reg3_min_sec1_inb_tight[3]; p4_min = reg3_min_sec1_inb_tight[4];
                p0_max = reg3_max_sec1_inb_tight[0]; p1_max = reg3_max_sec1_inb_tight[1]; p2_max = reg3_max_sec1_inb_tight[2]; p3_max = reg3_max_sec1_inb_tight[3]; p4_max = reg3_max_sec1_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 2 && inbending == true){
                p0_min = reg3_min_sec2_inb_tight[0]; p1_min = reg3_min_sec2_inb_tight[1]; p2_min = reg3_min_sec2_inb_tight[2]; p3_min = reg3_min_sec2_inb_tight[3]; p4_min = reg3_min_sec2_inb_tight[4];
                p0_max = reg3_max_sec2_inb_tight[0]; p1_max = reg3_max_sec2_inb_tight[1]; p2_max = reg3_max_sec2_inb_tight[2]; p3_max = reg3_max_sec2_inb_tight[3]; p4_max = reg3_max_sec2_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 3 && inbending == true){
                p0_min = reg3_min_sec3_inb_tight[0]; p1_min = reg3_min_sec3_inb_tight[1]; p2_min = reg3_min_sec3_inb_tight[2]; p3_min = reg3_min_sec3_inb_tight[3]; p4_min = reg3_min_sec3_inb_tight[4];
                p0_max = reg3_max_sec3_inb_tight[0]; p1_max = reg3_max_sec3_inb_tight[1]; p2_max = reg3_max_sec3_inb_tight[2]; p3_max = reg3_max_sec3_inb_tight[3]; p4_max = reg3_max_sec3_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 4 && inbending == true){
                p0_min = reg3_min_sec4_inb_tight[0]; p1_min = reg3_min_sec4_inb_tight[1]; p2_min = reg3_min_sec4_inb_tight[2]; p3_min = reg3_min_sec4_inb_tight[3]; p4_min = reg3_min_sec4_inb_tight[4];
                p0_max = reg3_max_sec4_inb_tight[0]; p1_max = reg3_max_sec4_inb_tight[1]; p2_max = reg3_max_sec4_inb_tight[2]; p3_max = reg3_max_sec4_inb_tight[3]; p4_max = reg3_max_sec4_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 5 && inbending == true){
                p0_min = reg3_min_sec5_inb_tight[0]; p1_min = reg3_min_sec5_inb_tight[1]; p2_min = reg3_min_sec5_inb_tight[2]; p3_min = reg3_min_sec5_inb_tight[3]; p4_min = reg3_min_sec5_inb_tight[4];
                p0_max = reg3_max_sec5_inb_tight[0]; p1_max = reg3_max_sec5_inb_tight[1]; p2_max = reg3_max_sec5_inb_tight[2]; p3_max = reg3_max_sec5_inb_tight[3]; p4_max = reg3_max_sec5_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 6 && inbending == true){
                p0_min = reg3_min_sec6_inb_tight[0]; p1_min = reg3_min_sec6_inb_tight[1]; p2_min = reg3_min_sec6_inb_tight[2]; p3_min = reg3_min_sec6_inb_tight[3]; p4_min = reg3_min_sec6_inb_tight[4];
                p0_max = reg3_max_sec6_inb_tight[0]; p1_max = reg3_max_sec6_inb_tight[1]; p2_max = reg3_max_sec6_inb_tight[2]; p3_max = reg3_max_sec6_inb_tight[3]; p4_max = reg3_max_sec6_inb_tight[4];
            }
            if(part->trk(DC)->getSector() == 1 && outbending == true){
                p0_min = reg3_min_sec1_outb_tight[0]; p1_min = reg3_min_sec1_outb_tight[1]; p2_min = reg3_min_sec1_outb_tight[2]; p3_min = reg3_min_sec1_outb_tight[3]; p4_min = reg3_min_sec1_outb_tight[4];
                p0_max = reg3_max_sec1_outb_tight[0]; p1_max = reg3_max_sec1_outb_tight[1]; p2_max = reg3_max_sec1_outb_tight[2]; p3_max = reg3_max_sec1_outb_tight[3]; p4_max = reg3_max_sec1_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 2 && outbending == true){
                p0_min = reg3_min_sec2_outb_tight[0]; p1_min = reg3_min_sec2_outb_tight[1]; p2_min = reg3_min_sec2_outb_tight[2]; p3_min = reg3_min_sec2_outb_tight[3]; p4_min = reg3_min_sec2_outb_tight[4];
                p0_max = reg3_max_sec2_outb_tight[0]; p1_max = reg3_max_sec2_outb_tight[1]; p2_max = reg3_max_sec2_outb_tight[2]; p3_max = reg3_max_sec2_outb_tight[3]; p4_max = reg3_max_sec2_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 3 && outbending == true){
                p0_min = reg3_min_sec3_outb_tight[0]; p1_min = reg3_min_sec3_outb_tight[1]; p2_min = reg3_min_sec3_outb_tight[2]; p3_min = reg3_min_sec3_outb_tight[3]; p4_min = reg3_min_sec3_outb_tight[4];
                p0_max = reg3_max_sec3_outb_tight[0]; p1_max = reg3_max_sec3_outb_tight[1]; p2_max = reg3_max_sec3_outb_tight[2]; p3_max = reg3_max_sec3_outb_tight[3]; p4_max = reg3_max_sec3_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 4 && outbending == true){
                p0_min = reg3_min_sec4_outb_tight[0]; p1_min = reg3_min_sec4_outb_tight[1]; p2_min = reg3_min_sec4_outb_tight[2]; p3_min = reg3_min_sec4_outb_tight[3]; p4_min = reg3_min_sec4_outb_tight[4];
                p0_max = reg3_max_sec4_outb_tight[0]; p1_max = reg3_max_sec4_outb_tight[1]; p2_max = reg3_max_sec4_outb_tight[2]; p3_max = reg3_max_sec4_outb_tight[3]; p4_max = reg3_max_sec4_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 5 && outbending == true){
                p0_min = reg3_min_sec5_outb_tight[0]; p1_min = reg3_min_sec5_outb_tight[1]; p2_min = reg3_min_sec5_outb_tight[2]; p3_min = reg3_min_sec5_outb_tight[3]; p4_min = reg3_min_sec5_outb_tight[4];
                p0_max = reg3_max_sec5_outb_tight[0]; p1_max = reg3_max_sec5_outb_tight[1]; p2_max = reg3_max_sec5_outb_tight[2]; p3_max = reg3_max_sec5_outb_tight[3]; p4_max = reg3_max_sec5_outb_tight[4];
            }
            if(part->trk(DC)->getSector() == 6 && outbending == true){
                p0_min = reg3_min_sec6_outb_tight[0]; p1_min = reg3_min_sec6_outb_tight[1]; p2_min = reg3_min_sec6_outb_tight[2]; p3_min = reg3_min_sec6_outb_tight[3]; p4_min = reg3_min_sec6_outb_tight[4];
                p0_max = reg3_max_sec6_outb_tight[0]; p1_max = reg3_max_sec6_outb_tight[1]; p2_max = reg3_max_sec6_outb_tight[2]; p3_max = reg3_max_sec6_outb_tight[3]; p4_max = reg3_max_sec6_outb_tight[4];
            }
        }
    }
    
    if(loose == true){
        if(region == 1){
            if(part->trk(DC)->getSector() == 1 && inbending == true){
                p0_min = reg1_min_sec1_inb_loose[0]; p1_min = reg1_min_sec1_inb_loose[1]; p2_min = reg1_min_sec1_inb_loose[2]; p3_min = reg1_min_sec1_inb_loose[3]; p4_min = reg1_min_sec1_inb_loose[4];
                p0_max = reg1_max_sec1_inb_loose[0]; p1_max = reg1_max_sec1_inb_loose[1]; p2_max = reg1_max_sec1_inb_loose[2]; p3_max = reg1_max_sec1_inb_loose[3]; p4_max = reg1_max_sec1_inb_loose[4];
            }
            if(part->trk(DC)->getSector()== 2 && inbending == true){
                p0_min = reg1_min_sec2_inb_loose[0]; p1_min = reg1_min_sec2_inb_loose[1]; p2_min = reg1_min_sec2_inb_loose[2]; p3_min = reg1_min_sec2_inb_loose[3]; p4_min = reg1_min_sec2_inb_loose[4];
                p0_max = reg1_max_sec2_inb_loose[0]; p1_max = reg1_max_sec2_inb_loose[1]; p2_max = reg1_max_sec2_inb_loose[2]; p3_max = reg1_max_sec2_inb_loose[3]; p4_max = reg1_max_sec2_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 3 && inbending == true){
                p0_min = reg1_min_sec3_inb_loose[0]; p1_min = reg1_min_sec3_inb_loose[1]; p2_min = reg1_min_sec3_inb_loose[2]; p3_min = reg1_min_sec3_inb_loose[3]; p4_min = reg1_min_sec3_inb_loose[4];
                p0_max = reg1_max_sec3_inb_loose[0]; p1_max = reg1_max_sec3_inb_loose[1]; p2_max = reg1_max_sec3_inb_loose[2]; p3_max = reg1_max_sec3_inb_loose[3]; p4_max = reg1_max_sec3_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 4 && inbending == true){
                p0_min = reg1_min_sec4_inb_loose[0]; p1_min = reg1_min_sec4_inb_loose[1]; p2_min = reg1_min_sec4_inb_loose[2]; p3_min = reg1_min_sec4_inb_loose[3]; p4_min = reg1_min_sec4_inb_loose[4];
                p0_max = reg1_max_sec4_inb_loose[0]; p1_max = reg1_max_sec4_inb_loose[1]; p2_max = reg1_max_sec4_inb_loose[2]; p3_max = reg1_max_sec4_inb_loose[3]; p4_max = reg1_max_sec4_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 5 && inbending == true){
                p0_min = reg1_min_sec5_inb_loose[0]; p1_min = reg1_min_sec5_inb_loose[1]; p2_min = reg1_min_sec5_inb_loose[2]; p3_min = reg1_min_sec5_inb_loose[3]; p4_min = reg1_min_sec5_inb_loose[4];
                p0_max = reg1_max_sec5_inb_loose[0]; p1_max = reg1_max_sec5_inb_loose[1]; p2_max = reg1_max_sec5_inb_loose[2]; p3_max = reg1_max_sec5_inb_loose[3]; p4_max = reg1_max_sec5_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 6 && inbending == true){
                p0_min = reg1_min_sec6_inb_loose[0]; p1_min = reg1_min_sec6_inb_loose[1]; p2_min = reg1_min_sec6_inb_loose[2]; p3_min = reg1_min_sec6_inb_loose[3]; p4_min = reg1_min_sec6_inb_loose[4];
                p0_max = reg1_max_sec6_inb_loose[0]; p1_max = reg1_max_sec6_inb_loose[1]; p2_max = reg1_max_sec6_inb_loose[2]; p3_max = reg1_max_sec6_inb_loose[3]; p4_max = reg1_max_sec6_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 1 && outbending == true){
                p0_min = reg1_min_sec1_outb_loose[0]; p1_min = reg1_min_sec1_outb_loose[1]; p2_min = reg1_min_sec1_outb_loose[2]; p3_min = reg1_min_sec1_outb_loose[3]; p4_min = reg1_min_sec1_outb_loose[4];
                p0_max = reg1_max_sec1_outb_loose[0]; p1_max = reg1_max_sec1_outb_loose[1]; p2_max = reg1_max_sec1_outb_loose[2]; p3_max = reg1_max_sec1_outb_loose[3]; p4_max = reg1_max_sec1_outb_loose[4];
            }
            if(part->trk(DC)->getSector() ==  2 && outbending == true){
                p0_min = reg1_min_sec2_outb_loose[0]; p1_min = reg1_min_sec2_outb_loose[1]; p2_min = reg1_min_sec2_outb_loose[2]; p3_min = reg1_min_sec2_outb_loose[3]; p4_min = reg1_min_sec2_outb_loose[4];
                p0_max = reg1_max_sec2_outb_loose[0]; p1_max = reg1_max_sec2_outb_loose[1]; p2_max = reg1_max_sec2_outb_loose[2]; p3_max = reg1_max_sec2_outb_loose[3]; p4_max = reg1_max_sec2_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 3 && outbending == true){
                p0_min = reg1_min_sec3_outb_loose[0]; p1_min = reg1_min_sec3_outb_loose[1]; p2_min = reg1_min_sec3_outb_loose[2]; p3_min = reg1_min_sec3_outb_loose[3]; p4_min = reg1_min_sec3_outb_loose[4];
                p0_max = reg1_max_sec3_outb_loose[0]; p1_max = reg1_max_sec3_outb_loose[1]; p2_max = reg1_max_sec3_outb_loose[2]; p3_max = reg1_max_sec3_outb_loose[3]; p4_max = reg1_max_sec3_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 4 && outbending == true){
                p0_min = reg1_min_sec4_outb_loose[0]; p1_min = reg1_min_sec4_outb_loose[1]; p2_min = reg1_min_sec4_outb_loose[2]; p3_min = reg1_min_sec4_outb_loose[3]; p4_min = reg1_min_sec4_outb_loose[4];
                p0_max = reg1_max_sec4_outb_loose[0]; p1_max = reg1_max_sec4_outb_loose[1]; p2_max = reg1_max_sec4_outb_loose[2]; p3_max = reg1_max_sec4_outb_loose[3]; p4_max = reg1_max_sec4_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 5 && outbending == true){
                p0_min = reg1_min_sec5_outb_loose[0]; p1_min = reg1_min_sec5_outb_loose[1]; p2_min = reg1_min_sec5_outb_loose[2]; p3_min = reg1_min_sec5_outb_loose[3]; p4_min = reg1_min_sec5_outb_loose[4];
                p0_max = reg1_max_sec5_outb_loose[0]; p1_max = reg1_max_sec5_outb_loose[1]; p2_max = reg1_max_sec5_outb_loose[2]; p3_max = reg1_max_sec5_outb_loose[3]; p4_max = reg1_max_sec5_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 6 && outbending == true){
                p0_min = reg1_min_sec6_outb_loose[0]; p1_min = reg1_min_sec6_outb_loose[1]; p2_min = reg1_min_sec6_outb_loose[2]; p3_min = reg1_min_sec6_outb_loose[3]; p4_min = reg1_min_sec6_outb_loose[4];
                p0_max = reg1_max_sec6_outb_loose[0]; p1_max = reg1_max_sec6_outb_loose[1]; p2_max = reg1_max_sec6_outb_loose[2]; p3_max = reg1_max_sec6_outb_loose[3]; p4_max = reg1_max_sec6_outb_loose[4];
            }
        }
        
        if(region == 2){
            if(part->trk(DC)->getSector() == 1 && inbending == true){
                p0_min = reg2_min_sec1_inb_loose[0]; p1_min = reg2_min_sec1_inb_loose[1]; p2_min = reg2_min_sec1_inb_loose[2]; p3_min = reg2_min_sec1_inb_loose[3]; p4_min = reg2_min_sec1_inb_loose[4];
                p0_max = reg2_max_sec1_inb_loose[0]; p1_max = reg2_max_sec1_inb_loose[1]; p2_max = reg2_max_sec1_inb_loose[2]; p3_max = reg2_max_sec1_inb_loose[3]; p4_max = reg2_max_sec1_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 2 && inbending == true){
                p0_min = reg2_min_sec2_inb_loose[0]; p1_min = reg2_min_sec2_inb_loose[1]; p2_min = reg2_min_sec2_inb_loose[2]; p3_min = reg2_min_sec2_inb_loose[3]; p4_min = reg2_min_sec2_inb_loose[4];
                p0_max = reg2_max_sec2_inb_loose[0]; p1_max = reg2_max_sec2_inb_loose[1]; p2_max = reg2_max_sec2_inb_loose[2]; p3_max = reg2_max_sec2_inb_loose[3]; p4_max = reg2_max_sec2_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 3 && inbending == true){
                p0_min = reg2_min_sec3_inb_loose[0]; p1_min = reg2_min_sec3_inb_loose[1]; p2_min = reg2_min_sec3_inb_loose[2]; p3_min = reg2_min_sec3_inb_loose[3]; p4_min = reg2_min_sec3_inb_loose[4];
                p0_max = reg2_max_sec3_inb_loose[0]; p1_max = reg2_max_sec3_inb_loose[1]; p2_max = reg2_max_sec3_inb_loose[2]; p3_max = reg2_max_sec3_inb_loose[3]; p4_max = reg2_max_sec3_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 4 && inbending == true){
                p0_min = reg2_min_sec4_inb_loose[0]; p1_min = reg2_min_sec4_inb_loose[1]; p2_min = reg2_min_sec4_inb_loose[2]; p3_min = reg2_min_sec4_inb_loose[3]; p4_min = reg2_min_sec4_inb_loose[4];
                p0_max = reg2_max_sec4_inb_loose[0]; p1_max = reg2_max_sec4_inb_loose[1]; p2_max = reg2_max_sec4_inb_loose[2]; p3_max = reg2_max_sec4_inb_loose[3]; p4_max = reg2_max_sec4_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 5 && inbending == true){
                p0_min = reg2_min_sec5_inb_loose[0]; p1_min = reg2_min_sec5_inb_loose[1]; p2_min = reg2_min_sec5_inb_loose[2]; p3_min = reg2_min_sec5_inb_loose[3]; p4_min = reg2_min_sec5_inb_loose[4];
                p0_max = reg2_max_sec5_inb_loose[0]; p1_max = reg2_max_sec5_inb_loose[1]; p2_max = reg2_max_sec5_inb_loose[2]; p3_max = reg2_max_sec5_inb_loose[3]; p4_max = reg2_max_sec5_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 6 && inbending == true){
                p0_min = reg2_min_sec6_inb_loose[0]; p1_min = reg2_min_sec6_inb_loose[1]; p2_min = reg2_min_sec6_inb_loose[2]; p3_min = reg2_min_sec6_inb_loose[3]; p4_min = reg2_min_sec6_inb_loose[4];
                p0_max = reg2_max_sec6_inb_loose[0]; p1_max = reg2_max_sec6_inb_loose[1]; p2_max = reg2_max_sec6_inb_loose[2]; p3_max = reg2_max_sec6_inb_loose[3]; p4_max = reg2_max_sec6_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 1 && outbending == true){
                p0_min = reg2_min_sec1_outb_loose[0]; p1_min = reg2_min_sec1_outb_loose[1]; p2_min = reg2_min_sec1_outb_loose[2]; p3_min = reg2_min_sec1_outb_loose[3]; p4_min = reg2_min_sec1_outb_loose[4];
                p0_max = reg2_max_sec1_outb_loose[0]; p1_max = reg2_max_sec1_outb_loose[1]; p2_max = reg2_max_sec1_outb_loose[2]; p3_max = reg2_max_sec1_outb_loose[3]; p4_max = reg2_max_sec1_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 2 && outbending == true){
                p0_min = reg2_min_sec2_outb_loose[0]; p1_min = reg2_min_sec2_outb_loose[1]; p2_min = reg2_min_sec2_outb_loose[2]; p3_min = reg2_min_sec2_outb_loose[3]; p4_min = reg2_min_sec2_outb_loose[4];
                p0_max = reg2_max_sec2_outb_loose[0]; p1_max = reg2_max_sec2_outb_loose[1]; p2_max = reg2_max_sec2_outb_loose[2]; p3_max = reg2_max_sec2_outb_loose[3]; p4_max = reg2_max_sec2_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 3 && outbending == true){
                p0_min = reg2_min_sec3_outb_loose[0]; p1_min = reg2_min_sec3_outb_loose[1]; p2_min = reg2_min_sec3_outb_loose[2]; p3_min = reg2_min_sec3_outb_loose[3]; p4_min = reg2_min_sec3_outb_loose[4];
                p0_max = reg2_max_sec3_outb_loose[0]; p1_max = reg2_max_sec3_outb_loose[1]; p2_max = reg2_max_sec3_outb_loose[2]; p3_max = reg2_max_sec3_outb_loose[3]; p4_max = reg2_max_sec3_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 4 && outbending == true){
                p0_min = reg2_min_sec4_outb_loose[0]; p1_min = reg2_min_sec4_outb_loose[1]; p2_min = reg2_min_sec4_outb_loose[2]; p3_min = reg2_min_sec4_outb_loose[3]; p4_min = reg2_min_sec4_outb_loose[4];
                p0_max = reg2_max_sec4_outb_loose[0]; p1_max = reg2_max_sec4_outb_loose[1]; p2_max = reg2_max_sec4_outb_loose[2]; p3_max = reg2_max_sec4_outb_loose[3]; p4_max = reg2_max_sec4_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 5 && outbending == true){
                p0_min = reg2_min_sec5_outb_loose[0]; p1_min = reg2_min_sec5_outb_loose[1]; p2_min = reg2_min_sec5_outb_loose[2]; p3_min = reg2_min_sec5_outb_loose[3]; p4_min = reg2_min_sec5_outb_loose[4];
                p0_max = reg2_max_sec5_outb_loose[0]; p1_max = reg2_max_sec5_outb_loose[1]; p2_max = reg2_max_sec5_outb_loose[2]; p3_max = reg2_max_sec5_outb_loose[3]; p4_max = reg2_max_sec5_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 6 && outbending == true){
                p0_min = reg2_min_sec6_outb_loose[0]; p1_min = reg2_min_sec6_outb_loose[1]; p2_min = reg2_min_sec6_outb_loose[2]; p3_min = reg2_min_sec6_outb_loose[3]; p4_min = reg2_min_sec6_outb_loose[4];
                p0_max = reg2_max_sec6_outb_loose[0]; p1_max = reg2_max_sec6_outb_loose[1]; p2_max = reg2_max_sec6_outb_loose[2]; p3_max = reg2_max_sec6_outb_loose[3]; p4_max = reg2_max_sec6_outb_loose[4];
            }
        }
        
        if(region == 3){
            if(part->trk(DC)->getSector() == 1 && inbending == true){
                p0_min = reg3_min_sec1_inb_loose[0]; p1_min = reg3_min_sec1_inb_loose[1]; p2_min = reg3_min_sec1_inb_loose[2]; p3_min = reg3_min_sec1_inb_loose[3]; p4_min = reg3_min_sec1_inb_loose[4];
                p0_max = reg3_max_sec1_inb_loose[0]; p1_max = reg3_max_sec1_inb_loose[1]; p2_max = reg3_max_sec1_inb_loose[2]; p3_max = reg3_max_sec1_inb_loose[3]; p4_max = reg3_max_sec1_inb_loose[4];
            }
            if(part->trk(DC)->getSector() ==  2 && inbending == true){
                p0_min = reg3_min_sec2_inb_loose[0]; p1_min = reg3_min_sec2_inb_loose[1]; p2_min = reg3_min_sec2_inb_loose[2]; p3_min = reg3_min_sec2_inb_loose[3]; p4_min = reg3_min_sec2_inb_loose[4];
                p0_max = reg3_max_sec2_inb_loose[0]; p1_max = reg3_max_sec2_inb_loose[1]; p2_max = reg3_max_sec2_inb_loose[2]; p3_max = reg3_max_sec2_inb_loose[3]; p4_max = reg3_max_sec2_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 3 && inbending == true){
                p0_min = reg3_min_sec3_inb_loose[0]; p1_min = reg3_min_sec3_inb_loose[1]; p2_min = reg3_min_sec3_inb_loose[2]; p3_min = reg3_min_sec3_inb_loose[3]; p4_min = reg3_min_sec3_inb_loose[4];
                p0_max = reg3_max_sec3_inb_loose[0]; p1_max = reg3_max_sec3_inb_loose[1]; p2_max = reg3_max_sec3_inb_loose[2]; p3_max = reg3_max_sec3_inb_loose[3]; p4_max = reg3_max_sec3_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 4 && inbending == true){
                p0_min = reg3_min_sec4_inb_loose[0]; p1_min = reg3_min_sec4_inb_loose[1]; p2_min = reg3_min_sec4_inb_loose[2]; p3_min = reg3_min_sec4_inb_loose[3]; p4_min = reg3_min_sec4_inb_loose[4];
                p0_max = reg3_max_sec4_inb_loose[0]; p1_max = reg3_max_sec4_inb_loose[1]; p2_max = reg3_max_sec4_inb_loose[2]; p3_max = reg3_max_sec4_inb_loose[3]; p4_max = reg3_max_sec4_inb_loose[4];
            }
            if(part->trk(DC)->getSector() ==  5 && inbending == true){
                p0_min = reg3_min_sec5_inb_loose[0]; p1_min = reg3_min_sec5_inb_loose[1]; p2_min = reg3_min_sec5_inb_loose[2]; p3_min = reg3_min_sec5_inb_loose[3]; p4_min = reg3_min_sec5_inb_loose[4];
                p0_max = reg3_max_sec5_inb_loose[0]; p1_max = reg3_max_sec5_inb_loose[1]; p2_max = reg3_max_sec5_inb_loose[2]; p3_max = reg3_max_sec5_inb_loose[3]; p4_max = reg3_max_sec5_inb_loose[4];
            }
            if(part->trk(DC)->getSector() ==  6 && inbending == true){
                p0_min = reg3_min_sec6_inb_loose[0]; p1_min = reg3_min_sec6_inb_loose[1]; p2_min = reg3_min_sec6_inb_loose[2]; p3_min = reg3_min_sec6_inb_loose[3]; p4_min = reg3_min_sec6_inb_loose[4];
                p0_max = reg3_max_sec6_inb_loose[0]; p1_max = reg3_max_sec6_inb_loose[1]; p2_max = reg3_max_sec6_inb_loose[2]; p3_max = reg3_max_sec6_inb_loose[3]; p4_max = reg3_max_sec6_inb_loose[4];
            }
            if(part->trk(DC)->getSector() == 1 && outbending == true){
                p0_min = reg3_min_sec1_outb_loose[0]; p1_min = reg3_min_sec1_outb_loose[1]; p2_min = reg3_min_sec1_outb_loose[2]; p3_min = reg3_min_sec1_outb_loose[3]; p4_min = reg3_min_sec1_outb_loose[4];
                p0_max = reg3_max_sec1_outb_loose[0]; p1_max = reg3_max_sec1_outb_loose[1]; p2_max = reg3_max_sec1_outb_loose[2]; p3_max = reg3_max_sec1_outb_loose[3]; p4_max = reg3_max_sec1_outb_loose[4];
            }
            if(part->trk(DC)->getSector() ==  2 && outbending == true){
                p0_min = reg3_min_sec2_outb_loose[0]; p1_min = reg3_min_sec2_outb_loose[1]; p2_min = reg3_min_sec2_outb_loose[2]; p3_min = reg3_min_sec2_outb_loose[3]; p4_min = reg3_min_sec2_outb_loose[4];
                p0_max = reg3_max_sec2_outb_loose[0]; p1_max = reg3_max_sec2_outb_loose[1]; p2_max = reg3_max_sec2_outb_loose[2]; p3_max = reg3_max_sec2_outb_loose[3]; p4_max = reg3_max_sec2_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 3 && outbending == true){
                p0_min = reg3_min_sec3_outb_loose[0]; p1_min = reg3_min_sec3_outb_loose[1]; p2_min = reg3_min_sec3_outb_loose[2]; p3_min = reg3_min_sec3_outb_loose[3]; p4_min = reg3_min_sec3_outb_loose[4];
                p0_max = reg3_max_sec3_outb_loose[0]; p1_max = reg3_max_sec3_outb_loose[1]; p2_max = reg3_max_sec3_outb_loose[2]; p3_max = reg3_max_sec3_outb_loose[3]; p4_max = reg3_max_sec3_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 4 && outbending == true){
                p0_min = reg3_min_sec4_outb_loose[0]; p1_min = reg3_min_sec4_outb_loose[1]; p2_min = reg3_min_sec4_outb_loose[2]; p3_min = reg3_min_sec4_outb_loose[3]; p4_min = reg3_min_sec4_outb_loose[4];
                p0_max = reg3_max_sec4_outb_loose[0]; p1_max = reg3_max_sec4_outb_loose[1]; p2_max = reg3_max_sec4_outb_loose[2]; p3_max = reg3_max_sec4_outb_loose[3]; p4_max = reg3_max_sec4_outb_loose[4];
            }
            if(part->trk(DC)->getSector() ==  5 && outbending == true){
                p0_min = reg3_min_sec5_outb_loose[0]; p1_min = reg3_min_sec5_outb_loose[1]; p2_min = reg3_min_sec5_outb_loose[2]; p3_min = reg3_min_sec5_outb_loose[3]; p4_min = reg3_min_sec5_outb_loose[4];
                p0_max = reg3_max_sec5_outb_loose[0]; p1_max = reg3_max_sec5_outb_loose[1]; p2_max = reg3_max_sec5_outb_loose[2]; p3_max = reg3_max_sec5_outb_loose[3]; p4_max = reg3_max_sec5_outb_loose[4];
            }
            if(part->trk(DC)->getSector() == 6 && outbending == true){
                p0_min = reg3_min_sec6_outb_loose[0]; p1_min = reg3_min_sec6_outb_loose[1]; p2_min = reg3_min_sec6_outb_loose[2]; p3_min = reg3_min_sec6_outb_loose[3]; p4_min = reg3_min_sec6_outb_loose[4];
                p0_max = reg3_max_sec6_outb_loose[0]; p1_max = reg3_max_sec6_outb_loose[1]; p2_max = reg3_max_sec6_outb_loose[2]; p3_max = reg3_max_sec6_outb_loose[3]; p4_max = reg3_max_sec6_outb_loose[4];
            }
        }
    }
    
    double phi_DC_min = p0_min + p1_min * log(theta_DC) + p2_min * theta_DC + p3_min * theta_DC*theta_DC + p4_min * theta_DC*theta_DC*theta_DC;
    double phi_DC_max = p0_max + p1_max * log(theta_DC) + p2_max * theta_DC + p3_max * theta_DC*theta_DC + p4_max * theta_DC*theta_DC*theta_DC;
    
    if(phi_DC_min < -25.5) phi_DC_min = -25.5;
    if(phi_DC_max > +25.5) phi_DC_max = +25.5;
    
    if(phi_DC > phi_DC_min && phi_DC < phi_DC_max) return true;
    else return false;
}

/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// 2.2 DC fiducial cuts for the 3 regions based on the chi2/NDF value of the tracks (for elctrons and hadrons)
///     - The cut checks at which point the chi2/NDF value starts to increase
///     - works well for hadrons but for electrons no clear increase can be observed in many bins, which leads to soem uncertainty


bool DC_fiducial_cut_chi2(region_part_ptr part,int region){
    
    double Pival=TMath::Pi();
    
    
    double  part_DC_c1x=part->traj(DC,6)->getX();
    double  part_DC_c1y=part->traj(DC,6)->getY();
    double  part_DC_c1z=part->traj(DC,6)->getZ();
    double  part_DC_c2x=part->traj(DC,18)->getX();
    double  part_DC_c2y=part->traj(DC,18)->getY();
    double  part_DC_c2z=part->traj(DC,18)->getZ();
    double  part_DC_c3x=part->traj(DC,36)->getX();
    double  part_DC_c3y=part->traj(DC,36)->getY();
    double  part_DC_c3z=part->traj(DC,36)->getY();
    
    
    //fitted values
    const double maxparams[6][6][3][4] = {
        {{{-35.1716, 25.102, -0.750281, 5.34679e-05},
            {-39.1633, 28.5551, -1.13429, 0.00419047},
            {-33.7705, 24.8068, -0.811239, 0.00138345}},
            {{-36.2389, 26.7979, -1.08147, 0.0050898},
                {-43.643, 31.6783, -1.49203, 0.00872922},
                {-54.4042, 40.6516, -2.52393, 0.0205649}},
            {{-38.3238, 26.1667, -0.777077, 0.000264835},
                {-34.2011, 24.2843, -0.696392, 3.75866e-12},
                {-36.4636, 25.8712, -0.786592, 2.24421e-10}},
            {{-31.8019, 23.154, -0.653992, 2.69968e-05},
                {-34.6637, 24.6043, -0.714901, 2.02675e-10},
                {-36.7209, 26.2469, -0.828638, 0.000340435}},
            {{-33.4016, 24.6901, -0.779889, 0.000430557},
                {-35.4583, 24.7491, -0.707953, 2.18559e-10},
                {-37.7335, 28.1547, -1.1986, 0.00582395}},
            {{-34.7808, 24.6988, -0.719936, 5.73299e-10},
                {-54.5797, 40.9138, -2.57493, 0.0213354},
                {-38.4972, 28.3142, -1.21741, 0.00640373}}},
        {{{-2.25358e-08, 12.631, -0.767619, 0.00739811},
            {-8.09501, 15.9098, -0.844083, 0.00667995},
            {-1.48113e-06, 12.2061, -0.73167, 0.0074309}},
            {{-2.10872e-07, 12.6689, -0.765156, 0.00720044},
                {-4.88862, 14.0376, -0.687202, 0.00506307},
                {-4.59793e-06, 11.5553, -0.591469, 0.00536957}},
            {{-1.13504e-08, 12.6011, -0.746025, 0.00687498},
                {-6.97884, 15.1788, -0.765889, 0.00570532},
                {-1.29468, 12.3844, -0.667561, 0.00619226}},
            {{-2.91953e-09, 13.883, -0.999624, 0.0104257},
                {-4.9855, 13.8864, -0.661348, 0.0048371},
                {-2.29438e-08, 11.8341, -0.668486, 0.00669247}},
            {{-2.02824e-08, 13.3855, -0.91158, 0.00926769},
                {-3.29092e-08, 10.8294, -0.382323, 0.00178367},
                {-4.59027e-06, 11.9414, -0.663872, 0.00625769}},
            {{-3.73322e-09, 12.6126, -0.723548, 0.0062217},
                {-4.56248, 14.1574, -0.727805, 0.00560108},
                {-2.39381e-08, 12.0663, -0.6651, 0.00602544}}},
        {{{-1.45923e-08, 13.0297, -0.828302, 0.00795271},
            {-5.41905, 13.2753, -0.503236, 0.00255607},
            {-3.67719, 12.1358, -0.462905, 0.00308219}},
            {{-9.953e-10, 11.549, -0.52816, 0.00378771},
                {-8.47154, 15.9863, -0.826166, 0.0062936},
                {-6.43715, 13.9081, -0.618535, 0.0046102}},
            {{-4.68458e-08, 12.9481, -0.781613, 0.00689754},
                {-3.46617, 12.2786, -0.440121, 0.00205448},
                {-4.43519, 10.9372, -0.210059, 3.69283e-10}},
            {{-4.18414e-07, 13.1542, -0.811251, 0.00714402},
                {-4.63166, 13.7769, -0.657207, 0.0047586},
                {-1.99278e-05, 11.3993, -0.575232, 0.00532141}},
            {{-7.07189e-10, 13.2814, -0.88476, 0.00874389},
                {-5.08373, 14.4384, -0.750795, 0.00586116},
                {-6.9642e-05, 9.50651, -0.189316, 3.07274e-06}},
            {{-5.85515e-08, 12.5116, -0.688741, 0.00557297},
                {-1.86306, 11.985, -0.482567, 0.00279836},
                {-4.94295e-07, 10.1342, -0.316715, 0.00176254}}},
        {{{-0.0157256, 11.1508, -0.415185, 0.00186904},
            {-13.6561, 19.4418, -1.15773, 0.00989432},
            {-6.24969e-07, 10.5776, -0.329325, 0.00103488}},
            {{-2.5686e-08, 11.4797, -0.476772, 0.00264288},
                {-0.0475099, 10.1207, -0.244786, 3.13032e-06},
                {-4.6875e-07, 12.019, -0.63598, 0.00543214}},
            {{-0.00702545, 11.1294, -0.407207, 0.00171263},
                {-7.27687, 15.5, -0.807858, 0.0062086},
                {-5.15078, 12.6368, -0.348584, 9.2687e-12}},
            {{-8.14106e-08, 13.28, -0.818164, 0.00703758},
                {-7.60722, 14.4871, -0.588662, 0.00326244},
                {-1.70764e-06, 12.0413, -0.63961, 0.00541784}},
            {{-1.09281, 11.5573, -0.41311, 0.00155228},
                {-3.71599, 12.8335, -0.521472, 0.00296792},
                {-0.000410815, 12.4833, -0.72999, 0.0066601}},
            {{-0.652641, 12.2766, -0.554202, 0.00314615},
                {-8.42824, 15.5087, -0.710609, 0.00447051},
                {-14.9692, 21.5885, -1.47528, 0.0136615}}},
        {{{-5.58945, 17.4004, -1.34516, 0.0142099},
            {-14.9585, 20.4538, -1.25118, 0.0106617},
            {-12.0069, 16.4545, -0.727162, 0.00495418}},
            {{-7.03048, 17.3519, -1.1831, 0.0111308},
                {-7.30641, 15.8503, -0.850952, 0.00648446},
                {-10.2549, 15.6139, -0.648352, 0.00380506}},
            {{-9.73111e-09, 13.498, -0.932479, 0.00939708},
                {-8.38053, 15.5588, -0.711323, 0.00433827},
                {-12.3097, 16.6403, -0.741362, 0.0050708}},
            {{-7.38905, 17.2652, -1.15517, 0.0109165},
                {-1.11835e-07, 10.4637, -0.301972, 0.000612754},
                {-12.2182, 17.4958, -0.919555, 0.00747512}},
            {{-0.492676, 14.4148, -1.0959, 0.0116708},
                {-5.34309, 14.3258, -0.691954, 0.00480109},
                {-12.5443, 16.1047, -0.59594, 0.00280171}},
            {{-4.08375e-07, 12.2846, -0.655278, 0.00525956},
                {-8.93101, 16.4266, -0.861853, 0.00644623},
                {-11.8406, 17.0417, -0.826301, 0.00596028}}},
        {{{-9.29415, 16.5566, -0.831923, 0.00562504},
            {-0.954483, 10.5813, -0.265766, 3.24615e-05},
            {-6.87423, 14.892, -0.76495, 0.00639603}},
            {{-18.8913, 19.3123, -0.711917, 0.00227889},
                {-13.9788, 18.5678, -0.940183, 0.00664397},
                {-11.7696, 18.3415, -1.04368, 0.0083506}},
            {{-3.82873, 12.7727, -0.425968, 0.000789835},
                {-9.81221, 14.6531, -0.471092, 0.00131406},
                {-14.2392, 15.9895, -0.430525, 2.20712e-12}},
            {{-1.76975e-07, 11.4006, -0.420134, 0.00141302},
                {-3.11764, 10.9707, -0.245823, 2.23044e-12},
                {-17.6005, 22.2881, -1.39992, 0.0117791}},
            {{-0.767518, 11.6824, -0.456275, 0.00214005},
                {-5.28047, 12.65, -0.350658, 9.80081e-05},
                {-0.0888832, 11.508, -0.49197, 0.00301269}},
            {{-4.72388, 15.8507, -1.00574, 0.00876768},
                {-2.80649, 11.4056, -0.301812, 0.000190262},
                {-13.0484, 18.665, -1.08614, 0.00960977}}}};
    const double minparams[6][6][3][4] = {
        {{{37.289, -27.5201, 1.12866, -0.00526111},
            {45.3103, -33.5226, 1.72923, -0.0114495},
            {61.5709, -47.6158, 3.4295, -0.0316429}},
            {{36.6259, -27.4064, 1.16617, -0.00604629},
                {50.3751, -37.5848, 2.19621, -0.0169241},
                {35.1563, -26.514, 1.09795, -0.00545864}},
            {{27.2367, -20.3068, 0.517752, -0.000335432},
                {39.0489, -28.6903, 1.24306, -0.0065226},
                {41.0208, -30.0339, 1.30776, -0.00626721}},
            {{29.261, -21.7041, 0.613556, -0.000774652},
                {39.5304, -29.1388, 1.34116, -0.00823818},
                {44.5313, -33.4056, 1.77581, -0.0123965}},
            {{36.5659, -25.119, 0.714074, -2.65397e-11},
                {31.6524, -22.6934, 0.613977, -5.46634e-10},
                {34.7312, -24.9901, 0.749061, -1.22922e-09}},
            {{33.154, -23.8803, 0.685794, -1.13236e-10},
                {42.6731, -31.0799, 1.40425, -0.00730816},
                {46.4732, -35.6988, 2.10144, -0.0164771}}},
        {{{2.40451, -15.0848, 1.05504, -0.0103356},
            {8.93825, -16.5995, 0.925874, -0.00767902},
            {7.23439e-08, -12.5963, 0.814574, -0.00864749}},
            {{6.2953e-07, -12.6365, 0.732206, -0.00639165},
                {12.6866, -18.7831, 1.0952, -0.00923029},
                {3.12805e-07, -12.5395, 0.795535, -0.00828991}},
            {{2.69495, -14.8778, 1.00751, -0.00975373},
                {6.05446, -14.6778, 0.767457, -0.00636729},
                {3.94741e-07, -11.1038, 0.524109, -0.00471514}},
            {{2.31558e-07, -11.5073, 0.494316, -0.00303611},
                {5.66995, -14.5948, 0.740956, -0.00561851},
                {4.40475e-06, -9.57062, 0.20354, -0.000213213}},
            {{2.74277e-08, -13.3573, 0.886651, -0.00857992},
                {9.98849e-05, -11.524, 0.531486, -0.00391441},
                {8.50811e-07, -9.72224, 0.240264, -0.000781498}},
            {{6.9021e-08, -11.8859, 0.53864, -0.00325092},
                {10.0169, -16.9153, 0.921593, -0.00752414},
                {9.90518e-07, -11.9578, 0.697029, -0.00717645}}},
        {{{6.87966e-10, -12.8497, 0.757379, -0.00651612},
            {16.2087, -19.3776, 0.951508, -0.00645029},
            {14.513, -18.8625, 1.05542, -0.00918985}},
            {{1.07197e-07, -12.5469, 0.703086, -0.00585238},
                {0.0871522, -9.22628, 0.159628, -0.000343326},
                {12.1181, -17.5575, 0.940249, -0.00788125}},
            {{2.10191e-09, -12.2782, 0.661926, -0.00555279},
                {12.5105, -17.9998, 0.951807, -0.00732845},
                {12.8043, -17.8322, 0.972401, -0.00841528}},
            {{8.11926e-10, -12.7225, 0.737941, -0.00647355},
                {7.50649, -15.987, 0.889398, -0.00729282},
                {0.174541, -10.0266, 0.306882, -0.00186093}},
            {{3.81202e-09, -12.0926, 0.598943, -0.00430458},
                {8.72368, -17.2511, 1.06348, -0.00953327},
                {1.5205, -9.86713, 0.183806, -6.40377e-12}},
            {{1.37378e-07, -12.9247, 0.769722, -0.00664936},
                {8.53877, -16.6167, 0.946138, -0.00788745},
                {8.47417, -14.3897, 0.581492, -0.00387111}}},
        {{{2.50079e-07, -12.5209, 0.678491, -0.00528954},
            {12.6171, -18.4205, 1.01802, -0.00807702},
            {10.4903, -18.0981, 1.10546, -0.00971519}},
            {{5.87069e-07, -12.0075, 0.585538, -0.00416654},
                {11.1348, -17.5468, 0.943652, -0.00729083},
                {0.949201, -10.5869, 0.267536, -6.04802e-05}},
            {{1.14857, -11.1478, 0.345528, -0.000841836},
                {10.9482, -17.1647, 0.909605, -0.00722404},
                {8.7569e-08, -10.4446, 0.316302, -0.00101964}},
            {{1.09759e-06, -11.5019, 0.48435, -0.00277852},
                {0.637937, -10.7065, 0.316211, -0.000801127},
                {5.67144e-07, -12.88, 0.831252, -0.00835441}},
            {{1.68853, -11.2582, 0.308152, -7.81686e-12},
                {9.44238, -17.1892, 1.00561, -0.00864837},
                {1.20713e-07, -12.2246, 0.669321, -0.0057622}},
            {{0.00217558, -10.8858, 0.347928, -0.000790679},
                {11.8583, -17.6423, 0.923581, -0.00703041},
                {3.24078, -13.4024, 0.668777, -0.00504175}}},
        {{{6.04158, -16.8155, 1.13335, -0.0105359},
            {8.24786, -17.0204, 1.05097, -0.00941875},
            {11.7617, -17.202, 0.864472, -0.00649032}},
            {{3.70947, -13.0663, 0.513818, -0.00222627},
                {16.7022, -21.9618, 1.42869, -0.012705},
                {6.8993, -14.8192, 0.740813, -0.00585407}},
            {{2.18472e-06, -11.9461, 0.583354, -0.00423414},
                {6.51489e-07, -10.5669, 0.353028, -0.00166977},
                {12.5113, -16.5038, 0.709888, -0.00471964}},
            {{0.812719, -11.3245, 0.390183, -0.00134086},
                {2.97251, -11.9374, 0.338592, -4.36096e-13},
                {13.8844, -17.5707, 0.818446, -0.00581811}},
            {{1.55496, -14.4569, 0.949497, -0.00857237},
                {0.34359, -10.5041, 0.286497, -0.000346977},
                {14.4141, -18.7457, 1.01652, -0.00845189}},
            {{1.26317e-08, -11.1424, 0.434251, -0.00236267},
                {6.58119, -15.8546, 0.930324, -0.00801288},
                {4.41865, -11.1991, 0.234652, -7.43723e-10}}},
        {{{6.87926, -12.8949, 0.334733, -6.38494e-06},
            {35.2336, -32.2007, 2.21489, -0.020555},
            {6.80949, -16.8945, 1.19056, -0.0127558}},
            {{0.95782, -12.4625, 0.599979, -0.00405342},
                {20.4051, -23.1936, 1.42408, -0.0120792},
                {10.277, -16.1457, 0.785186, -0.00612069}},
            {{0.236196, -11.6165, 0.458613, -0.002018},
                {12.8771, -19.6785, 1.26163, -0.0115917},
                {5.21194e-08, -12.551, 0.78718, -0.00794713}},
            {{8.40778, -14.9001, 0.534967, -0.00147246},
                {15.9376, -20.9945, 1.2908, -0.0110556},
                {10.4773, -16.2238, 0.783386, -0.00593478}},
            {{3.21187, -12.1221, 0.348938, -8.70415e-14},
                {13.8983, -19.1128, 1.04727, -0.00797426},
                {11.6342, -18.8428, 1.18853, -0.0107619}},
            {{3.7311, -12.4292, 0.419345, -0.00134704},
                {6.92884, -13.2494, 0.391862, -0.000767396},
                {5.5939, -14.4175, 0.729195, -0.00568477}}}};
    
    
    double theta_DCr = 5000;
    double phi_DCr_raw = 5000;
    
    switch (region)
    {
        case 1:
            theta_DCr = 180 / Pival * acos(part_DC_c1z / sqrt(pow(part_DC_c1x, 2) + pow(part_DC_c1y, 2) + pow(part_DC_c1z, 2)));
            phi_DCr_raw = 180 / Pival * atan2(part_DC_c1y / sqrt(pow(part_DC_c1x, 2) + pow(part_DC_c1y, 2) + pow(part_DC_c1z, 2)), part_DC_c1x / sqrt(pow(part_DC_c1x, 2) + pow(part_DC_c1y, 2) + pow(part_DC_c1z, 2)));
            break;
            
        case 2:
            theta_DCr = 180 / Pival * acos(part_DC_c2z / sqrt(pow(part_DC_c2x, 2) + pow(part_DC_c2y, 2) + pow(part_DC_c2z, 2)));
            phi_DCr_raw = 180 / Pival * atan2(part_DC_c2y / sqrt(pow(part_DC_c2x, 2) + pow(part_DC_c2y, 2) + pow(part_DC_c2z, 2)), part_DC_c2x / sqrt(pow(part_DC_c2x, 2) + pow(part_DC_c2y, 2) + pow(part_DC_c2z, 2)));
            break;
            
        case 3:
            theta_DCr = 180 / Pival * acos(part_DC_c3z / sqrt(pow(part_DC_c3x, 2) + pow(part_DC_c3y, 2) + pow(part_DC_c3z, 2)));
            phi_DCr_raw = 180 / Pival * atan2(part_DC_c3y / sqrt(pow(part_DC_c3x, 2) + pow(part_DC_c3y, 2) + pow(part_DC_c3z, 2)), part_DC_c3x / sqrt(pow(part_DC_c3x, 2) + pow(part_DC_c3y, 2) + pow(part_DC_c3z, 2)));
            break;
            
        default:
            return false;
            break;
    }
    
    double phi_DCr = 5000;
    if (part->trk(DC)->getSector() == 1) phi_DCr = phi_DCr_raw;
    if (part->trk(DC)->getSector() == 2) phi_DCr = phi_DCr_raw - 60;
    if (part->trk(DC)->getSector() == 3) phi_DCr = phi_DCr_raw - 120;
    if (part->trk(DC)->getSector() == 4 && phi_DCr_raw > 0) phi_DCr = phi_DCr_raw - 180;
    if (part->trk(DC)->getSector() == 4 && phi_DCr_raw < 0) phi_DCr = phi_DCr_raw + 180;
    if (part->trk(DC)->getSector() == 5) phi_DCr = phi_DCr_raw + 120;
    if (part->trk(DC)->getSector() == 6) phi_DCr = phi_DCr_raw + 60;
    
    int pid = 0;
    
    switch (part->getPid()){
            
        case 11:
            pid = 0;
            break;
        case 2212:
            pid = 1;
            break;
        case 211:
            pid = 2;
            break;
        case -211:
            pid = 3;
            break;
        case 321:
            pid = 4;
            break;
        case -321:
            pid = 5;
            if (part->trk(DC)->getSector() == 6 || (part->trk(DC)->getSector()== 5 && part->trk(DC)->getSector() == 3))  // use K+ cuts in some cases
            {
                pid = 4;
            }
            break;
            
        default:
            return false;
            break;
    }
    
    double calc_phi_min = minparams[pid][part->trk(DC)->getSector()][region][0] + minparams[pid][part->trk(DC)->getSector()][region][1] * std::log(theta_DCr)
    + minparams[pid][part->trk(DC)->getSector()][region][2] * theta_DCr + minparams[pid][part->trk(DC)->getSector()][region][3] * theta_DCr * theta_DCr;
    
    double calc_phi_max = maxparams[pid][part->trk(DC)->getSector()][region][0] + maxparams[pid][part->trk(DC)->getSector()][region][1] * std::log(theta_DCr)
    + maxparams[pid][part->trk(DC)->getSector()][region][2] * theta_DCr + maxparams[pid][part->trk(DC)->getSector()][region][3] * theta_DCr * theta_DCr;
    
    return ((phi_DCr > calc_phi_min) && (phi_DCr < calc_phi_max));
}

bool minimal_PCAL_energy_deposition(region_part_ptr part){
    
    bool tight = false;
    bool medium = true;
    bool loose = false;
    bool inbending=true;
    bool outbending=false;
    
    int PCAL_sector=part->cal(PCAL)->getSector();
    double PCAL_energy=part->cal(PCAL)->getEnergy();
    double ECAL_total= (part->cal(PCAL)->getEnergy())+(part->cal(ECIN)->getEnergy())+(part->cal(ECOUT)->getEnergy());
    
    double PCAL_loose[] = {0.06,  0.06,  0.06, 0.06,  0.06,  0.06 };
    double PCAL_medium[]   = {0.07,  0.07,  0.07, 0.07,  0.07,  0.07 };
    double PCAL_tight[]   = {0.09,  0.09,  0.09, 0.09,  0.09,  0.09 };
    
    double PCAL_energy_cut=0.0;
    for(Int_t k = 0; k < 6; k++){
        if(PCAL_sector-1 == k){
            if(tight==true){
                PCAL_energy_cut=PCAL_tight[k];
            }
            if(medium==true){
                PCAL_energy_cut=PCAL_medium[k];
            }
            
            if(loose==true){
                PCAL_energy_cut=PCAL_loose[k];
            }
        }
    }
    
    if(PCAL_energy > PCAL_energy_cut)return true;
    else return false;
    
}

//use Z-vetex cut between -13 < Z_vertex < 10 (for simplicity) independent of the sectors
bool z_vertex_cut(region_part_ptr part){
    
    bool tight = true;
    bool medium = false;
    bool loose = false;
    bool inbending=true;
    bool outbending=false;
    
    
    double Z_vertex= part->par()->getVz();
    
    double Z_vertex_low=-13.0;
    double Z_vertex_high=12.0;
    
    if(Z_vertex > Z_vertex_low && Z_vertex < Z_vertex_high)return true;
    else return false;
    
}


//Vertex differnce cuts for proton, pi+ and pi- (-20 < Vertex_differnce < 20)
/*
 bool vertex_differnce_cut(region_part_ptr part, int index){
 
 int pid=par->part()->getPid();
 double electronVz, protonVz,piplusVz,piminusVz;
 if(pid==11) electronVz=part->par()->getVz();
 if(pid==2212) protonVz=part->par()->getVz();
 if(pid==211) piplusVz=part->par()->getVz()
 if(pid==-211) piminusVz=part->par()->getVz();
 
 if(index==1){
 double proton_Vz_difference= electronVz-protonVz;
 }
 
 
 
 
 }*/
/*bool momentum_correction(region_part_ptr part){
 
 int sec=part->cal(PCAL)->getSector();
 
 float px = part->par()->getPx();
 float py = part->par()->getPy();
 float pz = part->par()->getPz();
 
 double xx[] = {0.0265368, 0.0448258, 0.0663685, 0.048579, 0.00255056, -0.0956177,
 0.0131988, 0.00650346, -0.00539524, -0.0678265, -0.0212479, 0.0311322,
 0.0117689, 0.0129722, 0.0143328, -0.0158383, -0.0544898, -0.0934764,
 0.020875, 0.0337656, 0.0502739, 0.0190154, -0.00195974, -0.0492416,
 0.00306965, 0.000586585, -0.00176532, -0.0247403, -0.010975, -0.00474736,
 -0.00417779, 0.00555585, 0.0136674, -0.0217805, -0.0619875, -0.096519,
 0.0434163, 0.0158384, 0.0141099, -0.0818153, -0.012961, -0.0108725,
 0.0135067, 0.0113164, 0.00857531, -0.0577991, -0.0468577, -0.0350849,
 0.00503738, 0.00956591, 0.0142982, -0.0922083, -0.0936707, -0.0968399,
 0.0199686, 0.00803408, 0.0228017, -0.0223379, 0.0305244, 0.0234411,
 0.01538, 0.0108718, 0.00322224, -0.0662346, -0.0379318, -0.0120417,
 0.00892442, 0.0131597, 0.0173936, -0.0485261, -0.0656887, -0.0952142,
 0.000587639, 0.0299187, 0.0490066, 0.0119492, 0.00704546, -0.0217605,
 0.0129687, 0.00812209, 0.000975356, -0.0835302, -0.0603761, -0.0256162,
 0.0168102, 0.018545, 0.0197326, -0.0646219, -0.0783352, -0.0921586,
 0.0328061, 0.0467324, 0.0579207, 0.00787372, -0.0098708, -0.0594644,
 0.0130907, 0.0114952, 0.00910189, -0.0366264, -0.0212249, -7.99685e-05,
 -0.000368554, 0.00125095, 0.00325831, 0.0157934, -0.0156403, -0.0523252};
 
 double pars[6][3][6];
 int ipar=0;
 for(int isec=0;isec<6;isec++)
 for(int ivec=0;ivec<3;ivec++) {
 double dp11=xx[ipar++], dp13=xx[ipar++], dp15=xx[ipar++];
 double dp91=xx[ipar++], dp93=xx[ipar++], dp95=xx[ipar++];
 
 pars[isec][ivec][0] = (9*dp11 - 18*dp13 + 9*dp15 - dp91 + 2*dp93 - dp95)/6400.;
 pars[isec][ivec][1] = (-9*dp11)/80. + (27*dp13)/160. - (9*dp15)/160. + dp91/80. - (3*dp93)/160. + dp95/160.;
 pars[isec][ivec][2] = (135*dp11)/64. - (45*dp13)/32. + (27*dp15)/64. - (15*dp91)/64. + (5*dp93)/32. - (3*dp95)/64.;
 pars[isec][ivec][3] = -dp11/6400. + dp13/3200. - dp15/6400. + dp91/6400. - dp93/3200. + dp95/6400.;
 pars[isec][ivec][4] = dp11/80. - (3*dp13)/160. + dp15/160. - dp91/80. + (3*dp93)/160. - dp95/160.;
 pars[isec][ivec][5] = (-15*dp11)/64. + (5*dp13)/32. - (3*dp15)/64. + (15*dp91)/64. - (5*dp93)/32. + (3*dp95)/64.;
 }
 
 auto dpp = [&](float px, float py, float pz, int sec, int ivec) {
 double pp = sqrt(px*px + py*py + pz*pz);
 double fi = TMath::RadToDeg()*atan2(py,px);
 fi = fi + (fi<0 && sec>1)*360 - (sec-1)*60;
 
 fi += ivec==2 ? 10 : 25 + ivec*25;
 
 double a0=pars[sec-1][ivec][0],
 b0=pars[sec-1][ivec][1],
 c0=pars[sec-1][ivec][2],
 a1=pars[sec-1][ivec][3],
 b1=pars[sec-1][ivec][4],
 c1=pars[sec-1][ivec][5];
 
 double dp = (a0*fi*fi + b0*fi + c0) + (a1*fi*fi + b1*fi + c1)*pp;
 return dp/pp;
 };
 
 if(pid==11) double fe = dpp(px,py,pz,sec,0) + 1;
 return fe;
 
 }*/




