#include "src/Constants.h"
double calc_Pt(TLorentzVector, TLorentzVector, TLorentzVector);

int buildDiphotons(const char * input_file=""){
    // Read the TFile
    TFile *f = new TFile(input_file,"UPDATE");
    // Read the TTree
    TTree *EventTree = (TTree*)f->Get("EventTree");
    
    double x, Q2, W,y;
    int Nmax=100;
    int run=0;
    double px[Nmax], py[Nmax], pz[Nmax], E[Nmax], theta[Nmax], phi[Nmax];
    double p_gamma[Nmax];
    int pid[Nmax];
    EventTree->SetBranchAddress("run",&run);
    EventTree->SetBranchAddress("x",&x);
    EventTree->SetBranchAddress("Q2",&Q2);
    EventTree->SetBranchAddress("W",&W);
    EventTree->SetBranchAddress("y",&y);
    EventTree->SetBranchAddress("Nmax",&Nmax);
    EventTree->SetBranchAddress("px",px);
    EventTree->SetBranchAddress("py",py);
    EventTree->SetBranchAddress("pz",pz);
    EventTree->SetBranchAddress("E",E);
    EventTree->SetBranchAddress("theta",theta);
    EventTree->SetBranchAddress("phi",phi);
    EventTree->SetBranchAddress("p_gamma",p_gamma);
    EventTree->SetBranchAddress("pid",pid);
    
    
    //Create the new TTree
    TString treename="diphoton";
    //If dihadron tree already exists, remove it
    if (f->Get(treename)) f->Delete("diphoton*;*");
    TTree *outtree = new TTree(treename,"Diphoton-by-Diphoton info");
    
    double z,pT,M_gg,p_gamma_1,p_gamma_2,E_gamma_1,E_gamma_2,E_ele,th_ele,phi_ele,E_gg,th_gg,phi_gg;
    
    outtree->Branch("run",&run,"run/D");
    outtree->Branch("x",&x,"x/D");
    outtree->Branch("Q2",&Q2,"Q2/D");
    outtree->Branch("y",&y,"y/D");
    outtree->Branch("W",&W,"W/D");
    outtree->Branch("z",&z,"z/D");
    outtree->Branch("pT",&pT,"pT/D");
    outtree->Branch("M_gg",&M_gg,"M_gg/D");
    outtree->Branch("p_gamma_1",&p_gamma_1,"p_gamma_1/D");
    outtree->Branch("p_gamma_2",&p_gamma_2,"p_gamma_2/D");
    outtree->Branch("E_gamma_1",&E_gamma_1,"E_gamma_1/D");
    outtree->Branch("E_gamma_2",&E_gamma_2,"E_gamma_2/D");
    outtree->Branch("E_ele",&E_ele,"E_ele/D");
    outtree->Branch("th_ele",&th_ele,"th_ele/D");
    outtree->Branch("phi_ele",&phi_ele,"phi_ele/D");
    outtree->Branch("E_gg",&E_gg,"E_gg/D");
    outtree->Branch("th_gg",&th_gg,"th_gg/D");
    outtree->Branch("phi_gg",&phi_gg,"phi_gg/D");
    
    // Initial particles
    TLorentzVector init_electron(0,0,0,0); // To be set one run is found
    TLorentzVector init_target(0,0,0,0.938272);
    
    // Photons
    TLorentzVector gamma1,gamma2,diphoton;
    
    // for loop over all events
    int N = EventTree->GetEntries();
    
    for (int ev=0; ev<N; ev++){
        EventTree->GetEntry(ev);
        init_electron.SetE(runBeamEnergy(run));
        init_electron.SetPz(sqrt(init_electron.E()*init_electron.E()-0.000511*0.000511));
        
        //Loop over all particles in the event to find electron
        TLorentzVector electron;
        TLorentzVector q; // virtual photon
        // Maximum energy electron is assumed to be scattered e- (see hipo2tree.C)
        int idx_e=-1;
        double max_e=-1;
        for (int i=0; i<Nmax; i++){
            if(pid[i]==11){
                if(E[i]>max_e){
                    idx_e=i;
                    max_e=E[i];
                }
            }
        }
        electron.SetPxPyPzE(px[idx_e],py[idx_e],pz[idx_e],E[idx_e]);
        q=init_electron-electron;
        
        E_ele = electron.E();
        th_ele = electron.Theta();
        phi_ele = electron.Phi();
        
        // Loop over first set of photons
        for(int i = 0; i<Nmax; i++){
            if(pid[i]!=22) continue;
            gamma1.SetPxPyPzE(px[i],py[i],pz[i],E[i]);
            p_gamma_1=p_gamma[i]; // <--- classifier output
            E_gamma_1=E[i];
            // Loop over second set of photons
            for(int j = i+1; j<Nmax; j++){
                if(pid[j]!=22) continue;
                gamma2.SetPxPyPzE(px[j],py[j],pz[j],E[j]);
                p_gamma_2=p_gamma[j]; // <--- classifier output
                E_gamma_2=E[j];
                
                // Now we have a diphoton, so we will be saving the info the the output TTree "diphoton"
                diphoton=gamma1+gamma2;
                z = (init_target*diphoton)/(init_target*q);
                pT = calc_Pt(q,diphoton,init_target);
                M_gg = diphoton.M();
                E_gg = diphoton.E();
                th_gg = diphoton.Theta();
                phi_gg = diphoton.Phi();
                outtree->Fill();
            }
        }
    }
    f->cd();
    outtree->Write();
    f->Close();
    return 0;
}

double calc_Pt(TLorentzVector q, TLorentzVector p, TLorentzVector init_target){
  TLorentzVector com = q+init_target;
  TVector3 comBOOST = com.BoostVector();
  TLorentzVector qq = q;
  TLorentzVector pp = p;
  qq.Boost(-comBOOST);
  pp.Boost(-comBOOST);
  return pp.Pt(qq.Vect());
}