#include "src/Structs.h"
#include "src/HipoBankInterface.C"
#include "src/CutManager.C"
#include "src/Constants.h"

int find_mc_match(std::vector<part> rec_parts, part mc_part);

//int hipo2tree(const char * input_hipo_file = "",
//              const char * output_root_file = ""){
int hipo2tree(const char * input_hipo_file = "mc_rga_inbending/45nA_job_3051_0.hipo",
              const char * output_root_file = "mc_3051_0_test.root"){
  
  // Open TTree and declare branches
  // -------------------------------------
  TFile *outfile = new TFile(output_root_file,"RECREATE");
  TTree* tree = new TTree("EventTree","EventTree");
    
  int Nmax = 100; // maximum number of particles
  int run = 0;
  int torus = 0;
  float beamE = 0; // To be set in event loop
  double x,Q2,y,W,s;
  int pindex[Nmax], status[Nmax], pid[Nmax], truepid[Nmax],trueparentpid[Nmax],trueparentid[Nmax];
  double px[Nmax], py[Nmax], pz[Nmax], p[Nmax], E[Nmax], m[Nmax], vz[Nmax], beta[Nmax],chi2[Nmax];
  double theta[Nmax], phi[Nmax];
    
  int pcal_sector[Nmax], ecin_sector[Nmax], ecout_sector[Nmax];
  double pcal_x[Nmax], pcal_y[Nmax], pcal_z[Nmax];
  double ecin_x[Nmax], ecin_y[Nmax], ecin_z[Nmax];
  double ecout_x[Nmax], ecout_y[Nmax], ecout_z[Nmax];
  double pcal_e[Nmax], pcal_lu[Nmax], pcal_lv[Nmax], pcal_lw[Nmax], pcal_m2u[Nmax], pcal_m2v[Nmax], pcal_m2w[Nmax];
  double ecin_e[Nmax], ecin_lu[Nmax], ecin_lv[Nmax], ecin_lw[Nmax], ecin_m2u[Nmax], ecin_m2v[Nmax], ecin_m2w[Nmax];
  double ecout_e[Nmax], ecout_lu[Nmax], ecout_lv[Nmax], ecout_lw[Nmax], ecout_m2u[Nmax], ecout_m2v[Nmax], ecout_m2w[Nmax];
    
  tree->Branch("Nmax",&Nmax,"Nmax/I");
  tree->Branch("run",&run,"run/I");
  tree->Branch("x",&x,"x/D");
  tree->Branch("Q2",&Q2,"Q2/D");
  tree->Branch("y",&y,"y/D");
  tree->Branch("W",&W,"W/D");
  tree->Branch("pindex", pindex, "pindex[Nmax]/I");
  tree->Branch("px", px, "px[Nmax]/D");
  tree->Branch("py", py, "py[Nmax]/D");
  tree->Branch("pz", pz, "pz[Nmax]/D");
  tree->Branch("p", p, "p[Nmax]/D");
  tree->Branch("E", E, "E[Nmax]/D");
  tree->Branch("m", m, "m[Nmax]/D");
  tree->Branch("pid", pid, "pid[Nmax]/I");
  tree->Branch("truepid", truepid, "truepid[Nmax]/I");
  tree->Branch("trueparentpid", trueparentpid, "trueparentpid[Nmax]/I");
  tree->Branch("trueparentid", trueparentid, "trueparentid[Nmax]/I");
  tree->Branch("theta", theta, "theta[Nmax]/D");
  tree->Branch("phi", phi, "phi[Nmax]/D");
    
  tree->Branch("pcal_sector", pcal_sector, "pcal_sector[Nmax]/I");
  tree->Branch("ecin_sector", ecin_sector, "ecin_sector[Nmax]/I");
  tree->Branch("ecout_sector", ecout_sector, "ecout_sector[Nmax]/I");
  tree->Branch("pcal_x", pcal_x, "pcal_x[Nmax]/D");
  tree->Branch("pcal_y", pcal_y, "pcal_y[Nmax]/D");
  tree->Branch("pcal_z", pcal_z, "pcal_z[Nmax]/D");
  tree->Branch("ecin_x", ecin_x, "ecin_x[Nmax]/D");
  tree->Branch("ecin_y", ecin_y, "ecin_y[Nmax]/D");
  tree->Branch("ecin_z", ecin_z, "ecin_z[Nmax]/D");
  tree->Branch("ecout_x", ecout_x, "ecout_x[Nmax]/D");
  tree->Branch("ecout_y", ecout_y, "ecout_y[Nmax]/D");
  tree->Branch("ecout_z", ecout_z, "ecout_z[Nmax]/D");
  tree->Branch("pcal_e", pcal_e, "pcal_e[Nmax]/D");
  tree->Branch("pcal_lu", pcal_lu, "pcal_lu[Nmax]/D");
  tree->Branch("pcal_lv", pcal_lv, "pcal_lv[Nmax]/D");
  tree->Branch("pcal_lw", pcal_lw, "pcal_lw[Nmax]/D");
  tree->Branch("pcal_m2u", pcal_m2u, "pcal_m2u[Nmax]/D");
  tree->Branch("pcal_m2v", pcal_m2v, "pcal_m2v[Nmax]/D");
  tree->Branch("pcal_m2w", pcal_m2w, "pcal_m2w[Nmax]/D");
  tree->Branch("ecin_e", ecin_e, "ecin_e[Nmax]/D");
  tree->Branch("ecin_lu", ecin_lu, "ecin_lu[Nmax]/D");
  tree->Branch("ecin_lv", ecin_lv, "ecin_lv[Nmax]/D");
  tree->Branch("ecin_lw", ecin_lw, "ecin_lw[Nmax]/D");
  tree->Branch("ecin_m2u", ecin_m2u, "ecin_m2u[Nmax]/D");
  tree->Branch("ecin_m2v", ecin_m2v, "ecin_m2v[Nmax]/D");
  tree->Branch("ecin_m2w", ecin_m2w, "ecin_m2w[Nmax]/D");
  tree->Branch("ecout_e", ecout_e, "ecout_e[Nmax]/D");
  tree->Branch("ecout_lu", ecout_lu, "ecout_lu[Nmax]/D");
  tree->Branch("ecout_lv", ecout_lv, "ecout_lv[Nmax]/D");
  tree->Branch("ecout_lw", ecout_lw, "ecout_lw[Nmax]/D");
  tree->Branch("ecout_m2u", ecout_m2u, "ecout_m2u[Nmax]/D");
  tree->Branch("ecout_m2v", ecout_m2v, "ecout_m2v[Nmax]/D");
  tree->Branch("ecout_m2w", ecout_m2w, "ecout_m2w[Nmax]/D");
    
  // Configure CLAS12 Reader and HipoChain
  // -------------------------------------
  clas12root::HipoChain _chain;
  clas12::clas12reader *_config_c12{nullptr};

  _chain.Add(input_hipo_file);
  _config_c12=_chain.GetC12Reader();
    
  // Configure PIDs for final state
  // -------------------------------------
  _config_c12->addAtLeastPid(11,1); // At least 1 electron
  _config_c12->addAtLeastPid(22,2); // At least 2 photons
 
  // Add RUN::config bank
  // -------------------------------------
  int _idx_RUNconfig = _config_c12->addBank("RUN::config");
  int _irun = _config_c12->getBankOrder(_idx_RUNconfig,"run");
  int _itorus = _config_c12->getBankOrder(_idx_RUNconfig,"torus");  
  
  // Establish CLAS12 event parser
  // -------------------------------------
  auto &_c12=_chain.C12ref();
    
  // Create Objects/Handlers
  // -------------------------------------
  HipoBankInterface _hipoInterface = HipoBankInterface(_c12);
  CutManager _cm = CutManager();  
    
  // Create particle struct vector
  // -------------------------------------
  std::vector<part> vec_particles;
  
  // Create init electron, init target, q, and scattered electron vectors
  // -------------------------------------
  TLorentzVector init_electron; // Has Pz and E set in event loop
  TLorentzVector init_target;
  TLorentzVector qVector, eVector;
  init_target.SetPxPyPzE(0,0,0,0.938272);

  // Loop over events in the hipo file
  // -------------------------------------
  int while_loop_index = 0;
  int tree_entries = 0;
  bool is_monte_carlo = false; // Monte Carlo Runs assumed to have RUN::config run==11
    
  while(_chain.Next()==true && while_loop_index < 10000){
    if(while_loop_index%10000==0 && while_loop_index!=0){
      std::cout << while_loop_index << " events read | " << tree_entries << " passed cuts (" << tree_entries*100.0/while_loop_index << "%)" << std::endl;
    }
    while_loop_index++;
      
    // Clear vectors
    vec_particles.clear();
      
    // Get run number and store in CutManager "_cm"
    // Important if one wants to define run dependent cuts
    // ----------------------------------------------------
    run = _c12->getBank(_idx_RUNconfig)->getInt(_irun,0);
    if(run==11) is_monte_carlo=true;
    torus = _c12->getBank(_idx_RUNconfig)->getFloat(_itorus,0);
    
    _cm.set_run(run);
    _cm.set_torus(torus);
      
    // Set the init_electron TLorentzVector based off run period
    // ----------------------------------------------------
    beamE = runBeamEnergy(run);
    init_electron.SetPxPyPzE(0,0,sqrt(beamE*beamE-0.000511*0.000511),beamE);
    s = pow(0.938272,2)+pow(0.000511,2)+2*0.938272*beamE;
      
    // Loop over reconstructed particles
    // -------------------------------------------------------
    auto particles=_c12->getDetParticles();
    for(unsigned int idx = 0 ; idx < particles.size() ; idx++){
      // Create new part struct
      part partstruct;
      // Extract each particle from event one-at-a-time
      // -------------------------------------------------------
      auto particle = particles.at(idx);
      partstruct.pid = particle->getPid(); 
      partstruct.chi2 = particle->getChi2Pid();
      partstruct.theta = particle->getTheta();
      partstruct.phi = particle->getPhi();
      partstruct.p = particle->getP();
      partstruct.px = partstruct.p*cos(partstruct.phi)*sin(partstruct.theta);
      partstruct.py = partstruct.p*sin(partstruct.phi)*sin(partstruct.theta);
      partstruct.pz = partstruct.p*cos(partstruct.theta);

      if(partstruct.pid!=22)
        partstruct.m = particle->getPdgMass();
      else
        partstruct.m = 0;
      partstruct.beta = particle->getBeta();
      partstruct.pindex = particle->getIndex();
      partstruct.vz = particle->par()->getVz();
      partstruct.status = particle->getStatus();
      partstruct.E = sqrt(partstruct.p*partstruct.p+partstruct.m*partstruct.m);
    
      // Ensure hadrons are not in CD
      if (partstruct.pid == 2212 || partstruct.pid == -2212 ||
	  partstruct.pid == 2112 ||
	  partstruct.pid == -321 || partstruct.pid == -211 ||
	  partstruct.pid == 211 || partstruct.pid == 321) {
	if(partstruct.status>=4000 && partstruct.status<5000)
	  continue;
      }
        
      // Add particle to the list
      _hipoInterface.loadBankData(_c12,partstruct);
      vec_particles.push_back(partstruct);
    }
      
    // Code for determine the scattered electron from REC::Particle
    // --> Find pid==11 particle with largest energy
    // -->   If no electron is found, skip
    // --> Check if the status of the maximum energy electron is in FD
    // -->   Skip if not (i.e. always skip events if the max energy electron was not in FD)
    // --> Set that particle as the scattered electron
    int idx_e=-1;
    double max_energy = -1; 
    for (int i = 0; i < vec_particles.size(); i++) {
      part partstruct = vec_particles[i];
      // check if the particle is an electron
      if (partstruct.pid == 11) {
        // compare energy with the current maximum and update if necessary
        if (partstruct.E > max_energy) {
          max_energy = partstruct.E;
          idx_e=i;
         }
       }
    }
    
    if(idx_e==-1){
        continue; // No scattered electron passing the above conditions found
    }
    
    // Get the scattered electron
    part scattered_electron = vec_particles[idx_e];
    if((scattered_electron.status <= -3000 || scattered_electron.status > -2000)) continue; // Max E electron has bad status
    vec_particles[idx_e].is_scattered_electron=true;  
    
    // Get the scattered electron 4-vector
    eVector.SetPxPyPzE(scattered_electron.px,scattered_electron.py,scattered_electron.pz,scattered_electron.E);
    // Get the virtual photon vector
    qVector = init_electron-eVector;
    
    // Calculate event kinematics
    y = (init_electron.E()-eVector.E())/init_electron.E();
    if(y>0.8) continue;
    Q2 = abs(qVector.M2());
    W = (qVector+init_target).M();
    x = Q2/s/y;
    
    // Apply the cuts to make a new vec_particles
    // Calls the CutManager which parses through the vector
    // and makes relevant cuts for each particle
    // 
    // Return a list of filtered particles
    // -------------------------------------------------
      
    vec_particles = _cm.filter_particles(vec_particles);
      
    // Check list of cut particles to see if event is worth keeping
    // -------------------------------------------------------------
    int num_scattered_e=0;
    int num_gamma=0;
    for(part particle : vec_particles){
      if(particle.pid==11 && particle.is_scattered_electron==1) num_scattered_e++;
      if(particle.pid==22) num_gamma++;
    }
    if(num_scattered_e==0 || num_gamma<2){ continue; } // Skip events w/o e- or 2+ gammas
    
    // Loop over all Monte Carlo particles
    // -------------------------------------
    if(is_monte_carlo){
        auto mcparticles=_c12->mcparts();
        for(int idx = 0 ; idx < mcparticles->getRows(); idx++){
            if(mcparticles->getType(idx)!=1) // Reject non-final state
                {continue;} 
            part mcpartstruct;
            mcpartstruct.truepid = mcparticles->getPid(idx);
            mcpartstruct.px = mcparticles->getPx(idx);
            mcpartstruct.py = mcparticles->getPy(idx);
            mcpartstruct.pz = mcparticles->getPz(idx);
            mcpartstruct.m = mcparticles->getMass(idx);
            mcpartstruct.E = sqrt( pow(mcpartstruct.px,2) + pow(mcpartstruct.py,2) + pow(mcpartstruct.pz,2) + pow(mcpartstruct.m,2) ); 
            mcpartstruct.theta = abs(atan(sqrt(pow(mcpartstruct.px,2)+pow(mcpartstruct.py,2))/mcpartstruct.pz));
            mcpartstruct.phi = atan2(mcpartstruct.py,mcpartstruct.px);
            mcpartstruct.trueparentid = mcparticles->getParent(idx)-1;
            mcpartstruct.trueparentpid = mcparticles->getPid(mcpartstruct.trueparentid);
            
            // Does this Monte Carlo particle's final state kinematics match with a reconstructed particle?
            int irecpart = find_mc_match(vec_particles, mcpartstruct);
            if(irecpart==-1) continue; // No match found
            
            // If it does match, then save some Monte Carlo info
            vec_particles.at(irecpart).truepid=mcpartstruct.truepid;
            vec_particles.at(irecpart).trueparentpid=mcpartstruct.trueparentpid;
            vec_particles.at(irecpart).trueparentid=mcpartstruct.trueparentid;
        }
    }
      
    // Loop over all particles and fill variables for the TTree
    // --------------------------------------------------------
    Nmax=vec_particles.size();
    for(int i = 0; i < vec_particles.size(); i++) {
      part par = vec_particles[i];
        
      pindex[i] = par.pindex;
      status[i] = par.status;
      px[i] = par.px;
      py[i] = par.py;
      pz[i] = par.pz;
      p[i] = par.p;
      E[i] = par.E;
      pid[i] = par.pid;
      truepid[i] = par.truepid;
      trueparentpid[i] = par.trueparentpid;
      trueparentid[i] = par.trueparentid;
      vz[i] = par.vz;
      chi2[i] = par.chi2;
      beta[i] = par.beta;
      m[i] = par.m;
      theta[i] = par.theta;
      phi[i] = par.phi;
      pcal_sector[i] = par.pcal_sector;
      pcal_e[i] = par.pcal_e;
      pcal_x[i] = par.pcal_x;
      pcal_y[i] = par.pcal_y;
      pcal_z[i] = par.pcal_z;
      pcal_lu[i] = par.pcal_lu;
      pcal_lv[i] = par.pcal_lv;
      pcal_lw[i] = par.pcal_lw;
      pcal_m2u[i] = par.pcal_m2u;
      pcal_m2v[i] = par.pcal_m2v;
      pcal_m2w[i] = par.pcal_m2w;
      ecin_sector[i] = par.ecin_sector;
      ecin_e[i] = par.ecin_e;
      ecin_x[i] = par.ecin_x;
      ecin_y[i] = par.ecin_y;
      ecin_z[i] = par.ecin_z;
      ecin_lu[i] = par.ecin_lu;
      ecin_lv[i] = par.ecin_lv;
      ecin_lw[i] = par.ecin_lw;
      ecin_m2u[i] = par.ecin_m2u;
      ecin_m2v[i] = par.ecin_m2v;
      ecin_m2w[i] = par.ecin_m2w;
      ecout_sector[i] = par.ecout_sector;
      ecout_e[i] = par.ecout_e;
      ecout_x[i] = par.ecout_x;
      ecout_y[i] = par.ecout_y;
      ecout_z[i] = par.ecout_z;
      ecout_lu[i] = par.ecout_lu;
      ecout_lv[i] = par.ecout_lv;
      ecout_lw[i] = par.ecout_lw;
      ecout_m2u[i] = par.ecout_m2u;
      ecout_m2v[i] = par.ecout_m2v;
      ecout_m2w[i] = par.ecout_m2w;
    }
    tree->Fill();
    tree_entries++;
  }
  
  outfile->cd();
  tree->Write();
  outfile->Close();
  return 0;

}


// Given the list of REC::Particle's, determine if the provided mc_part matches based on certain criteria
// Return the index of that reconstructed particle
int find_mc_match(std::vector<part> rec_parts, part mc_part){
    int ipar=0;
    for(auto rec_part: rec_parts){
        double rec_theta = rec_part.theta;
        double rec_phi = rec_part.phi;
        double rec_E = rec_part.E;
        
        double mc_theta = mc_part.theta;
        double mc_phi = mc_part.phi;
        double mc_E = mc_part.E;
        
        double dth = abs(rec_theta-mc_theta)*180/3.14159265;
        double dphi = abs(rec_phi-mc_phi)*180/3.14159265;
        double dE = abs(rec_E-mc_E);
        
        if (dth<2 && (dphi<4 || abs(dphi-2*3.14159265)<4) && dE<1) return ipar;
        ipar++;
    }
    
    // no match found
    return -1;
}
