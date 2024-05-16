#include "src/Structs.h"
#include "src/HipoBankInterface.C"
#include "src/CutManager.C"
#include "src/Constants.h"

int find_mc_match(std::vector<part> rec_parts, part mc_part);

int hipo2tree(const char * input_hipo_file = "",
              const char * output_root_file = "",
              const bool fill_EventTree = true){
  
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
    
    
  //Define the variables "m_g" , "m_ch" , "m_nh" 
  //Should not be changed because the model was trained with this specific set of inputs
  int m_g = 3; // Number of neighboring gammas
  int m_ch = 2; // Number of neighboring charged hadrons (protons, pions, kaons)
  int m_nh = 2; // Number of neighboring neutral hadrons (neutrons)
    
  
    TTree *MLInput = new TTree("MLinput","Nearest neighbor information");

    //Define the branches in MLInput: POI (photon-of-interest)
    //                                Nearest neighbor gammas, charged hadrons, neutral hadrons, electron
    int photon_has_match = 0; // 0 if bkg, 1 if signal
    double gE=0; // Photon energy
    double gEpcal=0; // Photon pcal energy
    double gTheta=0; // Photon angle
    double gm2u=0; // Photon shower shape
    double gm2v=0; // Photon shower shape
    
    double R_e;
    double dE_e;
    
    double R_gamma[m_g]; // Angular distance between calo shower centers
    double dE_gamma[m_g]; // Energy difference
    double Epcal_gamma[m_g]; // Energy deposited in the pcal
    double m2u_gamma[m_g]; // Shower shape variables
    double m2v_gamma[m_g]; // Shower shape variables
    
    double R_ch[m_ch]; // Angular distance between calo shower centers
    double dE_ch[m_ch]; // Energy difference
    double Epcal_ch[m_ch]; // Energy deposited in the pcal
    double m2u_ch[m_ch]; // Shower shape variables
    double m2v_ch[m_ch]; // Shower shape variables
    
    double R_nh[m_nh]; // Angular distance between calo shower centers
    double dE_nh[m_nh]; // Energy difference
    double Epcal_nh[m_nh]; // Energy deposited in the pcal
    double m2u_nh[m_nh]; // Shower shape variables
    double m2v_nh[m_nh]; // Shower shape variables
    
    double num_photons_0_1, num_photons_0_2, num_photons_0_35;
    
    
    // Place TTree Branches
    MLInput->Branch("photon_has_match",&photon_has_match,"photon_has_match/I");

    MLInput->Branch("m_g",&m_g,"m_g/I");
    MLInput->Branch("m_ch",&m_ch,"m_ch/I");
    MLInput->Branch("m_nh",&m_nh,"m_nh/I");
    
    MLInput->Branch("gE",&gE,"gE/D");
    MLInput->Branch("gEpcal",&gEpcal,"gEpcal/D");
    MLInput->Branch("gTheta",&gTheta,"gTheta/D");
    MLInput->Branch("gm2u",&gm2u,"gm2u/D");
    MLInput->Branch("gm2v",&gm2v,"gm2v/D");
    
    MLInput->Branch("R_e",&R_e,"R_e/D");
    MLInput->Branch("dE_e",&dE_e,"dE_e/D");
    
    MLInput->Branch("R_gamma",R_gamma,"R_gamma[m_g]/D");
    MLInput->Branch("dE_gamma",dE_gamma,"dE_gamma[m_g]/D");
    MLInput->Branch("Epcal_gamma",Epcal_gamma,"Epcal_gamma[m_g]/D");
    MLInput->Branch("m2u_gamma",m2u_gamma,"m2u_gamma[m_g]/D");
    MLInput->Branch("m2v_gamma",m2v_gamma,"m2v_gamma[m_g]/D");
    
    MLInput->Branch("R_ch",R_ch,"R_ch[m_ch]/D");
    MLInput->Branch("dE_ch",dE_ch,"dE_ch[m_ch]/D");
    MLInput->Branch("Epcal_ch",Epcal_ch,"Epcal_ch[m_ch]/D");
    MLInput->Branch("m2u_ch",m2u_ch,"m2u_ch[m_ch]/D");
    MLInput->Branch("m2v_ch",m2v_ch,"m2v_ch[m_ch]/D");
    
    MLInput->Branch("R_nh",R_nh,"R_nh[m_nh]/D");
    MLInput->Branch("dE_nh",dE_nh,"dE_nh[m_nh]/D");
    MLInput->Branch("Epcal_nh",Epcal_nh,"Epcal_nh[m_nh]/D");
    MLInput->Branch("m2u_nh",m2u_nh,"m2u_nh[m_nh]/D");
    MLInput->Branch("m2v_nh",m2v_nh,"m2v_nh[m_nh]/D");
    
    MLInput->Branch("num_photons_0_1",&num_photons_0_1,"num_photons_0_1/D");
    MLInput->Branch("num_photons_0_2",&num_photons_0_2,"num_photons_0_2/D");
    MLInput->Branch("num_photons_0_35",&num_photons_0_35,"num_photons_0_35/D");
    
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
    
  while(_chain.Next()==true){
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
      // As of 5/16/2024 we are training the models without this cut
//       if (partstruct.pid == 2212 || partstruct.pid == -2212 ||
// 	  partstruct.pid == 2112 ||
// 	  partstruct.pid == -321 || partstruct.pid == -211 ||
// 	  partstruct.pid == 211 || partstruct.pid == 321) {
// 	if(partstruct.status>=4000 && partstruct.status<5000)
// 	  continue;
//       }
        
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
    
    // if(y>0.8) continue;
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
      
    //Loop over the particles in the event
      for (int ipart=0; ipart<Nmax; ++ipart) {
        //Check if the particle is a photon
        if (pid[ipart] == 22) {
          R_e = 0;
          dE_e = 0;
          //Initialize the arrays
          for (int i=0; i<m_g; ++i) {
            R_gamma[i] = 0;
            dE_gamma[i] = 0;
            Epcal_gamma[i] = 0;
            m2u_gamma[i] = 0;
            m2v_gamma[i] = 0;
          }
          for (int i=0; i<m_ch; ++i){
            R_ch[i] = 0;
            dE_ch[i] = 0;
            Epcal_ch[i] = 0;
            m2u_ch[i] = 0;
            m2v_ch[i] = 0;
          }
          for (int i=0; i<m_nh; ++i){
            R_nh[i] = 0;
            dE_nh[i] = 0;
            Epcal_nh[i] = 0;
            m2u_nh[i] = 0;
            m2v_nh[i] = 0;
          }
    
          //Set the number of photons within R<0.1, R<0.2, R<0.35
          num_photons_0_1 = 0;
          num_photons_0_2 = 0;
          num_photons_0_35 = 0;
        
          // Set vars
	  gE = E[ipart];
          gEpcal = pcal_e[ipart];
          gTheta = theta[ipart];
          gm2u = pcal_m2u[ipart];
          gm2v = pcal_m2v[ipart];
	  photon_has_match = (truepid[ipart]==22); // truepid[ipart]=-999 if no mc match found 
          //Loop over the particles in the event
          for (int jpart=0; jpart<Nmax; ++jpart) {

            //Skip the same particle
            if (ipart == jpart) continue;

            //Calculate R
            //Declare two TVector3 objects
            TVector3 v_1;
            TVector3 v_2;
            //Declare and set x, y, and z coordinates for each TVector3
            double x1,x2,y1,y2,z1,z2;
            //Check if the particle has values in pcal
            if(pcal_x[ipart]==-999){
                //Check if the particle has values in ecin
                if(ecin_x[ipart]==-999){
                    //If not, set coordinates for v_1 to values in ecout
                    x1=ecout_x[ipart]; 
                    y1=ecout_y[ipart]; 
                    z1=ecout_z[ipart];
                } 
                //Otherwise set coordinates for v_1 to values in ecin
                else{
                    x1=ecin_x[ipart]; 
                    y1=ecin_y[ipart]; 
                    z1=ecin_z[ipart];
                }
            } 
            //Otherwise set coordinates for v_1 to values in pcal
            else{
                x1=pcal_x[ipart]; 
                y1=pcal_y[ipart]; 
                z1=pcal_z[ipart];  
            }
            //Repeat for jpart
            if(pcal_x[jpart]==-999){
                if(ecin_x[jpart]==-999){
                    x2=ecout_x[jpart]; 
                    y2=ecout_y[jpart]; 
                    z2=ecout_z[jpart];
                } else{
                    x2=ecin_x[jpart]; 
                    y2=ecin_y[jpart]; 
                    z2=ecin_z[jpart];
                }
            } else{
                x2=pcal_x[jpart]; 
                y2=pcal_y[jpart]; 
                z2=pcal_z[jpart];  
            }
            //Set coordinates for each TVector3
            v_1.SetXYZ(x1,y1,z1);
            v_2.SetXYZ(x2,y2,z2);
              
            //Calculate the angle between the two TVector3 objects and store in R
            float R = v_1.Angle(v_2);
            
              
            if(pid[jpart]==22){
                //Count the number of photons within R<0.1, R<0.2, R<0.35
                if (R < 0.1) num_photons_0_1++;
                if (R < 0.2) num_photons_0_2++;
                if (R < 0.35) num_photons_0_35++;
                //Find the photon's nearest photon neighbors
                //Specifically, sort by smallest R value, of which there are "m_g" in the sorted list
                //The code is a bit longer because if a smaller R is found, we must bump forward all other elements in the list
                for (int i=0; i<m_g; ++i) {
                  if (R < R_gamma[i] || R_gamma[i] == 0) {
                    int j = m_g - 1;
                    while (j > i) {
                        R_gamma[j] = R_gamma[j - 1];
                        dE_gamma[j] = dE_gamma[j - 1];
                        Epcal_gamma[j] = Epcal_gamma[j - 1];
                        m2u_gamma[j] = m2u_gamma[j - 1];
                        m2v_gamma[j] = m2v_gamma[j - 1];
                        j--;
                    }
                    R_gamma[i] = R;
                    dE_gamma[i] = E[ipart] - E[jpart];
                    Epcal_gamma[i] = pcal_e[jpart];
                    m2u_gamma[i] = pcal_m2u[jpart];
                    m2v_gamma[i] = pcal_m2v[jpart];
                    break;
                  }
                }
            }
            else if(pid[jpart]==211 || pid[jpart]==-211 || pid[jpart]==2212 || pid[jpart]==-2212 || pid[jpart]==321 || pid[jpart]==-321){
                //Find the photon's nearest charged hadrons neighbors
                for (int i=0; i<m_ch; ++i) {
                  if (R < R_ch[i] || R_ch[i] == 0) {
                    int j = m_ch - 1;
                    while (j > i) {
                        R_ch[j] = R_ch[j - 1];
                        dE_ch[j] = dE_ch[j - 1];
                        Epcal_ch[j] = Epcal_ch[j - 1];
                        m2u_ch[j] = m2u_ch[j - 1];
                        m2v_ch[j] = m2v_ch[j - 1];
                        j--;
                    }
                    R_ch[i] = R;
                    dE_ch[i] = E[ipart] - E[jpart];
                    Epcal_ch[i] = pcal_e[jpart];
                    m2u_ch[i] = pcal_m2u[jpart];
                    m2v_ch[i] = pcal_m2v[jpart];
                    break;
                  }
                }
            }
            else if(pid[jpart]==2112 || pid[jpart]==-2112){
                //Find the photon's nearest neutral hadrons neighbors
                for (int i=0; i<m_nh; ++i) {
                  if (R < R_nh[i] || R_nh[i] == 0) {
                    int j = m_nh - 1;
                    while (j > i) {
                        R_nh[j] = R_nh[j - 1];
                        dE_nh[j] = dE_nh[j - 1];
                        Epcal_nh[j] = Epcal_nh[j - 1];
                        m2u_nh[j] = m2u_nh[j - 1];
                        m2v_nh[j] = m2v_nh[j - 1];
                        j--;
                    }
                    R_nh[i] = R;
                    dE_nh[i] = E[ipart] - E[jpart];
                    Epcal_nh[i] = pcal_e[jpart];
                    m2u_nh[i] = pcal_m2u[jpart];
                    m2v_nh[i] = pcal_m2v[jpart];
                    break;
                  }
                }
            }
            else if(pid[jpart]==11){
                if(R<R_e || R_e==0){
                    R_e = R;
                    dE_e = E[ipart] - E[jpart];
                }
            }
            else{
                continue;
            }
              
            
          }

          //Fill the MLInput TTree for each photon found
          MLInput->Fill();
        }
      }
    if (fill_EventTree)
        tree->Fill();
    tree_entries++;
  }
  
  outfile->cd();
  if (fill_EventTree)
      tree->Write();
  MLInput->Write();
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
