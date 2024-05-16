#include "CutManager.h"

// Constructors
CutManager::CutManager(){
  _run = 0;
}

// Public member functions
// Set the run number
void CutManager::set_run(int run){
  _run=run;
}
// Set the torus
void CutManager::set_torus(int torus){
  _torus=torus;
}

// Return a vector of particles that passes the cuts
std::vector<part> CutManager::filter_particles(std::vector<part> particles){
    
  std::vector<part> filtered_particles;
  // Get scattered electron
  part electron;
  for (auto particle : particles) {
    if(particle.is_scattered_electron==true) electron=particle;
  }
    
  for (auto particle : particles) {
    bool pass = false;
    int pid = particle.pid;
        
    switch(pid){
    case 11:
      pass = apply_electron_cuts(particle);
      break;
    case 22:
      pass = apply_photon_cuts(particle,electron);
      break;
    default:
      pass = true;
      break;
    }
    
    // Forward Detector Cut
    if(particle.theta*180/3.14159265<5 || particle.theta*180/3.14159265>35)
        pass = false;
    
    if (pass==true)
      filtered_particles.push_back(particle);
  }

  return filtered_particles;
}





// Protected member functions
// Apply all relevant electron cuts
bool CutManager::apply_electron_cuts(part particle){
  // if(VzCut(particle)==false) return false;
  // if(minEpcal(particle)==false) return false;
  // if(caloEdges(particle)==false) return false;
  return true;
}


// Apply all relevant photon cuts
bool CutManager::apply_photon_cuts(part particle, part electron){
  if(minEpcal(particle)==false) return false;
  if(photonMinEtot(particle)==false) return false;
  // if(photonElectronAngle(particle,electron)==false) return false;
  // if(photonBetaCut(particle)==false) return false;
  // if(caloEdges(particle)==false) return false;
  return true;
}




// Vz in target cell cut
bool CutManager::VzCut(part particle){
  if(particle.pid==11)
  { 
    if(particle.vz<-13||particle.vz>12) return false; 
  }
  return true;
}


// Minimum energy deposited in PCAL cut
bool CutManager::minEpcal(part particle){
  int pid=particle.pid;
  if(pid==11) return particle.pcal_e>0.07;
  else if(pid==22) return particle.pcal_e>0.00;
  return true;
}


// Fiducial cut on calorimeter edges
bool CutManager::caloEdges(part particle){
  int pid = particle.pid;
  double min_lu,max_lu,min_lv,max_lv,min_lw,max_lw;
  if(pid==11){
    min_lu=9;
    max_lu=400;
    min_lv=9;
    max_lv=400;
    min_lw=9;
    max_lw=400;
  }
  else if(pid==22){
    min_lu=14;
    max_lu=400;
    min_lv=14;
    max_lv=400;
    min_lw=14;
    max_lw=400;
  }
  else{
      return true;
  }
    
  if(particle.pcal_lv < min_lv || particle.pcal_lv > max_lv) return false;
  if(particle.pcal_lw < min_lw || particle.pcal_lw > max_lw) return false;
  return true;
}


bool CutManager::photonMinEtot(part particle){
  if(particle.E<0.2 && particle.pid==22) return false;
  return true;
}



bool CutManager::photonElectronAngle(part particle_A, part particle_B){

  TLorentzVector lv_A, lv_B;
  lv_A.SetPx(particle_A.px);
  lv_A.SetPy(particle_A.py);
  lv_A.SetPz(particle_A.pz);
  lv_A.SetE(particle_A.E);

  lv_B.SetPx(particle_B.px);
  lv_B.SetPy(particle_B.py);
  lv_B.SetPz(particle_B.pz);
  lv_B.SetE(particle_B.E);

  double angle = lv_A.Angle(lv_B.Vect());

  if (angle < 8 * TMath::DegToRad())
    return false;
  else
    return true;
}


bool CutManager::photonBetaCut(part particle){
  if(particle.beta<0.9||particle.beta>1.1) return false;
  return true;
}
