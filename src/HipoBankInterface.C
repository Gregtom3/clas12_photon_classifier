#include "HipoBankInterface.h"

HipoBankInterface::HipoBankInterface(){}

HipoBankInterface::HipoBankInterface(const std::unique_ptr<clas12::clas12reader>& _c12){
  
  // Add REC::Cal info
  _idx_RECCal = _c12->addBank("REC::Calorimeter");
  _ipindex_RECCal = _c12->getBankOrder(_idx_RECCal,"pindex");
  _ix_RECCal  = _c12->getBankOrder(_idx_RECCal,"x");
  _iy_RECCal  = _c12->getBankOrder(_idx_RECCal,"y");
  _iz_RECCal  = _c12->getBankOrder(_idx_RECCal,"z");
  _ilu_RECCal = _c12->getBankOrder(_idx_RECCal,"lu");
  _ilv_RECCal = _c12->getBankOrder(_idx_RECCal,"lv");
  _ilw_RECCal = _c12->getBankOrder(_idx_RECCal,"lw");
  _im2u_RECCal = _c12->getBankOrder(_idx_RECCal,"m2u");
  _im2v_RECCal = _c12->getBankOrder(_idx_RECCal,"m2v");
  _im2w_RECCal = _c12->getBankOrder(_idx_RECCal,"m2w");
  _ilayer_RECCal = _c12->getBankOrder(_idx_RECCal,"layer");
  _isector_RECCal = _c12->getBankOrder(_idx_RECCal,"sector");
  _itime_RECCal = _c12->getBankOrder(_idx_RECCal,"time");
  _ipath_RECCal = _c12->getBankOrder(_idx_RECCal,"path");
  _ienergy_RECCal = _c12->getBankOrder(_idx_RECCal,"energy");
}





bool HipoBankInterface::loadBankData(const std::unique_ptr<clas12::clas12reader>& _c12 , part &particle){
  clear();

  // Grab necessary particle info
  // -------------------------------------------------------------
  int pindex = particle.pindex;

  // -------------------------------------------------------------
  // Parse the REC::Calorimeter
  // -------------------------------------------------------------

  for(auto i = 0 ; i < _c12->getBank(_idx_RECCal)->getRows() ; i++){
    // Continue loop if the pindex in the calo bank does not match
    if(_c12->getBank(_idx_RECCal)->getInt(_ipindex_RECCal,i)!=pindex)
      continue;
    
    int sectorCal = _c12->getBank(_idx_RECCal)->getInt(_isector_RECCal,i);
    float timeCal = _c12->getBank(_idx_RECCal)->getFloat(_itime_RECCal,i);
    float pathCal = _c12->getBank(_idx_RECCal)->getFloat(_ipath_RECCal,i);
    float lu = _c12->getBank(_idx_RECCal)->getFloat(_ilu_RECCal,i);
    float lv = _c12->getBank(_idx_RECCal)->getFloat(_ilv_RECCal,i);
    float lw = _c12->getBank(_idx_RECCal)->getFloat(_ilw_RECCal,i);
    float m2u = _c12->getBank(_idx_RECCal)->getFloat(_im2u_RECCal,i);
    float m2v = _c12->getBank(_idx_RECCal)->getFloat(_im2v_RECCal,i);
    float m2w = _c12->getBank(_idx_RECCal)->getFloat(_im2w_RECCal,i);
    float x = _c12->getBank(_idx_RECCal)->getFloat(_ix_RECCal,i);
    float y = _c12->getBank(_idx_RECCal)->getFloat(_iy_RECCal,i);
    float z = _c12->getBank(_idx_RECCal)->getFloat(_iz_RECCal,i);
    int layerCal = _c12->getBank(_idx_RECCal)->getInt(_ilayer_RECCal,i);
    int calidx = -1; //Array index for lu, lv, lw

    switch(layerCal){
    case 1: //PCal
      calidx = 0;
      _Ele_PCAL_e = _c12->getBank(_idx_RECCal)->getFloat(_ienergy_RECCal,i);
      break;
    case 4: //ECIN
      calidx = 1;
      _Ele_ECIN_e = _c12->getBank(_idx_RECCal)->getFloat(_ienergy_RECCal,i);
      break;
    case 7: //ECOUT
      calidx = 2;
      _Ele_ECOUT_e = _c12->getBank(_idx_RECCal)->getFloat(_ienergy_RECCal,i);
      break;   
    }
    
    // If there was a pindex attached to one of the calo layers...
    if(calidx!=-1){
      _x_Cal[calidx]=x;
      _y_Cal[calidx]=y;
      _z_Cal[calidx]=z;
    
      _lu_Cal[calidx]=lu;
      _lv_Cal[calidx]=lv;
      _lw_Cal[calidx]=lw;

      _m2u_Cal[calidx]=m2u;
      _m2v_Cal[calidx]=m2v;
      _m2w_Cal[calidx]=m2w;
        
      _sector_Cal[calidx]=sectorCal;
      _time_Cal[calidx]=timeCal;
      _path_Cal[calidx]=pathCal;
    
    }
  }

  importDataToParticle(particle);
  
  return true;
}

bool HipoBankInterface::importDataToParticle(part &particle)
{
  
  // -------------------------------------------------------------
  // Import the REC::Calorimeter data
  // -------------------------------------------------------------
  particle.pcal_sector = _sector_Cal[0];
  particle.ecin_sector = _sector_Cal[1];
  particle.ecout_sector = _sector_Cal[2];

  particle.pcal_e = _Ele_PCAL_e;
  particle.ecin_e = _Ele_ECIN_e;
  particle.ecout_e = _Ele_ECOUT_e;
  
  particle.pcal_m2u = _m2u_Cal[0];
  particle.pcal_m2v = _m2v_Cal[0];
  particle.pcal_m2w = _m2w_Cal[0];
    
  particle.pcal_x = _x_Cal[0];
  particle.ecin_x = _x_Cal[1];
  particle.ecout_x = _x_Cal[2];

  particle.pcal_y = _y_Cal[0];
  particle.ecin_y = _y_Cal[1];
  particle.ecout_y = _y_Cal[2];

  particle.pcal_z = _z_Cal[0];
  particle.ecin_z = _z_Cal[1];
  particle.ecout_z = _z_Cal[2];
  particle.pcal_lu = _lu_Cal[0];
  particle.ecin_lu = _lu_Cal[1];
  particle.ecout_lu = _lu_Cal[2];

  particle.pcal_lv = _lv_Cal[0];
  particle.ecin_lv = _lv_Cal[1];
  particle.ecout_lv = _lv_Cal[2];

  particle.pcal_lw = _lw_Cal[0];
  particle.ecin_lw = _lw_Cal[1];
  particle.ecout_lw = _lw_Cal[2];
  return true;
}

void HipoBankInterface::clear(){
  _Ele_PCAL_e = 0.0;
  _Ele_ECIN_e = 0.0;
  _Ele_ECOUT_e = 0.0;
  for(int i = 0 ; i < 3 ; i++){
    _sector_Cal[i]=0;
    _time_Cal[i]=0;
    _path_Cal[i]=0;
    _x_Cal[i]=0;
    _y_Cal[i]=0;
    _z_Cal[i]=0;
    _m2u_Cal[i]=0;
    _m2v_Cal[i]=0;
    _m2w_Cal[i]=0;
    _lu_Cal[i]=0;
    _lv_Cal[i]=0;
    _lw_Cal[i]=0;
  }
}
