struct part{
  int pindex=-999;
  int status=-999;
  double px=-999;
  double py=-999;
  double pz=-999;
  double p=-999;
  double E=-999;
  double m=-999;
  double theta=-999;
  double phi=-999;
  int pid=-999;
  int truepid=-999; // Monte Carlo
  int trueparentpid=-999; // Monte Carlo
  int trueparentid=-999;  // Monte Carlo
  double vz=-999;
  double chi2=-999;
  double beta=-999;
  bool is_scattered_electron=false;

  // Calorimeter Info
  int    pcal_sector=-999;
  double pcal_e=-999;
  double pcal_x=-999;
  double pcal_y=-999;
  double pcal_z=-999;
  double pcal_lu=-999;
  double pcal_lv=-999;
  double pcal_lw=-999;
  double pcal_m2u=-999;
  double pcal_m2v=-999;
  double pcal_m2w=-999;
    
  int    ecin_sector=-999;
  double ecin_e=-999;
  double ecin_x=-999;
  double ecin_y=-999;
  double ecin_z=-999;
  double ecin_lu=-999;
  double ecin_lv=-999;
  double ecin_lw=-999;
  double ecin_m2u=-999;
  double ecin_m2v=-999;
  double ecin_m2w=-999;

  int    ecout_sector=-999;
  double ecout_e=-999;
  double ecout_x=-999;
  double ecout_y=-999;
  double ecout_z=-999;
  double ecout_lu=-999;
  double ecout_lv=-999;
  double ecout_lw=-999;
  double ecout_m2u=-999;
  double ecout_m2v=-999;
  double ecout_m2w=-999;
};

