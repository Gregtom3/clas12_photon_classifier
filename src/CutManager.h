#ifndef CutManager_h
#define CutManager_h
class CutManager{
 public:
  // Constructors
  CutManager();
  // Public member variables and functions
  void set_run(int);
  void set_torus(int);
  std::vector<part> filter_particles(std::vector<part>);

 protected:
  // Protected member functions
  bool apply_electron_cuts(part);
  bool apply_photon_cuts(part,part);
    
  bool VzCut(part);
  bool minEpcal(part);
  bool caloEdges(part);
  bool photonMinEtot(part);
  bool photonElectronAngle(part,part);
  bool photonBetaCut(part);
    
  private:
        
  // Private member variables
  int _run=0;
  int _torus=0;      
};
#endif
