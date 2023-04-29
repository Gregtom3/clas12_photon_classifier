// Returns beam energy based on run
inline float runBeamEnergy(int run){
  if     (run>= 5032 && run<= 5666) return 10.6041; // rga fall 18
  else if(run>= 6616 && run<= 6783) return 10.1998; // rga spring 19
  else if(run>= 6120 && run<= 6399) return 10.5986; // rgb spring 19
  else if(run>= 6409 && run<= 6604) return 10.1998; // rgb spring 19
  else if(run>=11093 && run<=11283) return 10.4096; // rgb fall 19
  else if(run>=11284 && run<=11300) return 4.17179; // rgb fall BAND_FT 19
  else if(run>=11323 && run<=11571) return 10.3894; // rgb winter 20 (RCDB may still be incorrect)
  else if(run>=12000) return 10.559; // rgc 
  else if(run==11 || run==-11)                  return 10.6041; // MC for RGA inbending
  else {
    return 0.0;
  };
};