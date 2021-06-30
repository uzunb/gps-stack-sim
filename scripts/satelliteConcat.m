function [satpAll, orbitparametersAll] = satelliteConcat(satpEph, satpAlm, orbitparametersEph, orbitparametersAlm)
  
  satpAll = [satpEph satpAlm];
  
  orbitparametersAll.svid  = [orbitparametersEph.svid  orbitparametersAlm.svid ];
  orbitparametersAll.A     = [orbitparametersEph.A     orbitparametersAlm.A ];
  orbitparametersAll.I     = [orbitparametersEph.I     orbitparametersAlm.I ];
  orbitparametersAll.Omega = [orbitparametersEph.Omega orbitparametersAlm.Omega];
  orbitparametersAll.e     = [orbitparametersEph.e     orbitparametersAlm.e ];
  orbitparametersAll.r     = [orbitparametersEph.r     orbitparametersAlm.r ];
  orbitparametersAll.u     = [orbitparametersEph.u     orbitparametersAlm.u ];
end