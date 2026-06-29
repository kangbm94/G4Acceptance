bsub -q p -n 2 -R "rusage[mem=8000]" root -b -q G4XiAcceptanceCheck.cc
