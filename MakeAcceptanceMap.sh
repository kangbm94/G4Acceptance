bsub -q p -n 4 -R "rusage[mem=16000]" root -b -q G4XiAcceptance.cc 
