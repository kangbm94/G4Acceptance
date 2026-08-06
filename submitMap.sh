#sleep 3600
for i in {0..29}
do
  bsub -q s root -b -q 'G4XiAcceptance.cc('$i')'
  #for frac in {1..10}
  #do
    #bsub -q s root -b -q 'G4XiAcceptance.cc('$i','$frac',10)'
  #done
done
#bsub -q a root -b -q G4XiAcceptanceCheckCH2.cc
