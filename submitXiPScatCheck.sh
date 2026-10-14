#sleep 3600
ndiv=10
for i in {0..9}
do
  for frac in $(seq 1 "$ndiv")
  do
    #echo "bsub -q s root -b -q 'G4XiPScatAcceptanceCheck.cc($i,$frac,$ndiv)'"
    bsub -q s root -b -q 'G4XiPScatAcceptanceCheck.cc('$i','$frac','$ndiv')'
  done
done
#bsub -q a root -b -q G4XiAcceptanceCheckCH2.cc
