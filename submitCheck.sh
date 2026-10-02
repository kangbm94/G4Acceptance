#sleep 3600
ndiv=5
for i in {0..29}
do
  for frac in $(seq 1 "$ndiv")
  do
    echo "bsub -q s root -b -q 'G4XiAcceptanceCheck.cc($i,$frac,$ndiv)'"
    bsub -q s root -b -q 'G4XiAcceptanceCheck.cc('$i','$frac','$ndiv')'
  done
done
#bsub -q a root -b -q G4XiAcceptanceCheckCH2.cc
