#!/bin/bash
set -eu 

pdbs=`ls -1 ../*.pdb`
for pdb in $pdbs; do
  echo $pdb
  tmp=${pdb#../}
  prefix=${tmp%.pdb}
  pdbfixer $pdb \
      --output=${prefix}_solv_neutral.pdb \
      --water-box=2.5 2.5 2.5 \
      --positive-ion=Na+ \
      --negative-ion=Cl- 
      #--ionic-strength=0.15
done
