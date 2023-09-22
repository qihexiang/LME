#!/bin/bash

find -name "*.mol2" -exec obabel {} -O {}.gjf \;

for gjf in $(find -name "*.gjf");
do
sed -i "s@#@%chk="${gjf%.*}".chk\n%nprocshared=10\n%mem=4GB\n# opt freq wb97xd/genecp scrf=(smd,solvent=toluene)@g" $gjf
cat << EOF >> $gjf
C H O P N Si 0
6-31G*
****
Mn 0
lanl2dz
****

Mn 0
lanl2dz

EOF
done