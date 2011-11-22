#!/bin/bash

###
### You DO NOT NEED TO EDIT this script to create seeds.
###
### You simply need to create a seedlist.txt file that contains 
### the following information for each seed you want to create:
###
###  x y z name radius resolution outputdirectory 
###
### Where:
### x, y and z are the coordinates IN MNI SPACE of the center of the ROI
### name is the name you want to call the output file
### radius is the radius of the sphere (in mm)
### resolution is the resolution (e.g., 2mm) of the standard brain used in preprocessing
### outputdirectory is the name (with full path) of where you want the seeds to be saved
###
###
### For example, a seedlist.txt file should look as follows:
###	
###  27 -21 -12 righthippo 4 2mm /fmri/RSFC/seeds
### -27 -21 -12 lefthippo 4 2mm /fmri/RSFC/seeds
###
### Then, to run the script, simply type: ./make_spherical_seeds.sh seedlist.txt
###
###




############-------------------------------------############
############-------------------------------------############
############-------------------------------------############
############!!!!!DO NOT EDIT BELOW THIS LINE!!!!!############
############-------------------------------------############
############-------------------------------------############
############-------------------------------------############



while read line
do

###Enter x, y and z coords (center of sphere)
x_coord=$( echo $line | cut -d ' ' -f 1 )
y_coord=$( echo $line | cut -d ' ' -f 2 )
z_coord=$( echo $line | cut -d ' ' -f 3 )

###What should output (seed) be called?
seed_name=$( echo $line | cut -d ' ' -f 4 )

###What radius will be seed be (typical = 4mm) NOTE that this is in mm. 
###Do NOT append "mm" to the number entered here
seed_radius=$( echo $line | cut -d ' ' -f 5 )

###Where should output (seed) be saved - what directory?
output_dir=$( echo $line | cut -d ' ' -f 7 )

###What resolution will be seed be (e.g., 1mm, 2mm, 3mm?) 
###NOTE that this MUST match the resolution of your single-subject standard-space residuals image
###i.e.,  the resolution of your "standard brain"
seed_resolution=$( echo $line | cut -d ' ' -f 6 )

###Where (in what directory) is your standard brain stored? e.g., /home/fmri/scripts/templates
standard_dir=${FSLDIR}/data/standard/



echo ${x_coord} ${y_coord} ${z_coord} |3dUndump -prefix ${output_dir}/${seed_name}.nii -master ${standard_dir}/MNI152_T1_${seed_resolution}_brain.nii.gz -srad ${seed_radius} -orient LPI -xyz -

done <${1}
