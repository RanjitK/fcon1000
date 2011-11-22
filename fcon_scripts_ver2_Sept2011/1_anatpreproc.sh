#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO PREPROCESS THE ANATOMICAL SCAN
## parameters are passed from 0_preprocess.sh
##
## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
## for more information see www.nitrc.org/projects/fcon_1000
##
##########################################################################################################################

## subject
subject=$1
## analysisdirectory/subject
dir=$2
## name of the anatomical scan
anat=$3
##name of the anatomical directory
dir_anat=$4


###SETUP

anat_dir=${dir}/${dir_anat}

##########################################################################################################################
##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
##########################################################################################################################

echo --------------------------------------
echo !!!! PREPROCESSING ANATOMICAL SCAN!!!!
echo --------------------------------------

cwd=$( pwd )
cd ${anat_dir}

## 1. Deoblique anat
echo "deobliquing ${subject} anatomical"
3drefit -deoblique ${anat}.nii.gz

if [ -f ${anat}_RPI.nii.gz ]; then
	echo "${subject} anatomical already reoriented, SKIPPING"
else
	## 2. Reorient to fsl-friendly space
	echo "Reorienting ${subject} anatomical"
	3dresample -orient RPI -inset ${anat}.nii.gz -prefix ${anat}_RPI.nii.gz
fi

## 3. skull strip
if [ -f ${anat}_brain.nii.gz ]; then

	echo ----------------------------------------------------------------------------------
	echo 		!!! ${subject} skullstripping already complete, SKIPPING !!!
	echo If you want to rerun skullstripping, delete ${anat_dir}/${anat}_brain.nii.gz first
	echo ----------------------------------------------------------------------------------
else
	if [ "${anat}" == "mprage_skullstripped" ]; then cp ${anat}_RPI.nii.gz ${anat}_brain.nii.gz; exit; fi
	echo "skull stripping ${subject} anatomical"
	3dSkullStrip -input ${anat}_RPI.nii.gz -o_ply ${anat}_surf.nii.gz
	3dcalc -a ${anat}_RPI.nii.gz -b ${anat}_surf.nii.gz -expr 'a*step(b)' -prefix ${anat}_brain.nii.gz
fi
cd ${cwd}
