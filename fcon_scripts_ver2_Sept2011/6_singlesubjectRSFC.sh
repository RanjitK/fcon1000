#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO CALCULATE SEED-BASED RESTING-STATE FUNCTIONAL CONNECTIVITY
##
## This script can be run on its own, by filling in the appropriate parameters
## Alternatively this script gets called from batch_process.sh, where you can use it to run N sites, one after the other.
##
## Written by Clare Kelly, Maarten Mennes & Michael Milham
## for more information see www.nitrc.org/projects/fcon_1000
##
## Edited by CK 12 August 2011 to use 2mm MNI brain, nonlinear registration, and to fix holes in 3dfim+ output
##
##########################################################################################################################

## Directory where data are located (give full path)
## e.g., analysisdirectory=/fmri/RSFC/subjects
analysisdirectory=/Users/sikkas01/NYU/

## Name of subject list text file (give full path to file; file contains subject names/numbers only in a single column)
## e.g., subject_list=/fmri/RSFC/scripts/subs.txt
subject_list=/Users/sikkas01/NYU/NYU_subjects.txt

## Name of anatomical directory
## e.g., dir_anat=anat
dir_anat=session_1/anat_1

## Name of resting-state scan (DO NOT INCLUDE EXTENSION - i.e., do not include .nii.gz)
## e.g., rest_name=rest
rest_name=rest

## Name of resting-state directory
## e.g., rest_dir=rest_on
rest_dir=session_1/rest_1

## Name of seed list text file (give full path to file; file contains list with full/path/to/seed.nii.gz)
## e.g., seed_list=/fmri/RSFC/scripts/seeds.txt
seed_list=/Users/sikkas01/scripts/seed_list.txt

## What is the resolution of the standard brain used in registration (MNI brain)?
## e.g., standard_res=2mm
standard_res=3mm


##########################################################################################################################
##---START OF SCRIPT------------DO NOT EDIT BELOW THIS LINE-------------------------------------------------------------##
##########################################################################################################################


## Get subjects to run
subjects=$( cat ${subject_list} )
## Get seeds to run
seeds=$( cat ${seed_list} )

## A. SUBJECT LOOP
for subject in $subjects
do

	## directory setup
	func_dir=${analysisdirectory}/${subject}/${rest_dir}
	anat_dir=${analysisdirectory}/${subject}/${dir_anat}
	
	##registration
	anat_reg_dir=${anat_dir}/reg
	func_reg_dir=${func_dir}/reg

	## seed_timeseries
	seed_ts_dir=${func_dir}/seed_ts
	
	## RSFC maps
	RSFC_dir=${func_dir}/RSFC
	
	
	echo --------------------------
	echo running subject ${subject}
	echo --------------------------
	

	mkdir -p ${seed_ts_dir}
	mkdir -p ${RSFC_dir}
	
	
	## B. SEED_LOOP
	for seed in $seeds
	do
	
		seed_name=$( echo ${seed##*/} | sed s/\.nii\.gz//g )
		echo \------------------------
		echo running seed ${seed_name}
		echo \------------------------
	
		## Check if seed might have already been run
		if [ -f ${RSFC_dir}/${seed_name}_Z_2standard.nii.gz ]; then echo final file for seed ${seed_name} already exists; else echo; fi
	
		## 1. Extract Timeseries
		echo "Extracting timeseries for seed ${seed_name}"
		3dROIstats -quiet -mask_f2short -mask ${seed} ${func_dir}/${rest_name}_res2standard.nii.gz > ${seed_ts_dir}/${seed_name}.1D
	
		## 2. Compute voxel-wise correlation with Seed Timeseries		
		echo "Computing Correlation for seed ${seed_name}"
		3dfim+ -input ${func_dir}/${rest_name}_res.nii.gz -ideal_file ${seed_ts_dir}/${seed_name}.1D -fim_thr 0.0009 -out Correlation -bucket ${RSFC_dir}/${seed_name}_corr.nii.gz
	
		## 3. Z-transform correlations		
		echo "Z-transforming correlations for seed ${seed_name}"
		3dcalc -a ${RSFC_dir}/${seed_name}_corr.nii.gz -expr 'log((a+1)/(a-1))/2' -prefix ${RSFC_dir}/${seed_name}_Z.nii.gz
	
		## 4. Register Z-transformed correlations to standard space
		echo "Registering Z-transformed map to standard space (NONLINEAR)"
		#flirt -in ${RSFC_dir}/${seed_name}_Z.nii.gz -ref ${standard} -applyxfm -init ${reg_dir}/example_func2standard.mat -out ${RSFC_dir}/${seed_name}_Z_2standard.nii.gz
		applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_${standard_res}.nii.gz --in=${RSFC_dir}/${seed_name}_Z.nii.gz \
		--out=${RSFC_dir}/${seed_name}_Z_2standard.nii.gz --warp=${anat_reg_dir}/highres2standard_warp --premat=${func_reg_dir}/example_func2highres.mat
	
	## END OF SEED LOOP
	done

## END OF SUBJECT LOOP
done
