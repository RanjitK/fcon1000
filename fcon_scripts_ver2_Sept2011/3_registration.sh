#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO DO LINEAR AND NONLINEAR IMAGE REGISTRATION
## parameters are passed from 0_preprocess.sh
##
## !!!!!*****ALWAYS CHECK YOUR REGISTRATIONS*****!!!!!
##
##
## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
## for more information see www.nitrc.org/projects/fcon_1000
##
## Edited by CK 12 August 2011 to use 2mm MNI brain, nonlinear registration
##
##########################################################################################################################

## subject
subject=$1
## analysisdir/subject
dir=$2
## name of anatomical scan
anat=$3
## name of anatomical directory
dir_anat=$4
##standard brain to be used in registration
standard_res=$5
#standard=${FSLDIR}/data/standard/MNI152_T1_3mm_brain.nii.gz
dir_rest=$6

## directory setup
anat_dir=${dir}/${dir_anat}
func_dir=${dir}/${dir_rest}
###reg_dir=${dir}/reg
anat_reg_dir=${anat_dir}/reg
func_reg_dir=${func_dir}/reg


##########################################################################################################################
##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
##########################################################################################################################

echo ------------------------------
echo !!!! RUNNING REGISTRATION !!!!
echo ------------------------------

if [ -f ${func_reg_dir}/example_func2standard_NL.nii.gz ]; then

	echo -----------------------------------------------------------------------------------
	echo 		!!! ${subject} registration already complete, SKIPPING !!!
	echo If you want to rerun registration, DELETE ${anat_reg_dir} AND ${func_reg_dir} FIRST
	echo -----------------------------------------------------------------------------------
else	
	
	cwd=$( pwd )
	
	mkdir ${anat_reg_dir}
	mkdir ${func_reg_dir}
	
	## 1. Copy required images into reg directory
	### copy anatomical
	cp ${anat_dir}/${anat}_brain.nii.gz ${anat_reg_dir}/highres.nii.gz
	cp ${anat_dir}/${anat}_RPI.nii.gz ${anat_reg_dir}/highres_head.nii.gz
	### copy standard
	cp ${FSLDIR}/data/standard/MNI152_T1_${standard_res}_brain.nii.gz ${anat_reg_dir}/standard.nii.gz
	### copy example func created earlier
	cp ${func_dir}/example_func.nii.gz ${func_reg_dir}/
	
	## 2. cd into reg directory
	##cd ${reg_dir}
	
	## 3. FUNC->T1
	## You may want to change some of the options
	echo "${subject}: Registering functional to highres"
	flirt -ref ${anat_reg_dir}/highres -in ${func_reg_dir}/example_func -out ${func_reg_dir}/example_func2highres \
	-omat ${func_reg_dir}/example_func2highres.mat -cost corratio -dof 6 -interp trilinear
	# Create mat file for conversion from subject's anatomical to functional
	convert_xfm -inverse -omat ${func_reg_dir}/highres2example_func.mat ${func_reg_dir}/example_func2highres.mat
	
	## 4. T1->STANDARD
	## Initial linear registration
	echo "${subject}: Registering highres to standard (INITIAL LINEAR REGISTRATION)"
	if [ -f ${anat_reg_dir}/highres2standard.nii.gz ]; then
		echo "Linear highres to standard already complete, skipping"
	else
		flirt -ref ${FSLDIR}/data/standard/MNI152_T1_${standard_res}_brain.nii.gz -in ${anat_reg_dir}/highres \
		-out ${anat_reg_dir}/highres2standard -omat ${anat_reg_dir}/highres2standard.mat -cost corratio -searchcost corratio -dof 12 -interp trilinear
		## Create mat file for conversion from standard to high res
		convert_xfm -inverse -omat ${anat_reg_dir}/standard2highres.mat ${anat_reg_dir}/highres2standard.mat
	fi
	
	## 5. FUNC->STANDARD
	## Create mat file for registration of functional to standard
	echo "${subject}: Registering functional to standard (INITIAL LINEAR REGISTRATION)"
	convert_xfm -omat ${func_reg_dir}/example_func2standard.mat -concat ${anat_reg_dir}/highres2standard.mat \
	${func_reg_dir}/example_func2highres.mat
	## apply registration
	flirt -ref ${FSLDIR}/data/standard/MNI152_T1_${standard_res}_brain.nii.gz -in ${func_reg_dir}/example_func \
	-out ${func_reg_dir}/example_func2standard -applyxfm -init ${func_reg_dir}/example_func2standard.mat -interp trilinear
	## Create inverse mat file for registration of standard to functional
	convert_xfm -inverse -omat ${func_reg_dir}/standard2example_func.mat ${func_reg_dir}/example_func2standard.mat
	
	## 6. T1->STANDARD NONLINEAR
	# Perform nonlinear registration (higres to standard)
	echo "${subject}: Registering highres to standard (NONLINEAR). This may take a while."
	if [ -f ${anat_reg_dir}/highres2standard_NL.nii.gz ]; then
		echo "Nonlinear highres to standard already complete, skipping"
	else
		fnirt --in=${anat_reg_dir}/highres_head --aff=${anat_reg_dir}/highres2standard.mat \
		--cout=${anat_reg_dir}/highres2standard_warp --iout=${anat_reg_dir}/highres2standard_NL \
		--jout=${anat_reg_dir}/highres2standard_jac --config=T1_2_MNI152_${standard_res} --ref=${FSLDIR}/data/standard/MNI152_T1_${standard_res}.nii.gz \
		--refmask=${FSLDIR}/data/standard/MNI152_T1_${standard_res}_brain_mask_dil.nii.gz --warpres=10,10,10
	fi
	
	## 7. Apply nonlinear registration (func to standard)
	echo "${subject}: Registering functional to standard (NONLINEAR)"
	applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_${standard_res}.nii.gz --in=${func_reg_dir}/example_func --out=${func_reg_dir}/example_func2standard_NL \
	--warp=${anat_reg_dir}/highres2standard_warp --premat=${func_reg_dir}/example_func2highres.mat
	
	###***** ALWAYS CHECK YOUR REGISTRATIONS!!! YOU WILL EXPERIENCE PROBLEMS IF YOUR INPUT FILES ARE NOT ORIENTED CORRECTLY (IE. RPI, ACCORDING TO AFNI) *****###
	
	cd ${cwd}
	
fi
	
	
	
	
