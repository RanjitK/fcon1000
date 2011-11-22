#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO DO SEGMENTATION OF ANTOMICAL SCAN
## parameters are passed from 0_preprocess.sh
##
## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
## for more information see www.nitrc.org/projects/fcon_1000
##
##########################################################################################################################

## subject
subject=0021002
## analysisdirectory/subject
dir=/Users/sikkas01/NYU/${subject}
## anatomical scan name (no extension)
anat=mprage
## anatomical scan directory
dir_anat=anat
## resting-state name (no extension)
rest=rest
## resting-state directory
dir_rest=func
##Where are your tissue priors? Should be ${scripts_dir}/tissuepriors/${standard_res}
prior_dir=/Users/sikkas01/scripts/tissuepriors/3mm
##set your desired spatial smoothing FWHM - it should match that used in 2_funcpreproc.sh
FWHM=6
sigma=`echo "scale=10 ; ${FWHM}/2.3548" | bc`



############################
## directory setup
anat_dir=${dir}/${dir_anat}
func_dir=${dir}/${dir_rest}

anat_reg_dir=${anat_dir}/reg
func_reg_dir=${func_dir}/reg

segment_dir=${dir}/segment
echo ${prior_dir}
PRIOR_WHITE=${prior_dir}/avg152T1_white_bin.nii.gz
PRIOR_CSF=${prior_dir}/avg152T1_csf_bin.nii.gz


##########################################################################################################################
##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
##########################################################################################################################

echo ------------------------------
echo !!!! RUNNING SEGMENTATION !!!!
echo ------------------------------

if [ -f ${segment_dir}/wm_mask ]; then

	echo ----------------------------------------------------------------------------------------
	echo 		!!! ${subject} segmentation already complete, SKIPPING !!!
	echo If you want to rerun segmentation, DELETE ${segment_dir} AND ${anat_dir}/segment_* FIRST
	echo ----------------------------------------------------------------------------------------
else

	## 1. Make segment dir
	mkdir ${segment_dir}

	## 2. Change to anat dir
	cwd=$( pwd )
	cd ${anat_dir}
	
	## 3. Segment the brain
	
	if [ -f ${anat_dir}/segment_prob_0.nii.gz ];
	then
		echo "Segmentation already complete for ${subject}; skipping"
	else
		echo "Segmenting brain for ${subject}"
		fast -t 1 -g -p -o segment ${anat}_brain.nii.gz
	fi

	## 4. Copy functional mask from FSLpreproc step 5 - this is the global signal mask
	echo "Creating global mask"
	if [ -f ${segment_dir}/global_mask.nii.gz ]; then rm ${segment_dir}/global_mask.nii.gz; else echo; fi
	3dcopy ${func_dir}/${rest}_pp_mask.nii.gz ${segment_dir}/global_mask.nii.gz
	
	## CSF
	## 5. Register csf to native space
	echo "Registering ${subject} csf to native (functional) space"
	flirt -in ${anat_dir}/segment_prob_0 -ref ${func_reg_dir}/example_func -applyxfm -init ${func_reg_dir}/highres2example_func.mat -out ${segment_dir}/csf2func
	## 6. Smooth image to match smoothing on functional
	echo "Smoothing ${subject} csf"
	fslmaths ${segment_dir}/csf2func -kernel gauss ${sigma} -fmean ${segment_dir}/csf_sm
	## 7. register to standard
	echo "Registering ${subject} csf to standard space"
	flirt -in ${segment_dir}/csf_sm -ref ${anat_reg_dir}/standard -applyxfm -init ${func_reg_dir}/example_func2standard.mat -out ${segment_dir}/csf2standard
	## 8. find overlap with prior
	echo "Finding overlap between ${subject} csf and prior"
	fslmaths ${segment_dir}/csf2standard -mas ${PRIOR_CSF} ${segment_dir}/csf_masked
	## 9. revert back to functional space
	echo "Registering ${subject} csf back to native space"
	flirt -in ${segment_dir}/csf_masked -ref ${func_reg_dir}/example_func -applyxfm -init ${func_reg_dir}/standard2example_func.mat -out ${segment_dir}/csf_native
	## 10. Threshold and binarize probability map of csf
	echo "Threshold and binarize ${subject} csf probability map"
	fslmaths ${segment_dir}/csf_native -thr 0.4 -bin ${segment_dir}/csf_bin
	## 11. Mask again by the subject's functional
	echo "Mask csf image by ${subject} functional"
	fslmaths ${segment_dir}/csf_bin -mas ${segment_dir}/global_mask ${segment_dir}/csf_mask
	
	## WM
	## 12. Register wm to native space
	echo "Registering ${subject} wm to native (functional) space"
	flirt -in ${anat_dir}/segment_prob_2 -ref ${func_reg_dir}/example_func -applyxfm -init ${func_reg_dir}/highres2example_func.mat -out ${segment_dir}/wm2func
	## 13. Smooth image to match smoothing on functional
	echo "Smoothing ${subject} wm"
	fslmaths ${segment_dir}/wm2func -kernel gauss ${sigma} -fmean ${segment_dir}/wm_sm
	## 14. register to standard
	echo "Registering ${subject} wm to standard space"
	flirt -in ${segment_dir}/wm_sm -ref ${anat_reg_dir}/standard -applyxfm -init ${func_reg_dir}/example_func2standard.mat -out ${segment_dir}/wm2standard
	## 15. find overlap with prior
	echo "Finding overlap between ${subject} wm and prior"
	fslmaths ${segment_dir}/wm2standard -mas ${PRIOR_WHITE} ${segment_dir}/wm_masked
	## 16. revert back to functional space
	echo "Registering ${subject} wm back to native space"
	flirt -in ${segment_dir}/wm_masked -ref ${func_reg_dir}/example_func -applyxfm -init ${func_reg_dir}/standard2example_func.mat -out ${segment_dir}/wm_native
	## 17. Threshold and binarize probability map of wm
	echo "Threshold and binarize ${subject} wm probability map"
	fslmaths ${segment_dir}/wm_native -thr 0.66 -bin ${segment_dir}/wm_bin
	## 18. Mask again by the subject's functional
	echo "Mask wm image by ${subject} functional"
	fslmaths ${segment_dir}/wm_bin -mas ${segment_dir}/global_mask ${segment_dir}/wm_mask
	
	##
	cd ${cwd}
fi
