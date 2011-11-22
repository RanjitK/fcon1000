#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO DO REGRESS OUT NUISANCE COVARIATES FROM RESTING_STATE SCAN
## nuisance covariates are: global signal, white matter (WM, WHITE), CSF, and
## 6 motion parameters obtained during motion correction step (see 2_funcpreproc.sh)
##
## parameters are passed from 0_preprocess.sh
##
## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
## for more information see www.nitrc.org/projects/fcon_1000
##
## Edited by CK 12 August 2011 to use 2mm MNI brain, nonlinear registration,
##
##########################################################################################################################
## subject
subject=0021002
## analysisdirectory/subject
dir=/Users/sikkas01/NYU/${subject}
## resting-state filename (no extension)
rest=rest
## resting-state directory
dir_rest=func
## anatomical scan directory
dir_anat=anat
## TR
TR=2.0
## number of timepoints in the resting-state scan
n_vols=180
## full path to template nuisance feat .fsf file; e.g. /full/path/to/template.fsf
nuisance_template=/Users/sikkas01/scripts/templates/nuisance.fsf
##Resolution of the standard brain (MNI brain) 
standard_res=3mm

## directory setup
anat_dir=${dir}/${dir_anat}
func_dir=${dir}/${dir_rest}

anat_reg_dir=${anat_dir}/reg
func_reg_dir=${func_dir}/reg

segment_dir=${dir}/segment
nuisance_dir=${func_dir}/nuisance


##########################################################################################################################
##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
##########################################################################################################################

echo --------------------------------------------
echo !!!! RUNNING NUISANCE SIGNAL REGRESSION !!!!
echo --------------------------------------------

if [ -f ${func_dir}/${rest}_res2standard.nii.gz ]; then

	echo ---------------------------------------------------------------------------------------------------------
	echo 			!!! ${subject} nuisance regression already complete, SKIPPING !!!
	echo If you want to rerun nuisance regression, DELETE ${nuisance_dir} AND ${func_dir}/${rest}_res.nii.gz FIRST
	echo ---------------------------------------------------------------------------------------------------------
else

	## 1. make nuisance directory
	mkdir -p ${nuisance_dir}

	# 2. Seperate motion parameters into seperate files
	echo "Splitting up ${subject} motion parameters"
	awk '{print $1}' ${func_dir}/${rest}_mc.1D > ${nuisance_dir}/mc1.1D
	awk '{print $2}' ${func_dir}/${rest}_mc.1D > ${nuisance_dir}/mc2.1D
	awk '{print $3}' ${func_dir}/${rest}_mc.1D > ${nuisance_dir}/mc3.1D
	awk '{print $4}' ${func_dir}/${rest}_mc.1D > ${nuisance_dir}/mc4.1D
	awk '{print $5}' ${func_dir}/${rest}_mc.1D > ${nuisance_dir}/mc5.1D
	awk '{print $6}' ${func_dir}/${rest}_mc.1D > ${nuisance_dir}/mc6.1D
	
	# Extract signal for global, csf, and wm
	## 3. Global
	echo "Extracting global signal for ${subject}"
	3dmaskave -mask ${segment_dir}/global_mask.nii.gz -quiet ${func_dir}/${rest}_pp.nii.gz > ${nuisance_dir}/global.1D

	## 4. csf
	echo "Extracting signal from csf for ${subject}"
	3dmaskave -mask ${segment_dir}/csf_mask.nii.gz -quiet ${func_dir}/${rest}_pp.nii.gz > ${nuisance_dir}/csf.1D
		
	## 5. wm
	echo "Extracting signal from white matter for ${subject}"
	3dmaskave -mask ${segment_dir}/wm_mask.nii.gz -quiet ${func_dir}/${rest}_pp.nii.gz > ${nuisance_dir}/wm.1D
	
	## 6. Generate mat file (for use later)
	## create fsf file
	echo "Modifying model file"
	sed -e s:nuisance_dir:"${nuisance_dir}":g <${nuisance_template} >${nuisance_dir}/temp1
	sed -e s:nuisance_model_outputdir:"${nuisance_dir}/residuals.feat":g <${nuisance_dir}/temp1 >${nuisance_dir}/temp2
	sed -e s:nuisance_model_TR:"${TR}":g <${nuisance_dir}/temp2 >${nuisance_dir}/temp3
	sed -e s:nuisance_model_numTRs:"${n_vols}":g <${nuisance_dir}/temp3 >${nuisance_dir}/temp4
	sed -e s:nuisance_model_input_data:"${func_dir}/${rest}_pp.nii.gz":g <${nuisance_dir}/temp4 >${nuisance_dir}/nuisance.fsf 
	
	rm ${nuisance_dir}/temp*
	
	echo "Running feat model"
	feat_model ${nuisance_dir}/nuisance
	
	minVal=`3dBrickStat -min -mask ${func_dir}/${rest}_pp_mask.nii.gz ${func_dir}/${rest}_pp.nii.gz`
	
	## 7. Get residuals
	echo "Running film to get residuals"
	if [ -d ${nuisance_dir}/stats ]; then rm -r ${nuisance_dir}/stats; else echo; fi
	film_gls -rn ${nuisance_dir}/stats -noest -sa -ms 5 ${func_dir}/${rest}_pp.nii.gz ${nuisance_dir}/nuisance.mat ${minVal}
	
	## 8. Demeaning residuals and ADDING 100
	3dTstat -mean -prefix ${nuisance_dir}/stats/res4d_mean.nii.gz ${nuisance_dir}/stats/res4d.nii.gz
	3dcalc -a ${nuisance_dir}/stats/res4d.nii.gz -b ${nuisance_dir}/stats/res4d_mean.nii.gz -expr '(a-b)+100' -prefix ${func_dir}/${rest}_res.nii.gz

	## 9. Resampling residuals to MNI space
	#flirt -ref ${reg_dir}/standard -in ${func_dir}/${rest}_res -out ${func_dir}/${rest}_res2standard -applyxfm -init ${reg_dir}/example_func2standard.mat -interp trilinear
	if [ -f ${func_dir}/${rest}_res2standard.nii.gz ]; then rm ${func_dir}/${rest}_res2standard.nii.gz; else echo; fi
	echo "Resampling residuals to ${standard_res} MNI space"	
	applywarp --ref=${FSLDIR}/data/standard/MNI152_T1_${standard_res}.nii.gz --in=${func_dir}/${rest}_res --out=${func_dir}/${rest}_res2standard \
	--warp=${anat_reg_dir}/highres2standard_warp --premat=${func_reg_dir}/example_func2highres.mat
	
fi
	
	
	
	
	
	
	
	
	
