#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO PREPROCESS THE FUNCTIONAL SCAN
## parameters are passed from 0_preprocess.sh
##
## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
## for more information see www.nitrc.org/projects/fcon_1000
##
## Edited by CK 12 August 2011 to use 2mm MNI brain, nonlinear registration, and to include despiking
##
##########################################################################################################################

## subject
subject=$1
## analysisdirectory/subject
dir=$2
## name of the resting-state scan
rest=$3
## first timepoint (remember timepoint numbering starts from 0)
TRstart=$4
## last timepoint
TRend=$5
## TR
TR=$6
## Spatial smoothing FWHM
FWHM=$7
sigma=`echo "scale=10 ; ${FWHM}/2.3548" | bc`
## High pass and low pass cutoffs for temporal filtering
hp=$8
lp=$9
shift
## directory setup
func_dir=${dir}/$9

##########################################################################################################################
##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
##########################################################################################################################

echo ---------------------------------------
echo !!!! PREPROCESSING FUNCTIONAL SCAN !!!!
echo ---------------------------------------

cwd=$( pwd )
cd ${func_dir}

if [ -f ${rest}_pp.nii.gz ]; then

	echo -----------------------------------------------------------------------------------------
	echo 		!!! ${subject} initial preprocessing already complete, SKIPPING !!!
	echo If you want to rerun functional preprocessing, DELETE ${func_dir}/${rest}_pp.nii.gz FIRST
	echo -----------------------------------------------------------------------------------------
else	
	
	## 1. Dropping first # TRS
	echo "Dropping first TRs"
	if [ -f ${rest}_dr.nii.gz ]; then rm ${rest}_dr.nii.gz; else echo; fi
	3dcalc -a ${rest}.nii.gz[${TRstart}..${TRend}] -expr 'a' -prefix ${rest}_dr.nii.gz
	
	##2. Deoblique
	echo "Deobliquing ${subject}"
	3drefit -deoblique ${rest}_dr.nii.gz
	
	##3. Reorient into fsl friendly space (what AFNI calls RPI)
	echo "Reorienting ${subject}"
	if [ -f ${rest}_ro.nii.gz ]; then rm ${rest}_ro.nii.gz; else echo; fi
	3dresample -orient RPI -inset ${rest}_dr.nii.gz -prefix ${rest}_ro.nii.gz
	
	##4. Motion correct to average of timeseries
	echo "Motion correcting ${subject}"
	if [ -f ${rest}_mc.nii.gz ]; then rm ${rest}_mc.nii.gz ${rest}_ro_mean.nii.gz ${rest}_mc.1D ${rest}_maxdisp.1D; else echo; fi
	3dTstat -mean -prefix ${rest}_ro_mean.nii.gz ${rest}_ro.nii.gz 
	3dvolreg -Fourier -twopass -base ${rest}_ro_mean.nii.gz -zpad 4 -prefix ${rest}_mc.nii.gz -1Dfile ${rest}_mc.1D -maxdisp1D ${rest}_maxdisp.1D ${rest}_ro.nii.gz
	
	##5. Remove skull/edge detect
	echo "Skull stripping ${subject}"
	if [ -f ${rest}_ss.nii.gz ]; then rm ${rest}_ss.nii.gz ${rest}_mask.nii.gz; else echo; fi
	3dAutomask -prefix ${rest}_mask.nii.gz -dilate 1 ${rest}_mc.nii.gz
	3dcalc -a ${rest}_mc.nii.gz -b ${rest}_mask.nii.gz -expr 'a*b' -prefix ${rest}_ss.nii.gz
	
	##6. Get eighth image for use in registration
	echo "Getting example_func for registration for ${subject}"
	if [ -f example_func.nii.gz ]; then rm example_func.nii.gz; else echo; fi
	3dcalc -a ${rest}_ss.nii.gz[7] -expr 'a' -prefix example_func.nii.gz
	
	##6.5 (ADDED 12th August 2011) Despiking
	if [ -f ${rest}_ds.nii.gz ]; then rm ${rest}_ds.nii.gz; else echo; fi
	echo "Despiking ${subject}"
	3dDespike -prefix  ${rest}_ds.nii.gz ${rest}_ss.nii.gz
	
	##7. Spatial smoothing
	echo "Smoothing ${subject}"
	if [ -f ${rest}_sm.nii.gz ]; then rm ${rest}_sm.nii.gz; else echo; fi
	fslmaths ${rest}_ds.nii.gz -kernel gauss ${sigma} -fmean -mas ${rest}_mask.nii.gz ${rest}_sm.nii.gz
	
	##8. Grandmean scaling
	echo "Grand-mean scaling ${subject}"
	if [ -f ${rest}_gms.nii.gz ]; then rm ${rest}_gms.nii.gz; else echo; fi
	fslmaths ${rest}_sm.nii.gz -ing 10000 ${rest}_gms.nii.gz -odt float
	
	##9. Temporal filtering
	echo "Band-pass filtering ${subject}"
	if [ -f ${rest}_filt.nii.gz ]; then rm ${rest}_filt.nii.gz; else echo; fi
	3dFourier -lowpass ${lp} -highpass ${hp} -retrend -prefix ${rest}_filt.nii.gz ${rest}_gms.nii.gz
	
	##10.Detrending
	echo "Removing linear and quadratic trends for ${subject}"
	if [ -f ${rest}_dt.nii.gz ]; then rm ${rest}_filt_mean.nii.gz ${rest}_dt.nii.gz; else echo; fi
	3dTstat -mean -prefix ${rest}_filt_mean.nii.gz ${rest}_filt.nii.gz
	3dDetrend -polort 2 -prefix ${rest}_dt.nii.gz ${rest}_filt.nii.gz
	3dcalc -a ${rest}_filt_mean.nii.gz -b ${rest}_dt.nii.gz -expr 'a+b' -prefix ${rest}_pp.nii.gz
	
	##11.Create Mask
	echo "Generating mask of preprocessed data for ${subject}"
	if [ -f ${rest}_pp_mask.nii.gz ]; then rm ${rest}_pp_mask.nii.gz; else echo; fi
	fslmaths ${rest}_pp.nii.gz -Tmin -bin ${rest}_pp_mask.nii.gz -odt char
	
fi

cd ${cwd}
