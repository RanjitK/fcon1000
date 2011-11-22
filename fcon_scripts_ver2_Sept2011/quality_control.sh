#!/bin/bash

#########################################################################################################################
##SCRIPT TO RUN THE QUALITY CONTROL
## Fill in the parameters below
##
## Written by Clare Kelly
## for more information see www.nitrc.org/projects/fcon_1000
##
##########################################################################################################################

## Directory where data are located (give full path)
## e.g., analysisdirectory=/fmri/RSFC/subjects
analysisdirectory=${1}

## Name of resting-state scan (DO NOT INCLUDE EXTENSION - i.e., do not include .nii.gz)
## e.g., rest_name=rest
rest_name=${2}

## Name of resting-state directory
## e.g., rest_dir=restingstate
rest_dir=${3}

## Name of anatomical scan (DO NOT INCLUDE EXTENSION - i.e., do not include .nii.gz)
## e.g., anat_name=mprage
anat_name=anat

## Name of resting-state directory
## e.g., rest_dir=restingstate
anat_dir=${3}

## What is the resolution of the standard brain (MNI brain) used in registration?
## e.g., standard_res=2mm
standard_res=2mm

##What images below do you want to check? Enter YES or NO to each of the images listed:
##NOTE!!!! The 3_registration script MUST have been run in order to check registration!!!!!
## e.g., check_rest=YES
##       check_registration=NO
check_rest=YES
check_registration=NO


##########################################################################################################################


##########################################################################################################################
##---START OF SCRIPT----DO NOT EDIT BELOW THIS LINE---------------------------------------------------------------------##
##########################################################################################################################


cd ${analysisdirectory}

for subject in `ls -d *`
do

	if [ "${check_rest}" = "YES" ];
	then
		if [ ! -f ${subject}/${rest_dir}/${rest_name}.nii.gz ];
		then
			echo "PLEASE CHECK PATH SPECIFIED, ${subject}/${rest_dir}/${rest_name}.nii.gz DOES NOT EXIST"
		else
			echo "LOOKING AT ${subject} ${rest_name}"
			afni ${subject}/${rest_dir}/${rest_name}.nii.gz
		fi
	else
		echo
	fi

	if [ "${check_registration}" = "YES" ];
	then
		if [ ! -f ${subject}/${rest_dir}/${rest_name}.nii.gz ];
		then
			echo "PLEASE CHECK PATH SPECIFIED, ${subject}/${rest_dir}/${rest_name}.nii.gz DOES NOT EXIST"
		else
			echo "LOOKING AT ${subject} REST_ON REGISTRATION TO HIGHRES"	
			fslview ${subject}/${anat_dir}/${anat_name}_RPI.nii.gz  ${subject}/${rest_dir}/reg/example_func2highres.nii.gz 	

			echo "LOOKING AT ${subject} NORMALIZATION TO STANDARD"		
			fslview ${FSLDIR}/data/standard/MNI152_T1_${standard_res}.nii.gz  ${subject}/${anat_dir}/reg/highres2standard_NL.nii.gz -b 0,500
		fi
	else
		if [ ! -f ${subject}/${anat_dir}/${anat_name}.nii.gz ];
		then
			echo "PLEASE CHECK PATH SPECIFIED, ${subject}/${anat_dir}/${anat_name}.nii.gz DOES NOT EXIST"
		else
			echo "LOOKING AT ${subject} ${anat_name}"
			afni ${subject}/${anat_dir}/${anat_name}.nii.gz
		fi
	fi
	
done
