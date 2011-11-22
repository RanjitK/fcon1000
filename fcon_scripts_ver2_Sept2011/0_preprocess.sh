#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO RUN GENERAL RESTING-STATE PREPROCESSING
##
## This script can be run on its own, by filling in the appropriate parameters
## Alternatively this script gets called from batch_process.sh, where you can use it to run N sites, one after the other.
##
## Written by the Underpants Gnomes (a.k.a Clare Kelly, Zarrar Shehzad, Maarten Mennes & Michael Milham)
## for more information see www.nitrc.org/projects/fcon_1000
##
## Edited by CK 12 August 2011 to use 2mm MNI brain, nonlinear registration, and to include despiking
##
##########################################################################################################################

##########################################################################################################################
## PARAMETERS
## 
## Set each parameter (up to and including $9, TR) below
##
##########################################################################################################################
## Directory where these scripts are located (give full path)
## e.g., scripts_dir=/fmri/RSFC/scripts
scripts_dir=/Users/sikkas01/scripts/fcon_scripts_ver2_Sept2011

## Directory where data are located (give full path)
## e.g., analysisdirectory=/fmri/RSFC/subjects
analysisdirectory=/Users/sikkas01/NYU/

## Name of subject list text file (give full path to file; file contains subject names/numbers only in a single column)
## e.g., subject_list=/fmri/RSFC/scripts/subs.txt
subject_list=/Users/sikkas01/NYU/NYU_subjects.txt

## Name of anatomical scan (DO NOT INCLUDE EXTENSION - i.e., do not include .nii.gz)
## e.g., anat_name=mprage
anat_name=mprage

## Name of anatomical directory
## e.g., rest_dir=rest_on
anat_dir=session_1/anat_1

## Name of resting-state scan (DO NOT INCLUDE EXTENSION - i.e., do not include .nii.gz)
## e.g., rest_name=rest
rest_name=rest

## Name of resting-state directory
## e.g., rest_dir=rest_on
rest_dir=session_1/rest_1

## What is the first volume/time point of the resting state scan (default is 0; should be >0 if you need to disgard initial volumes)?
## e.g., first_vol=4
first_vol=0

## What is the last volmue/time point of the resting state scan (Should equal # of volumes minus 1)
## e.g., last_vol=175
last_vol=175

## What is the total number of volumes/time points in the resting state scan?
## e.g., n_vols=176
n_vols=176

## What is resting state scan TR (in seconds)?
## e.g., TR=2
TR=2

## What is the resolution of the standard brain (MNI brain)  used in registration?
## e.g., standard_res=2mm
standard_res=3mm

## What FWHM for spatial smoothing (should be 1.5-2 times the voxel size; default is 6mm for acquisition voxel size 3x3x4mm)?
FWHM=6
sigma=$( echo "scale=10;${FWHM}/2.3548" | bc )

## What low-pass (LP) and high-pass (HP) cut offs (default HP (i.e., keeping only frequencies higher than this point) = 0.01 ; 
## default LP (i.e., keeping only frequencies lower than this point) = 0.1)?
## YOU WILL NEED TO ALTER THESE (PARTICULARLY THE HP) BASED ON YOUR SCAN LENGTH, FOR LP OF 0.01, SCAN LENGTH MUST BE AT LEAST 100SECONDS
HP=0.01 ; LP=0.1


##########################################################################################################################


##########################################################################################################################
##---START OF SCRIPT----DO NOT EDIT BELOW THIS LINE---------------------------------------------------------------------##
##########################################################################################################################


## Get subjects to run
subjects=$( cat ${subject_list} )

## SUBJECT LOOP
for subject in $subjects
do

echo preprocessing $subject

##################################################################################################################################
## preprocess anatomical scan
## subject - full path to subject directory, subject is paremeterized - name of anatomical scan (no extension)
##################################################################################################################################
${scripts_dir}/1_anatpreproc.sh ${subject} ${analysisdirectory}/${subject} ${anat_name} ${anat_dir}

##################################################################################################################################
## preprocess functional scan
## subject - full path to subject directory, subject is parameterized - name of functional scan (no extension) -
## first volume - last volume - TR
## set temporal filtering and spatial smoothing in 2_funcpreproc.sh (default is 0.005-0.1Hz and 6 FWHM)
##################################################################################################################################
${scripts_dir}/2_funcpreproc.sh ${subject} ${analysisdirectory}/${subject} ${rest_name} ${first_vol} ${last_vol} ${TR} ${FWHM} ${HP} ${LP} ${rest_dir}

##################################################################################################################################
## registration
## subject - full path to subject directory, subject is parameterized - name of anatomical scan (no extension) - standard brain resolution
##################################################################################################################################
${scripts_dir}/3_registration.sh ${subject} ${analysisdirectory}/${subject} ${anat_name} ${anat_dir} ${standard_res} ${rest_dir}

##################################################################################################################################
## segmentation
## subject - full path to subject directory, subject is parameterized - name of anatomical scan (no extension) -
## name of resting-state scan (no extension) - resolution of tissuepriors (must match standard brain resolution)
##################################################################################################################################
${scripts_dir}/4_segment.sh ${subject} ${analysisdirectory}/${subject} ${anat_name} ${anat_dir} ${rest_name} ${rest_dir} ${scripts_dir}/tissuepriors/${standard_res} ${FWHM}

##################################################################################################################################
## nuisance signal regression
## subject - full path to subject directory, subject is parameterized - name of resting-state scan (no extension) -
## TR - number of volumes (default = last_vol + 1; has to be +1 since volume numbering starts at 0) -
## full path to dir where feat model template fsf file is located - standard brain resolution
##################################################################################################################################
${scripts_dir}/5_nuisance.sh ${subject} ${analysisdirectory}/${subject} ${rest_name} ${rest_dir} ${anat_dir} ${TR} ${n_vols} ${scripts_dir}/templates/nuisance.fsf ${standard_res}

## END OF SUBJECT LOOP
done
