#!/bin/bash
#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO CREATE .CSV FILE CONTAINING MAXIMUM RELATIVE MOTION PARAMETERS FOR EACH SUBJECT
##
## Fill in the parameters below
##
## Written by the Clare Kelly 12 August 2011
## for more information see Van Dijk, K.R.A., et al., The influence of head motion on intrinsic functional connectivity MRI, NeuroImage (2011), doi:10.1016/j.neuroimage.2011.07.044
##
##########################################################################################################################

##########################################################################################################################
## PARAMETERS
## 
## Set each parameter below
##
##########################################################################################################################


## Directory where these scripts are located (give full path)
## e.g., scripts_dir=/fmri/RSFC/scripts
scripts_dir=/Volumes/chong3/tim/RSFC/scripts

## Directory where data are located (give full path)
## e.g., analysisdirectory=/fmri/RSFC/subjects
analysisdirectory=/Volumes/chong3/tim/RSFC/subjects

## Name of subject list text file (give full path to file; file contains subject names/numbers only in a single column)
## e.g., subject_list=/fmri/RSFC/scripts/subs.txt
subject_list=/Volumes/chong3/tim/RSFC/scripts/subs.txt

## Name of resting-state scan (DO NOT INCLUDE EXTENSION - i.e., do not include .nii.gz)
## e.g., rest_name=rest
rest_name=rest

## Name of resting-state stem
## e.g., if you use rest_1, rest_2, rest_on, etc. the resting-state stem is rest
## e.g., if you use func, func_2, rest_on, etc. the resting-state stem is func etc.
restdir_stem=rest

## Name of nuisance directory (either 4_nuisance or nuisance)
## e.g., rest_dir=rest_on
nuisance_dir=nuisance

## Name of csv file to output movement parameters
## e.g., csv_name=subject_movement_15Aug2011.csv
csv_name=subject_movement_19Aug2011.csv





##########################################################################################################################


##########################################################################################################################
##---START OF SCRIPT----DO NOT EDIT BELOW THIS LINE---------------------------------------------------------------------##
##########################################################################################################################

## Set up .csv file

if [ -f ${csv_name} ]; then rm ${csv_name}; else echo; fi

#awk 'BEGIN {print "Subject,Mean SNR,std SNR,Max Relative Maxdisp,Max Abs Maxdisp,Max Relative Roll,Max Relative Pitch,Max Relative Yaw,Max Relative dS-I,Max Relative dL-R,Max #Relative dP-A,Max Abs Roll,Max Abs Pitch,Max Abs Yaw,Max Abs dS-I,Max Abs dL-R,Max Abs dP-A ";}'>${csvname}.csv
awk 'BEGIN {print "Subject,Rest Scan,Mean Relative RMS Displacement,Max Relative RSM Displacement,# Movements >0.1mm,Mean Relative Mean Rotation,Mean Relative Maxdisp,Max Relative Maxdisp,Max Abs Maxdisp,Max Relative Roll,Max Relative Pitch,Max Relative Yaw,Max Relative dS-I,Max Relative dL-R,Max Relative dP-A,Mean Relative Roll,Mean Relative Pitch,Mean Relative Yaw,Mean Relative dS-I,Mean Relative dL-R,Mean Relative dP-A,Max Abs Roll,Max Abs Pitch,Max Abs Yaw,Max Abs dS-I,Max Abs dL-R,Max Abs dP-A ";}'>${csv_name}

## Get subjects to run
subjects=$( cat ${subject_list} )

## SUBJECT LOOP
for subject in $subjects
do

	cd ${analysisdirectory}/${subject}
	
	for rest_dir in `ls -d ${restdir_stem}*`
	do

	## directory setup
	func_dir=${analysisdirectory}/${subject}/${rest_dir}

	cwd=$scripts_dir
	
	cd ${func_dir}
	echo ${func_dir}

	if [ -f movement_params.txt ]; then rm movement_params.txt; else echo; fi

	
	if [ -f ${rest_name}_mc.1D ]; then
			
		printf "%s", $subject >> movement_params.txt	
		printf "%s", $rest_dir >> movement_params.txt	
	
		##Relative RMS of translation (params 4-6)
		MEANrms=`awk '{x=$4} {y=$5} {z=$6} {printf( "%.6f\n", sqrt((x*x)+(y*y)+(z*z)))}' ${rest_name}_mc.1D |
		awk '{t=$1} NR>1{printf( "%.6f\n", t-p)}{p=t}' |
		awk '{ if($1>=0) { print $1} else {print $1*-1}}' |
		awk '{a += $1} END {print a/NR}'`
		printf "%.3f," ${MEANrms} >> movement_params.txt
		
		MAXrms=`awk '{x=$4} {y=$5} {z=$6} {printf( "%.6f\n", sqrt((x*x)+(y*y)+(z*z)))}' ${rest_name}_mc.1D |
		awk '{t=$1} NR>1{printf( "%.6f\n", t-p)}{p=t}' |
		awk '{ if($1>=0) { print $1} else {print $1*-1}}' |
		awk 'max=="" || $1 > max {max=$1} END{ print max}'`
		printf "%.3f," ${MAXrms} >> movement_params.txt
		
		##NUMBER OF relative RMS movements >0.1mm
		NUMmove=`awk '{x=$4} {y=$5} {z=$6} {printf( "%.6f\n", sqrt((x*x)+(y*y)+(z*z)))}' ${rest_name}_mc.1D |
		awk '{t=$1} NR>1{printf( "%.6f\n", t-p)}{p=t}' |
		awk '{ if($1>=0) { print $1} else {print $1*-1}}' |
		awk '{ if($1>=0.1) {a += 1}} END {print a}'`
		printf "%.3f," ${NUMmove} >> movement_params.txt

		##Mean of mean relative rotation (params 1-3)
		MEANrot=`awk '{ if($1>=0) {x=$1} else {x=$1*-1}} {printf( "%.6f\t", x )} { if($2>=0) {y=$2} else {y=$2*-1}} {printf( "%.6f\t", y)} { if($3>=0) {z=$3} else {z=$3*-1}} {printf( "%.6f\n", z)}' ${rest_name}_mc.1D |
		awk '{ mean=($1+$2+$3)/3} {printf( "%.6f\n", mean)}'|
		awk '{t=$1} NR>1{printf( "%.6f\n", t-p)}{p=t}'| 
		awk '{ if($1>=0) { print $1} else {print $1*-1}}' | 
		awk '{a += $1} END {print a/NR}'`
		printf "%.3f," ${MEANrot} >> movement_params.txt
	
	else
	echo "SUBJECT ${subject} does not have a ${rest_name}_mc.1D file in ${func_dir}"
	fi
	
	if [ -f ${rest_name}_maxdisp.1D ]; then
	
		#snrMEAN=`cat 0_checks/1_func/snr.1D | cut -d " " -f 1`
		#snrSTD=`cat 0_checks/1_func/snr.1D | cut -d " " -f 2`
		#printf "%.3f," ${snrMEAN} ${snrSTD}>> movement_params.txt		
	
		mean=`awk '{t=$1} NR>1{printf( "%.3f\n", t-p)}{p=t}' ${rest_name}_maxdisp.1D |awk '{a += $1} END {print a/NR}'`
		printf "%.3f," ${mean} >> movement_params.txt		
		
		num1=`awk '{t=$1} NR>1{printf( "%.3f\n", t-p)}{p=t}' ${rest_name}_maxdisp.1D |awk 'max=="" || $1 > max {max=$1}; END {print max}'`
		num2=`awk '{t=$1} NR>1{printf( "%.3f\n", t-p)}{p=t}' ${rest_name}_maxdisp.1D |awk 'min=="" || $1 < min {min=$1}; END {print min}' |awk '{ if($1>=0) { print $1} else {print $1*-1}}'`
		max=`awk 'BEGIN{if ('$num1'>'$num2') { print '$num1'} else {print '$num2'}}'`
		printf "%.3f," ${max} >> movement_params.txt
		
		max=`awk 'max=="" || $1 > max {max=$1} END{ print max}' ${rest_name}_maxdisp.1D`
		printf "%.3f," ${max} >> movement_params.txt	
	else
	echo "SUBJECT ${subject} does not have a ${rest_name}_maxdisp.1D file in ${func_dir}"
	fi
	
	
	if [ -f ${nuisance_dir}/mc1.1D ]; then		
		for j in {1..6}
		do
		
			num1=`awk -F: '{t=$1} NR>1{printf( "%.6f\n", t-p)}{p=t}' ${nuisance_dir}/mc${j}.1D |awk 'max=="" || $1 > max {max=$1}; END {print max}'`
			num2=`awk -F: '{t=$1} NR>1{printf( "%.6f\n", t-p)}{p=t}' ${nuisance_dir}/mc${j}.1D |awk 'min=="" || $1 < min {min=$1}; END {print min}' |awk '{ if($1>=0) { print $1} else {print $1*-1}}'`
			max=`awk 'BEGIN{if ('$num1'>'$num2') { print '$num1'} else {print '$num2'}}'`
			printf "%.6f," ${max} >> movement_params.txt

		done
		
		for j in {1..6}
		do	       
			mean=`awk -F: '{t=$1} NR>1{printf( "%.6f\n", t-p)}{p=t}' ${nuisance_dir}/mc${j}.1D |awk '{a += $1} END {print a/NR}'`
			printf "%.6f," ${mean} >> movement_params.txt

		done
			
		
		for j in {1..6}
		do	       
		
			num1=`awk -F: 'max=="" || $1 > max {max=$1}; END {print max}' ${nuisance_dir}/mc${j}.1D`
       	       		num2=`awk -F: 'min=="" || $1 < min {min=$1}; END {print min}' ${nuisance_dir}/mc${j}.1D |awk '{ if($1>=0) { print $1} else {print $1*-1}}'`
			max=`awk 'BEGIN{if ('$num1'>'$num2') { print '$num1'} else {print '$num2'}}'`
			printf "%.6f," ${max} >> movement_params.txt

		done
		
	else
	echo "SUBJECT ${subject} does not have motion parameters in ${func_dir}/nuisance"
	fi
	
	
	if [ -f movement_params.txt ]; then
	awk '{print $0}' movement_params.txt >>${scripts_dir}/${csv_name}		
	else
	echo "UNABLE TO GENERATE MOTION PARAMS FOR SUBJECT ${subject} ${rest_dir}"
	fi
	
	
	cd $cwd
	
	done
done
	

	
	
	
	
	
