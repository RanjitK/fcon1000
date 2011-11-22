#!/usr/bin/env bash

##########################################################################################################################
## SCRIPT TO RUN GROUP-ANALYSES
##
## usage: ./groupAnalysis.sh <group_subject_list> <image_list> <model_list> <main_results_dir> 
## e.g. ./groupAnalysis.sh /home/me/scripts/my_group_subject_list.txt /home/me/scripts/my_group_image_list.txt /home/me/scripts/my_group_model_list /home/me/myresults
##
## Before running this script you have to set up and save a group model using FSL FEAT.
## For more information on how to do that, please refer to http://www.fmrib.ox.ac.uk/fsl/feat5/detail.html
## 
## WARNING: Keep you group_subject_list in the same order as you assinged subjects to your model!!!
## WARNING: Check your registrations to standard space before running a group analysis!
##
## Written by Maarten Mennes
## for more information see www.nitrc.org/projects/fcon_1000
##
##########################################################################################################################


## Required inputs

## 1. full/path/to/group_subject_list.txt containing the full path to the subjects you want to run
## The group_subject_list.txt should look as follows:
## /home/me/analysis/sub0001/func/RSFC
## /home/me/analysis/sub0002/func/RSFC
##
## e.g., group_subject_list=/fmri/RSFC/scripts/ADHD_controls.txt
##
group_subject_list=/fmri/RSFC/scripts/subjects.txt

## 2. full/path/to/image_list.txt containing the images you want to analyse
## The image_list.txt should look as follows:
## my_seed1_Z_2standard.nii.gz
## my_seed2_Z_2_standard.nii.gz
## ALL images must be located in the directory specified in the group_subject_list.txt
## The script will run the group analysis for each image
##
## e.g., image_list=/fmri/RSFC/scripts/group_seed_list.txt
##
image_list=/fmri/RSFC/scripts/seed_list.txt

## 3. full/path/to/model_list.txt containing the full path and the name of the model you want to run
## The model_list.txt should look as follows:
## /home/me/group_models/mymodel
## /home/me/group_models/mymodel_includingAge
## The script will run the group analysis for each model
##
## e.g., model_list=/fmri/RSFC/group_models/ADHD_vs_controls.txt
model_list=/fmri/RSFC/group_models/ADHD_vs_controls.txt

## 4. full/path/to/main_results_dir
## Enter the full path to the main directory where you want the results for each analysis located.
## The script will generate subdirectories for each analysis and put all results for that analysis in that directory
## e.g., /home/me/myresults
main_results_dir=/fmri/RSFC/group_results/ADHD_vs_controls

## 5. Additional parameters - Below are thresholding defaults. Do not alter unless you know what you're doing.
# 1. threshold parameters for easythresh
z_threshold=2.3
p_threshold=0.05


##########################################################################################################################
##---START OF SCRIPT----------------------------------------------------------------------------------------------------##
##########################################################################################################################
## Get included subjects
subjects=$( cat ${group_subject_list} )
## Get images to run
images=$( cat ${image_list} )
## Get models to run
models=$( cat ${model_list} )


## MODEL LOOP
for model in $models
do

	model_name=${model##*/}
	model_name=${model_name%*.fsf}
	model_dir=${model%*/*}

	## IMAGE LOOP
	for image in $images
	do
	
	current_model=${model_name}_${image%*.nii.gz}
	group_results_dir=${main_results_dir}/${current_model}
	# added to be able to rerun with the same inputs, if only one model or image fails
	# just comment if unwanted
	rm -r ${group_results_dir}
	
		echo
		echo
		echo ---------------------------------------------------
		echo "Running group analysis for ${current_model}"
		echo ---------------------------------------------------
				
		## check if all subjects have $image
		count=0
		for subject in $subjects
		do
			if [ ! -f $subject/$image ]
			then
			echo "$subject does not have $image, please correct"
			let "count +=1"
			fi
		done
		if [ $count -gt '0' ]
		then
			echo Not all subjects have the needed file, please correct.
			exit
		fi
	
	
		## create concatenated group file
		tmpfile=${main_results_dir}/temp.txt
		echo "Creating group file for ${current_model}"
		for subject in $subjects
		do
			echo ${subject}/${image} >> ${tmpfile}
		done	
		files=$( cat ${tmpfile} )
		fslmerge -t ${main_results_dir}/${current_model}_data.nii.gz $files
		rm ${tmpfile}
		## check whether data file was correctly made
		if [ ! -f ${main_results_dir}/${current_model}_data.nii.gz ]
		then
			echo data file for ${model_name} ${image} was not created, unable to complete analysis
			exit
		fi

		## create analysis specific mask
		if [ ! -f ${main_results_dir}/${model_name}_mask.nii.gz ]
		then
			echo making analysis specific mask
			fslmaths ${main_results_dir}/${current_model}_data.nii.gz -abs -Tmin -bin ${main_results_dir}/${current_model}_mask.nii.gz
		fi
		
		## check if all needed files are present
		for f in ${main_results_dir}/${current_model}_data.nii.gz ${main_results_dir}/${current_model}_mask.nii.gz ${model_dir}/${model_name}.mat ${model_dir}/${model_name}.con ${model_dir}/${model_name}.grp
		do
			if [ ! -f ${f} ]
			then
				echo ${f} was not found, please check
			exit
			fi
		done
		## check for ${model}.fts, but do not exit if absent as f-tests are optional
		if [ ! -f ${model_dir}/${model_name}.fts ]
		then
			echo no F-Tests will be calculated at this time
			ftest=no
		else
			echo F-Tests will be calculated
			ftest=yes
		fi
		
		## run flameo
		## without F-Test
		if [ ${ftest} == "no" ]
		then
			echo
			echo "running flameo for ${current_model}"
			flameo \
			--cope=${main_results_dir}/${current_model}_data.nii.gz \
			--mask=${main_results_dir}/${current_model}_mask.nii.gz \
			--dm=${model_dir}/${model_name}.mat \
			--tc=${model_dir}/${model_name}.con \
			--cs=${model_dir}/${model_name}.grp \
			--ld=${group_results_dir} \
			--runmode=ols
			
			# can't do this earlier because flameo create + directory if directory already exists
			mv ${main_results_dir}/${current_model}_data.nii.gz ${group_results_dir}/.
			mv ${main_results_dir}/${current_model}_mask.nii.gz ${group_results_dir}/.

		## with F-Test
		elif [ ${ftest} == "yes" ]
		then
			echo
			echo "running flameo for ${current_model}"
			flameo \
			--cope=${main_results_dir}/${current_model}_data.nii.gz \
			--mask=${main_results_dir}/${current_model}_mask.nii.gz \
			--dm=${model_dir}/${model_name}.mat \
			--tc=${model_dir}/${model_name}.con \
			--cs=${model_dir}/${model_name}.grp \
			--fc=${model_dir}/${model_name}.fts \
			--ld=${group_results_dir} \
			--runmode=ols

			# can't do this earlier because flameo create + directory if directory already exists
			mv ${main_results_dir}/${current_model}_data.nii.gz ${group_results_dir}/.
			mv ${main_results_dir}/${current_model}_mask.nii.gz ${group_results_dir}/.
		fi
		
		
		## Threshold zstat maps
		echo "Thresholding zstat-maps"
		curdir=$( pwd )
		cd ${group_results_dir}
		# get number of contrasts
		#cat ${model_dir}/${model_name}.con | grep '/NumContrasts' > temp3
		#numcon=`cat temp3 | cut -f 2`
		#rm temp3
		numcon=`ls -d zstat* 2> /dev/null | wc -l | sed 's/ //g'`
		echo "There are ${numcon} zstat maps to threshold"
		
		if [ ${numcon} -gt '0' ]
		then
			## Zstat Loop
			for ((z=1 ; z <= ${numcon} ; z++))
			do
				group_mm=` fslhd ${group_results_dir}/${current_model}_mask.nii.gz | grep pixdim3 | cut -f 9 -d ' ' | cut -f 1 -d'.'`
				echo "thresholding zstat${z}"
				easythresh zstat${z} \
				${group_results_dir}/${current_model}_mask.nii.gz \
				${z_threshold} ${p_threshold} \
				${FSLDIR}/data/standard/MNI152_T1_${group_mm}mm_brain.nii.gz \
				zstat${z}
			
				# threshold map at 0 to avoid negative values which show up for some reason in the thresholded maps
				fslmaths thresh_zstat${z}.nii.gz -thr 0 thresh_zstat${z}.nii.gz
			# end of zstat loop
			done
		fi
		
		## Threshold zfstat maps
		if [ ${ftest} == "yes" ]
		then
			echo "Thresholding zfstat-maps"
			curdir=$( pwd )
			cd ${group_results_dir}
			# get number of contrasts
			#cat ${model_dir}/${model_name}.con | grep '/NumContrasts' > temp3
			#numcon=`cat temp3 | cut -f 2`
			#rm temp3
			numcon=`ls -d zfstat* 2> /dev/null | wc -l | sed 's/ //g' `
			echo "There are ${numcon} zfstat maps to threshold"
			
			if [ ${numcon} -gt '0' ]
			then
				## Zstat Loop
				for ((z=1 ; z <= ${numcon} ; z++))
				do
					group_mm=` fslhd ${group_results_dir}/${current_model}_mask.nii.gz | grep pixdim3 | cut -f 9 -d ' ' | cut -f 1 -d '.'`
					echo "thresholding zfstat${z}"
					easythresh zfstat${z} \
					${group_results_dir}/${current_model}_mask.nii.gz \
					${z_threshold} ${p_threshold} \
					${FSLDIR}/data/standard/MNI152_T1_${group_mm}mm_brain.nii.gz \
					zfstat${z}
				
					# threshold map at 0 to avoid negative values which show up for some reason in the thresholded maps
					fslmaths thresh_zfstat${z}.nii.gz -thr 0 thresh_zfstat${z}.nii.gz
				# end of zfstat loop
				done
			fi
		# if F-test = yes
		fi

		## create group_reg file
		## this file can give you and idea of how well the subjects in your analysis overlay with each other and the MNI brain.
		## e.g., maybe there is one subject with limited coverage.
		echo Creating group_reg file
		echo "Use this command to check the group registration: fslview ${group_results_dir}/${current_model}_reg.nii.gz -l \"Random-Rainbow\""
		n_vols=` fslhd ${group_results_dir}/${current_model}_data.nii.gz | grep -w 'dim4' | sed 's/dim4 //g' `
		fslmaths ${group_results_dir}/${current_model}_data.nii.gz -abs -bin -Tmean -mul ${n_vols} ${group_results_dir}/${current_model}_reg.nii.gz

		## clean up
		cd ${curdir}

	## end of image loop
	done
## end of model loop
done




