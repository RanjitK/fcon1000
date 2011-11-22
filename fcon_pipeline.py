#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import e_afni
import fcon_nodes
import sys
import os
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
from nipype.interfaces.afni.preprocess import ThreedSkullStrip
from nipype.interfaces.afni.preprocess import Threedcalc
from nipype.interfaces.afni.base import Info, AFNITraitedSpec, AFNICommand
from nipype.interfaces.base import Bunch, TraitedSpec, File, Directory, InputMultiPath
from nipype.utils.filemanip import fname_presuffix, list_to_filename, split_filename
from nipype.utils.docparse import get_doc
from nipype.interfaces.fsl import MultiImageMaths
from nipype.interfaces.fsl import ExtractROI
from nipype.interfaces.fsl.maths import DilateImage
from multiprocessing import Process
from multiprocessing import Pool
from ConfigParser import SafeConfigParser

#*******************************************************************************************
FSLDIR = ''
scripts_dir = ''
prior_dir = ''
anat_name = ''
anat_dir = ''
rest_name = ''
rest_dir = ''
working_dir = ''
seed_list= ''
standard_brain = ''
batch_list = ''
standard_res = ''
logFile = ''
nuisance_template = ''
recon_subjects = ''
seeds_dir = ''
infosource = None
datasource = None
datasource_seeds = None
datasink = None
workflow = None
HP = None
LP = None
hp = None
lp = None
#*******************************************************************************************



def anatpreproc():

	global anat_name
	global infosource
	global datasource
	global workflow
	
	#workflow.connect( infosource,'subject_id',datasource, 'subject_id' )
	workflow.connect(datasource, 'anat',fcon_nodes.anat_refit,'in_file')
	workflow.connect(fcon_nodes.anat_refit,'out_file',fcon_nodes.anat_reorient,'in_file')
	workflow.connect(fcon_nodes.anat_reorient,'out_file',fcon_nodes.anat_skullstrip,'in_file')
	workflow.connect(fcon_nodes.anat_skullstrip,'out_file',fcon_nodes.anat_calc,'infile_b')
	workflow.connect(fcon_nodes.anat_reorient,'out_file',fcon_nodes.anat_calc,'infile_a')
	
def last_vol(vols):

	v = []
	for vol in vols:
		v.append(int(vol) - 1)

	return v

def funcpreproc():

	global rest_name
	global infosource
	global datasource
	global workflow


	workflow.connect(datasource, 'rest',fcon_nodes.TR,'in_files')
	workflow.connect(datasource, 'rest',fcon_nodes.NVOLS,'in_files')
	workflow.connect(fcon_nodes.NVOLS,('nvols',last_vol),fcon_nodes.func_calc,'stop_idx')
	workflow.connect(datasource, 'rest', fcon_nodes.func_calc,'infile_a')
	workflow.connect( fcon_nodes.func_calc , 'out_file', fcon_nodes.func_refit , 'in_file' )
	workflow.connect( fcon_nodes.func_refit , 'out_file' , fcon_nodes.func_reorient , 'in_file' )
	workflow.connect( fcon_nodes.func_reorient , 'out_file' , fcon_nodes.func_tstat , 'in_file' )
	workflow.connect( fcon_nodes.func_reorient , 'out_file' , fcon_nodes.func_volreg , 'in_file' )
	workflow.connect( fcon_nodes.func_tstat , 'out_file' , fcon_nodes.func_volreg , 'basefile' )
	workflow.connect( fcon_nodes.func_volreg , 'out_file', fcon_nodes.func_automask , 'in_file')
	workflow.connect( fcon_nodes.func_volreg,'out_file', fcon_nodes.func_calcR,'infile_a')
	workflow.connect( fcon_nodes.func_automask,'out_file', fcon_nodes.func_calcR,'infile_b')
	workflow.connect( fcon_nodes.func_calcR,'out_file',fcon_nodes.func_calcI,'infile_a')
	workflow.connect( fcon_nodes.func_calcR, 'out_file', fcon_nodes.func_despike,'in_file')
	workflow.connect( fcon_nodes.func_despike, 'out_file', fcon_nodes.func_smooth , 'in_file')
	workflow.connect( fcon_nodes.func_automask,'out_file',fcon_nodes.func_smooth,'operand_files')
	workflow.connect( fcon_nodes.func_smooth,'out_file' , fcon_nodes.func_scale , 'in_file')
	workflow.connect( fcon_nodes.func_scale, 'out_file', fcon_nodes.func_filter, 'in_file')
	workflow.connect( infosource,'hp',fcon_nodes.func_filter,'highpass')
	workflow.connect( infosource,'lp',fcon_nodes.func_filter,'lowpass')
	workflow.connect( fcon_nodes.func_filter, 'out_file', fcon_nodes.func_detrenda , 'in_file')
	workflow.connect( fcon_nodes.func_filter, 'out_file', fcon_nodes.func_detrendb, 'in_file')
	workflow.connect( fcon_nodes.func_detrenda, 'out_file' , fcon_nodes.func_detrendc , 'infile_a')
	workflow.connect( fcon_nodes.func_detrendb , 'out_file', fcon_nodes.func_detrendc , 'infile_b')
	workflow.connect( fcon_nodes.func_detrendc, 'out_file' , fcon_nodes.func_mask, 'in_file')



def registration():

	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global workflow
	global infosource

	
	workflow.connect(fcon_nodes.func_calcI, 'out_file', fcon_nodes.reg_flirt, 'in_file')
	workflow.connect(fcon_nodes.anat_calc,'out_file', fcon_nodes.reg_flirt, 'reference')
	workflow.connect(fcon_nodes.reg_flirt, 'out_matrix_file', fcon_nodes.reg_xfm1, 'in_file')
	workflow.connect(fcon_nodes.anat_calc,'out_file', fcon_nodes.reg_flirt1, 'in_file')
	workflow.connect(infosource,'standard_res_brain', fcon_nodes.reg_flirt1, 'reference')
	workflow.connect(fcon_nodes.reg_flirt1,'out_matrix_file',fcon_nodes.reg_xfm2, 'in_file')
	workflow.connect(fcon_nodes.reg_flirt,'out_matrix_file',fcon_nodes.reg_xfm3, 'in_file')
	workflow.connect(fcon_nodes.reg_flirt1, 'out_matrix_file',fcon_nodes.reg_xfm3, 'in_file2')
	workflow.connect(fcon_nodes.reg_xfm3,'out_file',fcon_nodes.reg_flirt2,'in_matrix_file')
	workflow.connect(infosource,'standard_res_brain', fcon_nodes.reg_flirt2, 'reference')
	workflow.connect(fcon_nodes.func_calcI, 'out_file', fcon_nodes.reg_flirt2, 'in_file')
	workflow.connect(fcon_nodes.reg_xfm3, 'out_file', fcon_nodes.reg_xfm4, 'in_file')
	workflow.connect(fcon_nodes.anat_reorient, 'out_file', fcon_nodes.reg_fnt , 'in_file')
	workflow.connect(fcon_nodes.reg_flirt1,'out_matrix_file', fcon_nodes.reg_fnt, 'affine_file')
	workflow.connect(infosource,'standard',fcon_nodes.reg_fnt,'ref_file')
	workflow.connect(infosource,'standard_brain_mask_dil',fcon_nodes.reg_fnt,'refmask_file' )
	workflow.connect(infosource,'config_file',fcon_nodes.reg_fnt,'config_file')
	workflow.connect(fcon_nodes.func_calcI, 'out_file',fcon_nodes.reg_warp,'in_file')
	workflow.connect(infosource,'standard',fcon_nodes.reg_warp, 'ref_file')
	workflow.connect(fcon_nodes.reg_fnt, 'fieldcoeff_file', fcon_nodes.reg_warp,'field_file')
	workflow.connect(fcon_nodes.reg_flirt,'out_matrix_file', fcon_nodes.reg_warp,'premat')

def pick_wm_0(probability_maps):

	import sys
	import os
	if(isinstance(probability_maps,list)):
		for file in probability_maps:
			print file
			if file.endswith("prob_0.nii.gz"):
				
				return file
	return None


def pick_wm_1(probability_maps):
	
	import sys
	if(isinstance(probability_maps,list)):
		for file in probability_maps:
			print file
			if file.endswith("prob_2.nii.gz"):
			
				return file
	return None

def segment():
	
	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource


	workflow.connect(fcon_nodes.anat_calc,'out_file',fcon_nodes.seg_segment, 'in_files' )
	workflow.connect( fcon_nodes.seg_segment, ('probability_maps',pick_wm_0), fcon_nodes.seg_flirt , 'in_file' )
	workflow.connect(fcon_nodes.func_calcI,'out_file',fcon_nodes.seg_flirt,'reference')
	workflow.connect(fcon_nodes.reg_xfm1, 'out_file', fcon_nodes.seg_flirt,'in_matrix_file')
	workflow.connect( fcon_nodes.seg_flirt , 'out_file' , fcon_nodes.seg_smooth , 'in_file' )
	workflow.connect( fcon_nodes.seg_smooth, 'out_file', fcon_nodes.seg_flirt1, 'in_file')
	workflow.connect(infosource,'standard_res_brain', fcon_nodes.seg_flirt1, 'reference')
	workflow.connect(fcon_nodes.reg_xfm3, 'out_file',fcon_nodes.seg_flirt1, 'in_matrix_file')
	workflow.connect( fcon_nodes.seg_flirt1, 'out_file', fcon_nodes.seg_smooth1, 'in_file')
	workflow.connect( infosource, 'PRIOR_CSF', fcon_nodes.seg_smooth1, 'operand_files')
	workflow.connect( fcon_nodes.seg_smooth1,'out_file', fcon_nodes.seg_flirt2,'in_file')
	workflow.connect( fcon_nodes.func_calcI,'out_file', fcon_nodes.seg_flirt2,'reference')
	workflow.connect( fcon_nodes.reg_xfm4, 'out_file', fcon_nodes.seg_flirt2,'in_matrix_file')
	workflow.connect( fcon_nodes.seg_flirt2, 'out_file', fcon_nodes.seg_thresh, 'in_file')
	workflow.connect( fcon_nodes.seg_thresh,'out_file', fcon_nodes.seg_mask,'in_file')
	workflow.connect(fcon_nodes.func_mask, 'out_file' ,fcon_nodes.seg_copy, 'in_file' )
	workflow.connect( fcon_nodes.seg_copy,'out_file',  fcon_nodes.seg_mask,'operand_files')
	workflow.connect(  fcon_nodes.seg_segment, ('probability_maps',pick_wm_1),  fcon_nodes.seg_flirt3 , 'in_file' )
	workflow.connect( fcon_nodes.func_calcI,'out_file', fcon_nodes.seg_flirt3 , 'reference' )
	workflow.connect( fcon_nodes.reg_xfm1, 'out_file', fcon_nodes.seg_flirt3 , 'in_matrix_file' )
	workflow.connect( fcon_nodes.seg_flirt3, 'out_file',  fcon_nodes.seg_smooth2, 'in_file')
	workflow.connect( fcon_nodes.seg_smooth2,'out_file', fcon_nodes.seg_flirt4, 'in_file')
	workflow.connect( infosource,'standard_res_brain', fcon_nodes.seg_flirt4, 'reference')
	workflow.connect( fcon_nodes.reg_xfm3, 'out_file', fcon_nodes.seg_flirt4, 'in_matrix_file')
	workflow.connect( fcon_nodes.seg_flirt4,'out_file', fcon_nodes.seg_prior1, 'in_file')
	workflow.connect( infosource, 'PRIOR_WHITE', fcon_nodes.seg_prior1, 'operand_files')
	workflow.connect( fcon_nodes.seg_prior1,'out_file', fcon_nodes.seg_flirt5, 'in_file')
	workflow.connect( fcon_nodes.func_calcI,'out_file', fcon_nodes.seg_flirt5, 'reference')
	workflow.connect( fcon_nodes.reg_xfm4,'out_file', fcon_nodes.seg_flirt5, 'in_matrix_file')
	workflow.connect( fcon_nodes.seg_flirt5,'out_file', fcon_nodes.seg_thresh1,'in_file')
	workflow.connect( fcon_nodes.seg_thresh1,'out_file', fcon_nodes.seg_mask1,'in_file')
	workflow.connect( fcon_nodes.seg_copy,'out_file', fcon_nodes.seg_mask1,'operand_files')


def signalExtraction():

	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource

	workflow.connect(fcon_nodes.func_detrendc,'out_file',fcon_nodes.nuisance_globalE,'in_file')
	workflow.connect(fcon_nodes.seg_copy,'out_file',fcon_nodes.nuisance_globalE,'mask')
	workflow.connect(fcon_nodes.func_detrendc,'out_file',fcon_nodes.nuisance_csf,'in_file')
	workflow.connect(fcon_nodes.seg_mask,'out_file',fcon_nodes.nuisance_csf,'mask')
	workflow.connect(fcon_nodes.func_detrendc,'out_file',fcon_nodes.nuisance_wm,'in_file')
	workflow.connect(fcon_nodes.seg_mask1,'out_file',fcon_nodes.nuisance_wm,'mask')
	workflow.connect(fcon_nodes.func_volreg , 'oned_file',fcon_nodes.MC,'in_file')
	workflow.connect(fcon_nodes.func_detrendc, 'out_file',fcon_nodes.MC,'pp')
	workflow.connect(infosource, 'nuisance_template',fcon_nodes.FSF,'nuisance_template')
	workflow.connect(fcon_nodes.func_detrendc, 'out_file',fcon_nodes.FSF,'rest_pp')
	workflow.connect(fcon_nodes.TR,'TR', fcon_nodes.FSF,'TR')
	workflow.connect(fcon_nodes.NVOLS,'nvols',fcon_nodes.FSF,'n_vols')
	
	#workflow.connect(fcon_nodes.MC,'EV_Lists',fcon_nodes.EV_List,'EV_Lists')
	#workflow.connect(fcon_nodes.nuisance_globalE,'out_file',fcon_nodes.EV_List,'global1Ds')
	#workflow.connect(fcon_nodes.nuisance_csf,'out_file',fcon_nodes.EV_List,'csf1Ds')
	#workflow.connect(fcon_nodes.nuisance_wm,'out_file',fcon_nodes.EV_List,'wm1Ds')
	workflow.connect(fcon_nodes.MC,'EV_Lists',fcon_nodes.copyS,'EV_Lists')
	workflow.connect(fcon_nodes.FSF,'nuisance_files',fcon_nodes.copyS,'nuisance_files')
	workflow.connect(fcon_nodes.nuisance_globalE,'out_file',fcon_nodes.copyS,'global1Ds')
	workflow.connect(fcon_nodes.nuisance_csf,'out_file',fcon_nodes.copyS,'csf1Ds')
	workflow.connect(fcon_nodes.nuisance_wm,'out_file',fcon_nodes.copyS,'wm1Ds')

def nuisance():

	

	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource


	signalExtraction()
	
	workflow.connect(fcon_nodes.FSF,'nuisance_files',fcon_nodes.nuisance_featM,'fsf_file')
	workflow.connect(fcon_nodes.copyS,'EV_final_lists_set',fcon_nodes.nuisance_featM,'ev_files')
	workflow.connect(fcon_nodes.func_detrendc,'out_file',fcon_nodes.nuisance_brick,'in_file')
	workflow.connect(fcon_nodes.func_mask, 'out_file' ,fcon_nodes.nuisance_brick,'mask')
	workflow.connect(fcon_nodes.nuisance_featM, 'design_file',fcon_nodes.nuisance_fgls,'design_file')
	workflow.connect(fcon_nodes.func_detrendc, ('out_file',fcon_nodes.getStatsDir),fcon_nodes.nuisance_fgls, 'results_dir' )
	workflow.connect(fcon_nodes.func_detrendc,'out_file',fcon_nodes.nuisance_fgls,'in_file')
	workflow.connect(fcon_nodes.nuisance_brick,'min_val',fcon_nodes.nuisance_fgls,'threshold')
	workflow.connect(fcon_nodes.nuisance_fgls,'residual4d',fcon_nodes.nuisance_stat,'in_file')
	workflow.connect(fcon_nodes.nuisance_fgls,'residual4d',fcon_nodes.nuisance_calc,'infile_a')
	workflow.connect(fcon_nodes.nuisance_stat,'out_file', fcon_nodes.nuisance_calc,'infile_b')
	workflow.connect(fcon_nodes.nuisance_calc,'out_file', fcon_nodes.nuisance_warp,'in_file')
	workflow.connect(infosource,'standard',fcon_nodes.nuisance_warp,'ref_file')
	workflow.connect(fcon_nodes.reg_fnt, 'fieldcoeff_file',fcon_nodes.nuisance_warp,'field_file')
	workflow.connect(fcon_nodes.reg_flirt,'out_matrix_file',fcon_nodes.nuisance_warp,'premat')


def despkingWF():
	
	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource


	workflow.connect(fcon_nodes.func_despike, 'out_file',fcon_nodes.alff_detrend,'in_file')
	workflow.connect(fcon_nodes.alff_detrend,'out_file',fcon_nodes.alff_smooth,'in_file')
	workflow.connect(fcon_nodes.func_automask,'out_file',fcon_nodes.alff_smooth,'operand_files')
	workflow.connect(fcon_nodes.alff_smooth,'out_file',fcon_nodes.alff_scale,'in_file')
	workflow.connect(fcon_nodes.alff_scale,'out_file',fcon_nodes.alff_cp,'in_file')
	workflow.connect(fcon_nodes.alff_cp,'out_file',fcon_nodes.alff_mean,'in_file')

tolist = lambda x: [x]

def takemod(nvols):

	print '>>>>' + str(nvols)
	decisions = []
	for vol in nvols:
		mod = int (int(vol) % 2)

		if mod == 1:
			decisions.append([0])
			#return [0]
		else:
			decisions.append([1])
			#return [1]

	return decisions
	

def set_op_str(n2):

	strs = []
	for n in n2:	
		str = "-Tmean -mul %f" %(n)
		strs.append(str)
	return strs


def set_op1_str(nvols):

	strs = []
	for vol in nvols:
		str = '-Tmean -mul %d -div 2' %(int(vol))
		strs.append(str)
	
	return strs

def FALFF_NORM_REGISTER_WF():

	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource


	workflow.connect(fcon_nodes.alff_sqrt,'out_file',fcon_nodes.alff_falff,'in_file')
	workflow.connect(fcon_nodes.NVOLS,('nvols',set_op1_str),fcon_nodes.alff_falff,'op_string')
	workflow.connect(fcon_nodes.alff_sum,'out_file',fcon_nodes.alff_falff1,'in_file')
	workflow.connect(fcon_nodes.alff_falff,'out_file',fcon_nodes.alff_falff1,'operand_files')

	workflow.connect(fcon_nodes.alff_sum,'out_file',fcon_nodes.alff_normM,'in_file')
	workflow.connect(fcon_nodes.func_automask,'out_file',fcon_nodes.alff_normM,'mask_file')
	workflow.connect(fcon_nodes.alff_sum,'out_file',fcon_nodes.alff_normS,'in_file')
	workflow.connect(fcon_nodes.func_automask,'out_file',fcon_nodes.alff_normS,'mask_file')
	workflow.connect(fcon_nodes.alff_falff1,'out_file',fcon_nodes.alff_normM1,'in_file')
	workflow.connect(fcon_nodes.func_automask,'out_file',fcon_nodes.alff_normM1,'mask_file')
	workflow.connect(fcon_nodes.alff_falff1,'out_file',fcon_nodes.alff_normS1,'in_file')
	workflow.connect(fcon_nodes.func_automask,'out_file',fcon_nodes.alff_normS1,'mask_file')

	workflow.connect(fcon_nodes.alff_normM,'out_stat', fcon_nodes.alff_op_string,'mean')
	workflow.connect(fcon_nodes.alff_normS,'out_stat', fcon_nodes.alff_op_string,'std_dev')
	workflow.connect(fcon_nodes.alff_op_string,'op_string', fcon_nodes.alff_Z_alff, 'op_string')
	workflow.connect(fcon_nodes.alff_sum,'out_file',fcon_nodes.alff_Z_alff,'in_file')
	workflow.connect(fcon_nodes.func_automask,'out_file',fcon_nodes.alff_Z_alff,'operand_files')

	workflow.connect(fcon_nodes.alff_normM1,'out_stat', fcon_nodes.alff_op_string1,'mean')
	workflow.connect(fcon_nodes.alff_normS1,'out_stat', fcon_nodes.alff_op_string1,'std_dev')
	workflow.connect(fcon_nodes.alff_op_string1,'op_string', fcon_nodes.alff_Z_falff, 'op_string')
	workflow.connect(fcon_nodes.alff_falff1,'out_file',fcon_nodes.alff_Z_falff,'in_file')
	workflow.connect(fcon_nodes.func_automask,'out_file',fcon_nodes.alff_Z_falff,'operand_files')

	workflow.connect(infosource,'standard', fcon_nodes.alff_warp_alff,'ref_file')
	workflow.connect(fcon_nodes.alff_Z_alff,'out_file',fcon_nodes.alff_warp_alff,'in_file')
	workflow.connect(fcon_nodes.reg_fnt, 'fieldcoeff_file',fcon_nodes.alff_warp_alff,'field_file')
	workflow.connect(fcon_nodes.reg_flirt,'out_matrix_file',fcon_nodes.alff_warp_alff,'premat')

	workflow.connect(infosource,'standard', fcon_nodes.alff_warp_falff,'ref_file')
	workflow.connect(fcon_nodes.alff_Z_falff,'out_file',fcon_nodes.alff_warp_falff,'in_file')
	workflow.connect(fcon_nodes.reg_fnt, 'fieldcoeff_file',fcon_nodes.alff_warp_falff,'field_file')
	workflow.connect(fcon_nodes.reg_flirt,'out_matrix_file',fcon_nodes.alff_warp_falff,'premat')



def falff():

	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource

	despkingWF()

	workflow.connect(fcon_nodes.alff_cp,'out_file',fcon_nodes.alff_roi,'in_file')
	workflow.connect(fcon_nodes.NVOLS,'nvols',fcon_nodes.alff_roi,'t_size')
	workflow.connect(fcon_nodes.alff_cp,'out_file',fcon_nodes.alff_cp1,'in_file')

	workflow.connect(fcon_nodes.alff_roi,'roi_file',fcon_nodes.alff_concatnode,'in1')
	workflow.connect(fcon_nodes.alff_cp1,'out_file',fcon_nodes.alff_concatnode,'in2')
	workflow.connect(fcon_nodes.alff_concatnode, 'out', fcon_nodes.alff_selectnode, 'inlist')
	workflow.connect(fcon_nodes.NVOLS,('nvols',takemod), fcon_nodes.alff_selectnode, 'index')
	workflow.connect(fcon_nodes.alff_selectnode,'out', fcon_nodes.alff_pspec,'in_file')
	workflow.connect(fcon_nodes.alff_pspec,'out_file',fcon_nodes.alff_sqrt,'in_file')

	workflow.connect(fcon_nodes.NVOLS,'nvols',fcon_nodes.calcN1,'nvols')
	workflow.connect(fcon_nodes.TR,'TR',fcon_nodes.calcN1,'TR')
	workflow.connect(infosource,'HP',fcon_nodes.calcN1,'HP')

	workflow.connect(fcon_nodes.NVOLS,'nvols',fcon_nodes.calcN2,'nvols')
	workflow.connect(fcon_nodes.TR,'TR',fcon_nodes.calcN2,'TR')
	workflow.connect(infosource,'LP',fcon_nodes.calcN2,'LP')
	workflow.connect(infosource,'HP',fcon_nodes.calcN2,'HP')

	workflow.connect(fcon_nodes.alff_sqrt,'out_file',fcon_nodes.alff_roi1,'in_file')
	workflow.connect(fcon_nodes.calcN1,'n1',fcon_nodes.alff_roi1,'t_min')
	workflow.connect(fcon_nodes.calcN2,'n2',fcon_nodes.alff_roi1,'t_size')
	workflow.connect(fcon_nodes.alff_roi1,'roi_file',fcon_nodes.alff_sum,'in_file')
	workflow.connect(fcon_nodes.calcN2,('n2',set_op_str),fcon_nodes.alff_sum,'op_string')

	FALFF_NORM_REGISTER_WF()

def extract_seed_name(seed):

	fname = os.path.basename(seed)

	seed_name = (fname.split('.'))[0]

	return seed_name

def printToFile(time_series,seed_ts_dir,seed_name):

	print seed_ts_dir + '/%s.1D'%(seed_name)	
	print time_series
	f = open(seed_ts_dir + '/%s.1D'%(seed_name),'w')

	for tr in time_series:
		print >>f,"%f" %(tr)
	
	f.close()


def RSFC():
	
	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource
	#global datasource_seeds

	workflow.connect(fcon_nodes.nuisance_warp,'out_file',fcon_nodes.RSFC_time_series,'in_file')
	#workflow.connect(datasource_seeds,'seeds',fcon_nodes.RSFC_time_series,'mask')
	workflow.connect(fcon_nodes.nuisance_warp,'out_file',fcon_nodes.RSFC_printToFile,'out_dir')
	workflow.connect(fcon_nodes.RSFC_time_series,'stats',fcon_nodes.RSFC_printToFile,'time_series')
	#workflow.connect(datasource_seeds,'seeds',fcon_nodes.RSFC_printToFile,'seed')
	workflow.connect(fcon_nodes.RSFC_printToFile,'ts_oneD',fcon_nodes.RSFC_corr,'ideal_file')
	workflow.connect(fcon_nodes.nuisance_calc,'out_file',fcon_nodes.RSFC_corr,'in_file')
	workflow.connect(fcon_nodes.RSFC_corr,'out_file',fcon_nodes.RSFC_z_trans,'infile_a')
	workflow.connect(fcon_nodes.RSFC_z_trans,'out_file',fcon_nodes.RSFC_register,'in_file')
	workflow.connect(infosource,'standard',fcon_nodes.RSFC_register,'ref_file')
	workflow.connect(fcon_nodes.reg_fnt, 'fieldcoeff_file',fcon_nodes.RSFC_register,'field_file')
	workflow.connect(fcon_nodes.reg_flirt,'out_matrix_file',fcon_nodes.RSFC_register,'premat')


def VMHC():


	global workflow
	global scripts_dir
	global standard_brain
	global rest_name
	global FSLDIR
	global standard_res
	global infosource

	workflow.connect(fcon_nodes.anat_calc,'out_file',fcon_nodes.VMHC_flirt,'in_file')
	workflow.connect(infosource,'brain_symmetric',fcon_nodes.VMHC_flirt,'reference')
	workflow.connect(fcon_nodes.anat_reorient, 'out_file',fcon_nodes.VMHC_fnt,'in_file')
	#workflow.connect(,fcon_nodes.VMHC_fnt,'affine_file')

def readRestingImageStats(subject,analysisdirectory):

	global rest_name
	global rest_dir

	img = analysisdirectory + '/' + subject + '/' + rest_dir + '/' + rest_name + '.nii.gz'
	cmd = "3dinfo %s | grep 'Time step = *' > %s/TRInfo_N_%s" %(img,scripts_dir,(subject.split('/'))[0])
	print cmd
	sys.stderr.write('\n' + commands.getoutput(cmd))
	
	file = open('%s/TRInfo_N_%s' %(scripts_dir,(subject.split('/'))[0]),'r')

	lines = file.readlines()

	TR = ""
	n_vols = ""
	for line in lines:
		line = line.rstrip('\r\n')
		
		fields = line.split(' ')
		
		i = fields.index('step')
		TR = (fields[i+2]).rstrip('s')

		i = fields.index('steps')
		n_vols = fields[i+2]
		
	file.close()
	
	cmd = "rm -f %s/TRInfo_N_%s" %(scripts_dir,(subject.split('/'))[0])
	print cmd
	sys.stderr.write('\n'+commands.getoutput(cmd))

	return TR+' '+n_vols

def readSeedList():

	global seed_list
	global seeds_dir

	f= open(seed_list,'r')

	lines = f.readlines()

	f.close()

	seeds = []
	for line in lines:
		line = line.rstrip('\r\n')
		seeds.append(os.path.basename(line))
		seeds_dir = os.path.dirname(line)
	
	return seeds

def gatherData(sublist, analysisdirectory):

	global datasource
	global datasource_seeds
	global working_dir
	global anat_name
	global anat_dir
	global rest_name
	global rest_dir
	global seed_list


	datasource = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'anat', 'rest' ]) , name= 'datasource')
	datasource.inputs.base_directory = analysisdirectory
	#datasource.inputs.template = '%s/*/%s.nii.gz'
	datasource.inputs.template = '%s/*/*/%s.nii.gz'
	datasource.inputs.template_args = dict( anat = [['subject_id',anat_name]] , rest = [['subject_id',rest_name]])
	datasource.iterables = ('subject_id', sublist)


	fcon_nodes.setSeedList(seed_list)

	#datasource_seeds = pe.Node( interface = nio.DataGrabber(infields=['seeds'],outfields = ['seeds']), name = 'datasource_seeds')
	#datasource_seeds.inputs.base_directory = seeds_dir
	#datasource_seeds.inputs.template = '%s'
	#datasource_seeds.inputs.seeds = seeds



def getFields():

	fields = ['subject_id','brain_symmetric','symm_standard','standard_res_brain','standard','standard_brain_mask_dil','config_file','PRIOR_CSF','PRIOR_WHITE','nuisance_template','HP','LP','hp','lp']
	return fields

def getInfoSource(sublist, analysisdirectory):
	
	global infosource
	global FSLDIR
	global standard_res
	global prior_dir
	global nuisance_template
	global HP
	global LP
	global hp
	global lp


	formula  = getFields()
	infosource = pe.Node(interface=util.IdentityInterface(fields= formula), name="infosource")

	infosource.inputs.standard_res_brain = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s_brain.nii.gz' %(standard_res))

	infosource.inputs.standard = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
	infosource.inputs.standard_brain_mask_dil = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' %(standard_res))
	infosource.inputs.config_file = os.path.abspath(FSLDIR+'/etc/flirtsch/T1_2_MNI152_%s.cnf'%(standard_res))

	infosource.inputs.PRIOR_CSF = os.path.abspath(prior_dir + '/avg152T1_csf_bin.nii.gz')
	infosource.inputs.PRIOR_WHITE = os.path.abspath(prior_dir + '/avg152T1_white_bin.nii.gz')
	infosource.inputs.HP = float(HP)
	infosource.inputs.LP = float(LP)
	infosource.inputs.hp = float(hp)
	infosource.inputs.lp = float(lp)

	infosource.inputs.brain_symmetric = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s_brain_symmetric.nii.gz' %(standard_res))
	infosource.inputs.symm_standard = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s_symmetric.nii.gz' %(standard_res))

	infosource.inputs.nuisance_template = os.path.abspath(nuisance_template)
	infosource.iterables = ('subject_id', sublist)




def getSink(analysisdirectory):

	global anat_name
	global rest_name
	global datasink
	global datasource
	sink = pe.Node(nio.DataSink(), name = 'sinker-anatomical')
	sink.inputs.base_directory = os.path.abspath(analysisdirectory)
	sink.inputs.regexp_substitutions = [(r"_subject_id_(\w|\d)+", ''),(r"%s.nii.gz/" %(anat_name),''),(r"_anat_(\w|\d)+",'')]


	datasink = pe.MapNode(nio.DataSink(), name = 'sinker-functional' , iterfield = ["container"])
	#datasink = pe.MapNode(nio.DataSink(), name = 'sinker-functional', iterfield = ["container"])
	datasink.inputs.base_directory = os.path.abspath(analysisdirectory)
	#sink.inputs.substitutions = [('_subject_id_','results_'),('_func_','')]
	#sink.inputs.regexp_substitutions = [(r"_subject_id_(\w|\d)+", ''),(r"%s.nii.gz/" %(rest_name),''),(r"_func_(\w|\d)+",'')]
	#datasink.inputs.regexp_substitutions = [(r"_subject_id_(\w|\d)+", ''),(r"%s.nii.gz/" %(rest_name),''),(r":",'/')]
	datasink.inputs.regexp_substitutions = [(r"_subject_id_(\w|\d)+", ''),(r"%s\.nii\.gz/" %(rest_name),''),(r":",'/')]


	return sink

def adjustPath(rest_path):

	import re
	result_path = []
	for r_p in rest_path:
		result_path.append(re.sub(r'\/',':',r_p))

	print str(result_path)
	return result_path

def adjustOnceMore(rest_path):

	import re
	result_path = []

	for r_p in rest_path:
		result_path.append(re.sub(r'(.)+lifeSaver(\d)+\/','',r_p))

	return result_path

def makeOutputConnectionsAnat(datasinkAnat):

	global workflow
	global datasink
	global infosource

	workflow.connect(datasource, 'anat', datasinkAnat, 'container')
	
	## Connect anatpreproc nodes to datasink
	workflow.connect(fcon_nodes.anat_refit,'out_file',datasinkAnat,'@mprage')
	workflow.connect(fcon_nodes.anat_reorient,'out_file',datasinkAnat,'@mprage_RPI')
	workflow.connect(fcon_nodes.anat_skullstrip,'out_file', datasinkAnat,'@mprage_surf')
	workflow.connect(fcon_nodes.anat_calc,'out_file',datasinkAnat,'@mprage_brain')
	

	## Connect funcpreproc nodes to datasink
	workflow.connect(datasource,'rest',fcon_nodes.lifeSaver,'in_file')
	workflow.connect(datasource,('rest',adjustPath),fcon_nodes.lifeSaver,'format_string')
	workflow.connect(fcon_nodes.lifeSaver,('out_file',adjustOnceMore),datasink,'container')

#	workflow.connect(datasource, 'rest', datasink, 'container')
	workflow.connect( fcon_nodes.func_calc , 'out_file', datasink , '@rest_dr' )
	workflow.connect( fcon_nodes.func_refit , 'out_file' , datasink, '@rest_dr_1' )
	workflow.connect( fcon_nodes.func_reorient , 'out_file' , datasink, '@rest_ro' )
	workflow.connect( fcon_nodes.func_tstat , 'out_file' , datasink, '@rest_ro_mean' )
	workflow.connect( fcon_nodes.func_volreg , 'oned_file', datasink, '@rest_mc_1D')
	workflow.connect( fcon_nodes.func_volreg , 'out_file', datasink, '@rest_mc')
	workflow.connect( fcon_nodes.func_automask,'out_file', datasink,'@rest_mask')
	workflow.connect( fcon_nodes.func_calcR,'out_file', datasink,'@rest_ss')
	workflow.connect( fcon_nodes.func_calcI,'out_file', datasink,'@example_func')
	workflow.connect( fcon_nodes.func_despike, 'out_file', datasink, '@rest_ds')
	workflow.connect( fcon_nodes.func_smooth,'out_file' , datasink, '@rest_sm')
	workflow.connect( fcon_nodes.func_scale, 'out_file', datasink, '@rest_gms')
	workflow.connect( fcon_nodes.func_filter, 'out_file', datasink, '@rest_filt')
	workflow.connect( fcon_nodes.func_detrenda, 'out_file' , datasink, '@rest_filt_mean')
	workflow.connect( fcon_nodes.func_detrendb , 'out_file', datasink, '@rest_dt')
	workflow.connect( fcon_nodes.func_detrendc, 'out_file' , datasink, '@rest_pp')
	workflow.connect( fcon_nodes.func_mask, 'out_file' , datasink, '@rest_pp_mask')

	## Connect registration nodes to datasink
	workflow.connect(fcon_nodes.reg_flirt, 'out_file',datasink,'reg.@example_func2highres')
	workflow.connect(fcon_nodes.reg_flirt, 'out_matrix_file',datasink,'reg.@example_func2highresmat')
	workflow.connect(fcon_nodes.reg_xfm1, 'out_file',datasink,'reg.@highres2example_funcmat')
	workflow.connect(fcon_nodes.reg_flirt1, 'out_file',datasinkAnat, 'reg.@highres2standard')
	workflow.connect(fcon_nodes.reg_flirt1, 'out_matrix_file',datasinkAnat, 'reg.@highres2standardmat')
	workflow.connect(fcon_nodes.reg_xfm2, 'out_file',datasinkAnat,'reg.@standard2highresmat')
	workflow.connect(fcon_nodes.reg_xfm3, 'out_file',datasink,'reg.@example_func2standardmat')
	workflow.connect(fcon_nodes.reg_flirt2, 'out_file',datasink,'reg.@example_func2standard')
	workflow.connect(fcon_nodes.reg_xfm4, 'out_file',datasink,'reg.@standard2example_funcmat')
	workflow.connect(fcon_nodes.reg_fnt , 'warped_file',datasinkAnat,'reg.@highres2standard_NL')
	workflow.connect(fcon_nodes.reg_fnt, 'jacobian_file',datasinkAnat,'reg.@highres2standard_jac')
	workflow.connect(fcon_nodes.reg_fnt, 'fieldcoeff_file',datasinkAnat,'reg.@highres2standard_warp')
	workflow.connect(fcon_nodes.reg_warp,'out_file',datasinkAnat,'reg.@example_func2standard_NL')

	## Connect segmentation nodes to datasink	

	workflow.connect( fcon_nodes.seg_segment, 'probability_maps',datasinkAnat,'segment.@seg_segment' )
	workflow.connect( fcon_nodes.seg_flirt , 'out_file', datasink, 'segment.@seg_flirt' )
	workflow.connect( fcon_nodes.seg_smooth, 'out_file', datasink, 'segment.@seg_smooth')
	workflow.connect( fcon_nodes.seg_flirt1, 'out_file', datasink, 'segment.@seg_flirt1')
	workflow.connect( fcon_nodes.seg_smooth1,'out_file', datasink,'segment.@seg_smooth1')
	workflow.connect( fcon_nodes.seg_flirt2, 'out_file', datasink, 'segment.@seg_flirt2')
	workflow.connect( fcon_nodes.seg_thresh,'out_file', datasink,'segment.@seg_thresh')
	workflow.connect( fcon_nodes.seg_mask,'out_file', datasink,'segment.@seg_mask')
	workflow.connect( fcon_nodes.seg_copy,'out_file',  datasink,'segment.@seg_copy')
	workflow.connect( fcon_nodes.seg_flirt3, 'out_file', datasink, 'segment.@seg_flirt3')
	workflow.connect( fcon_nodes.seg_smooth2,'out_file', datasink, 'segment.@seg_smooth2')
	workflow.connect( fcon_nodes.seg_flirt4,'out_file', datasink, 'segment.@seg_flirt4')
	workflow.connect( fcon_nodes.seg_prior1,'out_file', datasink, 'segment.@seg_prior1')
	workflow.connect( fcon_nodes.seg_flirt5,'out_file', datasink,'segment.@seg_flirt5')
	workflow.connect( fcon_nodes.seg_thresh1,'out_file', datasink,'segment.@seg_thresh1')
	workflow.connect( fcon_nodes.seg_mask1,'out_file', datasink,'segment.@seg_mask1')

	## Connect nuisance nodes to datasink	

	workflow.connect(fcon_nodes.nuisance_globalE,'out_file',datasink,'nuisance.@nuisance_globalE')
	workflow.connect(fcon_nodes.nuisance_csf,'out_file',datasink,'nuisance.@nuisance_csf')
	workflow.connect(fcon_nodes.nuisance_wm,'out_file',datasink,'nuisance.@nuisance_wm')
	workflow.connect(fcon_nodes.nuisance_featM,'design_file',datasink,'nuisance.@nuisance_featM')
	workflow.connect(fcon_nodes.nuisance_fgls,'results_dir',datasink,'nuisance.@nuisance_fgls')
	workflow.connect(fcon_nodes.nuisance_fgls,'residual4d',datasink,'nuisance.@nuisance_fgls_')
	workflow.connect(fcon_nodes.nuisance_stat,'out_file',datasink,'nuisance.@nuisance_stat')
	workflow.connect(fcon_nodes.nuisance_calc,'out_file',datasink,'nuisance.@nuisance_calc')
	workflow.connect(fcon_nodes.nuisance_warp,'out_file',datasink,'nuisance.@nuisance_warp')

	## Connect falff nodes to datasink	

	workflow.connect(fcon_nodes.alff_detrend,'out_file',datasink,'ALFF.@alff_detrend')
	workflow.connect(fcon_nodes.alff_smooth,'out_file',datasink,'ALFF.@alff_smooth')
	workflow.connect(fcon_nodes.alff_scale,'out_file',datasink,'ALFF.@alff_scale')
	workflow.connect(fcon_nodes.alff_cp,'out_file',datasink,'ALFF.@alff_cp')
	workflow.connect(fcon_nodes.alff_mean,'out_file',datasink,'ALFF.@alff_mean')
	workflow.connect(fcon_nodes.alff_pspec,'out_file',datasink,'ALFF.@alff_pspec')
	workflow.connect(fcon_nodes.alff_sqrt,'out_file',datasink,'ALFF.@alff_sqrt')
	workflow.connect(fcon_nodes.alff_roi1,'roi_file',datasink,'ALFF.@alff_roi1')
	workflow.connect(fcon_nodes.alff_sum,'out_file',datasink,'ALFF.@alff_sum')
	workflow.connect(fcon_nodes.alff_falff,'out_file',datasink,'ALFF.@alff_falff')
	workflow.connect(fcon_nodes.alff_falff1,'out_file',datasink,'ALFF.@alff_falff1')

	## connect RSFC nodes to datasink
	workflow.connect(fcon_nodes.RSFC_printToFile,'ts_oneD',datasink,'RSFC.@ts_oneD')
	workflow.connect(fcon_nodes.RSFC_corr,'out_file',datasink,'RSFC.@RSFC_corr')
	workflow.connect(fcon_nodes.RSFC_z_trans,'out_file',datasink,'RSFC.@RSFC_z_trans')
	workflow.connect(fcon_nodes.RSFC_register,'out_file',datasink,'RSFC.@RSFC_register')





	

def processS(sublist, analysisdirectory):


	global infosource
	global datasource
	global workflow
	global datasink

	numCores = int(sys.argv[1])
	
	getInfoSource(sublist, analysisdirectory)
	gatherData(sublist, analysisdirectory)
	
	wfname =  'fcon1000'
	workflow = pe.Workflow(name=wfname)
	workflow.base_dir = working_dir
	workflow.config = {'execution' : {'stop_on_first_crash' : True}}
	anatpreproc()
	funcpreproc()
	registration()
	segment()
	nuisance()
	falff()
	RSFC()
	datasinkAnat = getSink(analysisdirectory)
	makeOutputConnectionsAnat(datasinkAnat)
	#VMHC()
	workflow.run(plugin='MultiProc', plugin_args={'n_procs' : numCores})
		
#		f.write('finished : %s\n'%(line))
#		f.close()		
#	except:
#		f = open(logFile,'a')
#		f.write('failed : %s\n'%(line))
#		f.close()		

def checkProcessedSubject(dir,subject):


	checkList = ['mprage_brain.nii.gz']

	for file in checkList:
		if ('mprage' in file):
			if os.path.isdir(dir + '/' +subject + '/session_1/anat_1/'):

				fpath_1 = dir + '/' +subject + '/session_1/anat_1/' + file
				if os.path.isfile(fpath_1):
					return 1
				else:
					files = os.listdir(dir + '/' +subject + '/session_1/anat_1/')
					if(len(files) > 1):
						return 1
					else:
						return 0

			elif os.path.isdir(dir + '/' +subject + '/session_2/anat_1/'):	
				fpath_1 = dir + '/' +subject + '/session_2/anat_1/' + file
				if os.path.isfile(fpath_1):
					return 1
				else:
					files = os.listdir(dir + '/' +subject + '/session_2/anat_1/')
					if (len(files) > 1):
						return 1
					else:
						return 0
	return 1

def processSubjects(analysisdirectory,subject_list):

	#global subject_list

	
	try:
		fp = open(subject_list,'r')

		lines = fp.readlines()

		fp.close()
		line1 = []
		cnt = 0
		for line in lines:
			line = line.rstrip('\r\n')

			flag = 0
			print 'inside processSubjects ; before check processed subject ',line
			#flag = checkProcessedSubject(analysisdirectory, line)

			print 'returned %d' %(flag)
			if (flag == 0):
				if ( os.path.isdir(analysisdirectory + '/' + line + '/session_1/rest_1') or os.path.isdir(analysisdirectory + '/' + line + '/session_2/rest_1')):
					line1.append(line)
				
					cnt += 1


		if not (cnt == 0):
			#process =  pool.map(processS,line1)
			processS(line1,analysisdirectory)
			#pool.close()
			#pool.join()

	except:
		raise
def readSubjects():


	fp = open(batch_list,'r')

	lines = fp.readlines()
	line1 = []

	for line in lines:
		line = line.rstrip('\r\n')
		line1.append(line)
	
	fp.close()

	processes = [ Process(target=processSubjects, args=((l.split(' '))[0],(l.split(' '))[1])) for l in line1 ]
	processes[0].start()	
#	l = line1[0]
#	print l
#	processSubjects((l.split(' '))[0],(l.split(' '))[1])

def readDirSetup():

	global FSLDIR
	global scripts_dir
	global working_dir
	global seed_list
	global batch_list
	global rest_name
	global rest_dir
	global anat_name
	global anat_dir
	global standard_res
	global standard_brain
	global logFile
	global prior_dir
	global nuisance_template
	global recon_subjects
	global HP
	global LP
	global hp
	global lp
	parsermap = {}
	
	parser = SafeConfigParser()
	parser.read('dir_setup.ini')

	for section in parser.sections():

		for variable , value in parser.items(section):
			print variable + ' ' + value
			parsermap[variable] = value
	
	FSLDIR = parsermap['fsldir']
	scripts_dir = parsermap['scripts_dir']
	prior = parsermap['prior_dir']
	working_dir = parsermap['working_dir']
	seed_list = scripts_dir + '/' + parsermap['seed_file']
	batch_list = scripts_dir + '/' + parsermap['batch_file']
	rest_name = parsermap['rest_name']
	rest_dir = parsermap['rest_dir']
	anat_name = parsermap['anat_name']
	anat_dir = parsermap['anat_dir']
	standard_res = parsermap['standard_res']
	nuisance_template = parsermap['nuisance_template']
	prior_dir = prior + '/%s' %(standard_res)
	standard_brain = FSLDIR + '/data/standard/MNI152_T1_%s_brain.nii.gz' %(standard_res)
	logFile = parsermap['logfile']
	recon_subjects = parsermap['recon_subjects']
	HP = parsermap['_hp_']
	LP = parsermap['_lp_']
	hp = parsermap['hp']
	lp = parsermap['lp']


def main():

	if ( len(sys.argv) < 2 ):
		sys.stderr.write("./fcon_pipeline.py <Number_of_cores>")
	else:
	
		readDirSetup()
		readSubjects()

	

if __name__ == "__main__":

	sys.exit(main())
