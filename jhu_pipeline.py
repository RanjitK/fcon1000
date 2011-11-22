#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import e_afni
import jhu_nodes
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
rest_name = ''
working_dir = ''
standard_brain = ''
batch_list = ''
standard_res = ''
logFile = ''
nuisance_template = ''
recon_subjects = ''
infosource = None
datasource = None
datasink = None
workflow = None
#*******************************************************************************************



def anatpreproc():

	global anat_name
	global infosource
	global datasource
	global workflow
	
	#workflow.connect( infosource,'subject_id',datasource, 'subject_id' )
	workflow.connect(datasource, 'anat',jhu_nodes.anat_refit,'in_file')
	workflow.connect(jhu_nodes.anat_refit,'out_file',jhu_nodes.anat_reorient,'in_file')
	workflow.connect(jhu_nodes.anat_reorient,'out_file',jhu_nodes.anat_skullstrip,'in_file')
	workflow.connect(jhu_nodes.anat_skullstrip,'out_file',jhu_nodes.anat_calc,'infile_b')
	workflow.connect(jhu_nodes.anat_reorient,'out_file',jhu_nodes.anat_calc,'infile_a')
#	workflow.connect(jhu_nodes.anat_calc,'out_file', jhu_nodes.anat_reconall,'T1_files')
	workflow.connect(infosource,'subject_id',jhu_nodes.anat_reconall,'subject_id')
	workflow.connect(infosource,'recon_subjects',jhu_nodes.anat_reconall,'subjects_dir')
	
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


	workflow.connect(datasource, 'rest',jhu_nodes.TR,'in_files')
	workflow.connect(datasource, 'rest',jhu_nodes.NVOLS,'in_files')
	workflow.connect(jhu_nodes.NVOLS,('nvols',last_vol),jhu_nodes.func_calc,'stop_idx')
	workflow.connect(datasource, 'rest', jhu_nodes.func_calc,'infile_a')
	workflow.connect( jhu_nodes.func_calc , 'out_file', jhu_nodes.func_refit , 'in_file' )
	workflow.connect( jhu_nodes.func_refit , 'out_file' , jhu_nodes.func_slice_time_correction , 'in_file' )
	workflow.connect( jhu_nodes.TR, 'TR', jhu_nodes.func_slice_time_correction,'time_repetition') 
	#workflow.connect( jhu_nodes.func_refit , 'out_file' , jhu_nodes.func_reorient , 'in_file' )
	workflow.connect(jhu_nodes.func_slice_time_correction,'slice_time_corrected_file',jhu_nodes.func_reorient , 'in_file')
	workflow.connect( jhu_nodes.func_reorient , 'out_file' , jhu_nodes.func_tstat , 'in_file' )
	workflow.connect( jhu_nodes.func_reorient , 'out_file' , jhu_nodes.func_volreg , 'in_file' )
	workflow.connect( jhu_nodes.func_tstat , 'out_file' , jhu_nodes.func_volreg , 'basefile' )
	workflow.connect( jhu_nodes.func_volreg , 'out_file', jhu_nodes.func_bbreg, 'source_file')

	workflow.connect( jhu_nodes.anat_reconall,'subjects_dir',jhu_nodes.func_bbreg,'subjects_dir')
	workflow.connect( jhu_nodes.anat_reconall,'subject_id',jhu_nodes.func_bbreg,'subject_id')

	workflow.connect(jhu_nodes.func_bbreg,'out_reg_file',jhu_nodes.func_sampler_lh,'reg_file')
	workflow.connect(jhu_nodes.func_volreg,'out_file', jhu_nodes.func_sampler_lh,'source_file')


	workflow.connect(jhu_nodes.func_bbreg,'out_reg_file',jhu_nodes.func_sampler_rh,'reg_file')
	workflow.connect(jhu_nodes.func_volreg,'out_file', jhu_nodes.func_sampler_rh,'source_file')

#--------
#	workflow.connect( jhu_nodes.func_volreg , 'out_file', jhu_nodes.func_automask , 'in_file')
#	workflow.connect( jhu_nodes.func_volreg,'out_file', jhu_nodes.func_calcR,'infile_a')
#	workflow.connect( jhu_nodes.func_automask,'out_file', jhu_nodes.func_calcR,'infile_b')
#	workflow.connect( jhu_nodes.func_calcR,'out_file',jhu_nodes.func_calcI,'infile_a')
#	workflow.connect( jhu_nodes.func_calcR, 'out_file', jhu_nodes.func_despike,'in_file')
#	workflow.connect( jhu_nodes.func_despike, 'out_file', jhu_nodes.func_smooth , 'in_file')
#	workflow.connect( jhu_nodes.func_automask,'out_file',jhu_nodes.func_smooth,'operand_files')
#	workflow.connect( jhu_nodes.func_smooth,'out_file' , jhu_nodes.func_scale , 'in_file')
#	workflow.connect( jhu_nodes.func_scale, 'out_file', jhu_nodes.func_filter, 'in_file')
#	workflow.connect( jhu_nodes.func_filter, 'out_file', jhu_nodes.func_detrenda , 'in_file')
#	workflow.connect( jhu_nodes.func_filter, 'out_file', jhu_nodes.func_detrendb, 'in_file')
#	workflow.connect( jhu_nodes.func_detrenda, 'out_file' , jhu_nodes.func_detrendc , 'infile_a')
#	workflow.connect( jhu_nodes.func_detrendb , 'out_file', jhu_nodes.func_detrendc , 'infile_b')
#	workflow.connect( jhu_nodes.func_detrendc, 'out_file' , jhu_nodes.func_mask, 'in_file')



def gatherData(sublist, analysisdirectory):

	global datasource
	global working_dir
	global anat_name
	global anat_dir
	global rest_name
	global rest_dir


	datasource = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'anat', 'rest' ]) , name= 'datasource')
	datasource.inputs.base_directory = analysisdirectory
	datasource.inputs.template = '%s/*/*/%s.nii.gz'
	datasource.inputs.template_args = dict( anat = [['subject_id',anat_name]] , rest = [['subject_id',rest_name]])
	datasource.iterables = ('subject_id', sublist)


def getFields():

	fields = ['subject_id','standard_res_brain','standard','standard_brain_mask_dil','config_file','PRIOR_CSF','PRIOR_WHITE','nuisance_template','recon_subjects']
	return fields

def getInfoSource(sublist, analysisdirectory):
	
	global infosource
	global FSLDIR
	global standard_res
	global prior_dir
	global nuisance_template
	global recon_subjects

	formula  = getFields()
	infosource = pe.Node(interface=util.IdentityInterface(fields= formula), name="infosource")

	infosource.inputs.standard_res_brain = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s_brain.nii.gz' %(standard_res))

	infosource.inputs.standard = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
	infosource.inputs.standard_brain_mask_dil = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' %(standard_res))
	infosource.inputs.config_file = os.path.abspath(FSLDIR+'/etc/flirtsch/T1_2_MNI152_%s.cnf'%(standard_res))

	infosource.inputs.PRIOR_CSF = os.path.abspath(prior_dir + '/avg152T1_csf_bin.nii.gz')
	infosource.inputs.PRIOR_WHITE = os.path.abspath(prior_dir + '/avg152T1_white_bin.nii.gz')
	infosource.inputs.nuisance_template = os.path.abspath(nuisance_template)
	infosource.inputs.recon_subjects = os.path.abspath(recon_subjects)
	infosource.iterables = ('subject_id', sublist)




def getSink(analysisdirectory):

	sink = pe.Node(nio.DataSink(), name = 'sinker')
	sink.inputs.base_directory = os.path.abspath(analysisdirectory)
	sink.inputs.substitutions = [('_subject_id_','results_'),('_func_','')]
	return sink


def makeOutputConnections():

	global workflow
	global datasink
	global infosource

	workflow.connect(infosource, 'subject_id', datasink, 'container')
	
	## Connect anatpreproc nodes to datasink
	workflow.connect(jhu_nodes.anat_refit,'out_file',datasink,'anat.@refit')
	workflow.connect(jhu_nodes.anat_reorient,'out_file',datasink,'anat.@reorient')
	workflow.connect(jhu_nodes.anat_skullstrip,'out_file', datasink,'anat.@skullstrip')
	workflow.connect(jhu_nodes.anat_calc,'out_file',datasink,'anat.@brain')
#	workflow.connect(jhu_nodes.anat_reconall,'annot',datasink,'ReconAll.@annot')
#	workflow.connect(jhu_nodes.anat_reconall,'aparc_aseg',datasink,'ReconAll.@aparc_aseg')
#	workflow.connect(jhu_nodes.anat_reconall,'aseg',datasink,'ReconAll.@aseg')
#	workflow.connect(jhu_nodes.anat_reconall,'brain',datasink,'ReconAll.@brain')
#	workflow.connect(jhu_nodes.anat_reconall,'brainmask',datasink,'ReconAll.@brainmask')
#	workflow.connect(jhu_nodes.anat_reconall,'curv',datasink,'ReconAll.@curv')
#	workflow.connect(jhu_nodes.anat_reconall,'filled',datasink,'ReconAll.@filled')
#	workflow.connect(jhu_nodes.anat_reconall,'inflated',datasink,'ReconAll.@inflated')
#	workflow.connect(jhu_nodes.anat_reconall,'label',datasink,'ReconAll.@label')
#	workflow.connect(jhu_nodes.anat_reconall,'norm',datasink,'ReconAll.@norm')
#	workflow.connect(jhu_nodes.anat_reconall,'nu',datasink,'ReconAll.@nu')
#	workflow.connect(jhu_nodes.anat_reconall,'orig',datasink,'ReconAll.@orig')
#	workflow.connect(jhu_nodes.anat_reconall,'pial',datasink,'ReconAll.@pial')
#	workflow.connect(jhu_nodes.anat_reconall,'rawavg',datasink,'ReconAll.@rawavg')
#	workflow.connect(jhu_nodes.anat_reconall,'ribbon',datasink,'ReconAll.@ribbon')
#	workflow.connect(jhu_nodes.anat_reconall,'smoothwm',datasink,'ReconAll.@smoothwm')
#	workflow.connect(jhu_nodes.anat_reconall,'sphere',datasink,'ReconAll.@sphere')
#	workflow.connect(jhu_nodes.anat_reconall,'sphere_reg',datasink,'ReconAll.@sphere_reg')
#	workflow.connect(jhu_nodes.anat_reconall,'sulc',datasink,'ReconAll.@sulc')
#	workflow.connect(jhu_nodes.anat_reconall,'thickness',datasink,'ReconAll.@thickness')
#	workflow.connect(jhu_nodes.anat_reconall,'volume',datasink,'ReconAll.@volume')
#	workflow.connect(jhu_nodes.anat_reconall,'white',datasink,'ReconAll.@white')
#	workflow.connect(jhu_nodes.anat_reconall,'wm',datasink,'ReconAll.@wm')
#	workflow.connect(jhu_nodes.anat_reconall,'wmparc',datasink,'ReconAll.@wmparc')
	

	## Connect funcpreproc nodes to datasink

	workflow.connect( jhu_nodes.func_calc , 'out_file', datasink , 'rest.@rest_dr' )
	workflow.connect( jhu_nodes.func_refit , 'out_file' , datasink, 'rest.@rest_dr_1' )
	workflow.connect( jhu_nodes.func_reorient , 'out_file' , datasink, 'rest.@rest_ro' )
	workflow.connect( jhu_nodes.func_tstat , 'out_file' , datasink, 'rest.@rest_ro_mean' )
	workflow.connect( jhu_nodes.func_volreg , 'oned_file', datasink, 'rest.@rest_mc_1D')
	workflow.connect( jhu_nodes.func_volreg , 'out_file', datasink, 'rest.@rest_mc')

#	workflow.connect( jhu_nodes.func_bbreg, 'min_cost_file',datasink,'bbreg.@min_cost_file')
#	workflow.connect( jhu_nodes.func_bbreg, 'out_fsl_file', datasink, 'bbreg.@out_fsl_file')
#	workflow.connect( jhu_nodes.func_bbreg, 'out_reg_file', datasink, 'bbreg.@out_reg_file')
#	workflow.connect( jhu_nodes.func_bbreg, 'registered_file', datasink,'bbreg.@registered_file')
#
#	workflow.connect( jhu_nodes.func_sampler_lh, 'hits_file', datasink, 'fs_rg_surface_lh.@hits_file')
	workflow.connect( jhu_nodes.func_sampler_lh, 'out_file', datasink, 'fs_rg_surface_lh.@out_file')
	workflow.connect( jhu_nodes.func_sampler_lh, 'out_file', datasink, 'fs_rg_surface_rh.@out_file')
#	workflow.connect( jhu_nodes.func_sampler_lh, 'vox_file', datasink, 'fs_rg_surface_lh.@vox_file')
#----
#	workflow.connect( jhu_nodes.func_automask,'out_file', datasink,'@rest_mask')
#	workflow.connect( jhu_nodes.func_calcR,'out_file', datasink,'@rest_ss')
#	workflow.connect( jhu_nodes.func_calcI,'out_file', datasink,'@example_func')
#	workflow.connect( jhu_nodes.func_despike, 'out_file', datasink, '@rest_ds')
#	workflow.connect( jhu_nodes.func_smooth,'out_file' , datasink, '@rest_sm')
#	workflow.connect( jhu_nodes.func_scale, 'out_file', datasink, '@rest_gms')
#	workflow.connect( jhu_nodes.func_filter, 'out_file', datasink, '@rest_filt')
##	workflow.connect( jhu_nodes.func_detrenda, 'out_file' , datasink, '@rest_filt_mean')
#	workflow.connect( jhu_nodes.func_detrendb , 'out_file', datasink, '@rest_dt')
#	workflow.connect( jhu_nodes.func_detrendc, 'out_file' , datasink, '@rest_pp')
#	workflow.connect( jhu_nodes.func_mask, 'out_file' , datasink, '@rest_pp_mask')

	## Connect registration nodes to datasink
	
#	workflow.connect(jhu_nodes.reg_flirt, 'out_file',datasink,'@example_func2highres')
#	workflow.connect(jhu_nodes.reg_flirt, 'out_matrix_file',datasink,'@example_func2highresmat')
#	workflow.connect(jhu_nodes.reg_xfm1, 'out_file',datasink,'@highres2example_funcmat')
#	workflow.connect(jhu_nodes.reg_flirt1, 'out_file',datasink, '@highres2standard')
#	workflow.connect(jhu_nodes.reg_flirt1, 'out_matrix_file',datasink, '@highres2standardmat')
#	workflow.connect(jhu_nodes.reg_xfm2, 'out_file',datasink,'@standard2highresmat')
#	workflow.connect(jhu_nodes.reg_xfm3, 'out_file',datasink,'@example_func2standardmat')
#	workflow.connect(jhu_nodes.reg_flirt2, 'out_file',datasink,'@example_func2standard')
#	workflow.connect(jhu_nodes.reg_xfm4, 'out_file',datasink,'@standard2example_funcmat')
#	workflow.connect(jhu_nodes.reg_fnt , 'warped_file',datasink,'@highres2standard_NL')
#	workflow.connect(jhu_nodes.reg_fnt, 'jacobian_file',datasink,'@highres2standard_jac')
#	workflow.connect(jhu_nodes.reg_fnt, 'fieldcoeff_file',datasink,'@highres2standard_warp')
#	workflow.connect(jhu_nodes.reg_warp,'out_file',datasink,'@example_func2standard_NL')

	## Connect segmentation nodes to datasink	

#	workflow.connect( jhu_nodes.seg_segment, 'probability_maps',datasink,'@seg_segment' )
#	workflow.connect( jhu_nodes.seg_flirt , 'out_file', datasink, '@seg_flirt' )
#	workflow.connect( jhu_nodes.seg_smooth, 'out_file', datasink, '@seg_smooth')
#	workflow.connect( jhu_nodes.seg_flirt1, 'out_file', datasink, '@seg_flirt1')
#	workflow.connect( jhu_nodes.seg_smooth1,'out_file', datasink,'@seg_smooth1')
#	workflow.connect( jhu_nodes.seg_flirt2, 'out_file', datasink, '@seg_flirt2')
#	workflow.connect( jhu_nodes.seg_thresh,'out_file', datasink,'@seg_thresh')
#	workflow.connect( jhu_nodes.seg_mask,'out_file', datasink,'@seg_mask')
#	workflow.connect( jhu_nodes.seg_copy,'out_file',  datasink,'@seg_copy')
#	workflow.connect( jhu_nodes.seg_flirt3, 'out_file', datasink, '@seg_flirt3')
#	workflow.connect( jhu_nodes.seg_smooth2,'out_file', datasink, '@seg_smooth2')
#	workflow.connect( jhu_nodes.seg_flirt4,'out_file', datasink, '@seg_flirt4')
#	workflow.connect( jhu_nodes.seg_prior1,'out_file', datasink, '@seg_prior1')
#	workflow.connect( jhu_nodes.seg_flirt5,'out_file', datasink,'@seg_flirt5')
#	workflow.connect( jhu_nodes.seg_thresh1,'out_file', datasink,'@seg_thresh1')
#	workflow.connect( jhu_nodes.seg_mask1,'out_file', datasink,'@seg_mask1')


def processS(sublist, analysisdirectory):


	global infosource
	global datasource
	global workflow
	global datasink

	numCores = int(sys.argv[1])
	
	getInfoSource(sublist, analysisdirectory)
	gatherData(sublist, analysisdirectory)
	
	wfname =  'jhuTest'
	workflow = pe.Workflow(name=wfname)
	workflow.base_dir = working_dir
	anatpreproc()
	funcpreproc()
	datasink = getSink(analysisdirectory)
	makeOutputConnections()
	workflow.run(plugin='MultiProc', plugin_args={'n_procs' : numCores})
		

def processSubjects(analysisdirectory,subject_list):

	try:
		fp = open(subject_list,'r')

		lines = fp.readlines()

		fp.close()
		line1 = []
		cnt = 0
		for line in lines:
			line = line.rstrip('\r\n')

			line1.append(line)

		processS(line1,analysisdirectory)

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
	batch_list = scripts_dir + '/' + parsermap['batch_file']
	rest_name = parsermap['rest_name']
	anat_name = parsermap['anat_name']
	standard_res = parsermap['standard_res']
	nuisance_template = parsermap['nuisance_template']
	prior_dir = prior + '/%s' %(standard_res)
	standard_brain = FSLDIR + '/data/standard/MNI152_T1_%s_brain.nii.gz' %(standard_res)
	recon_subjects = parsermap['recon_subjects']
	logFile = parsermap['logfile']

	cmd = 'export SUBJECTS_DIR="%s"'%(recon_subjects)
	print cmd
	sys.stderr.write('>>><<<\n'+commands.getoutput(cmd))
	#os.system('export SUBJECTS_DIR="%s"'%(recon_subjects))


def main():

	if ( len(sys.argv) < 2 ):
		sys.stderr.write("./fcon_pipeline.py <Number_of_cores>")
	else:
	
		readDirSetup()
		readSubjects()

	

if __name__ == "__main__":

	sys.exit(main())
