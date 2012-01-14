#!/usr/local/epd-7.0-2-rh5-x86_64/bin/python
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


def extract_subjectID(out_file):

	import sys
	outs = (out_file.split('_subject_id_'))[1]

	out_file = (outs.split('/'))[0]

	sys.stderr.write('\n>>>>>D<<<<< ' + out_file)

	return out_file

def anatpreproc():

	global anat_name
	global infosource
	global datasource
	global workflow
	
	workflow.connect(datasource, 'anat',jhu_nodes.anat_refit,'in_file')
	workflow.connect(jhu_nodes.anat_refit,'out_file',jhu_nodes.anat_reorient,'in_file')
	workflow.connect(jhu_nodes.anat_reorient,'out_file',jhu_nodes.anat_skullstrip,'in_file')
	workflow.connect(jhu_nodes.anat_skullstrip,'out_file',jhu_nodes.anat_calc,'infile_b')
	workflow.connect(jhu_nodes.anat_reorient,'out_file',jhu_nodes.anat_calc,'infile_a')
	workflow.connect(jhu_nodes.anat_calc,'out_file', jhu_nodes.anat_reconall,'T1_files')
	workflow.connect(jhu_nodes.anat_reorient,('out_file',extract_subjectID),jhu_nodes.anat_reconall,'subject_id')
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




def gatherData(sublist, analysisdirectory):

	global datasource
	global working_dir
	global anat_name
	global anat_dir
	global rest_name
	global rest_dir


	datasource = pe.Node( interface = nio.DataGrabber(infields=['subject_id'], outfields = [ 'anat', 'rest' ]) , name= 'datasource')
	datasource.inputs.base_directory = analysisdirectory
	datasource.inputs.template = '%s/*/%s.nii.gz'
	datasource.inputs.template_args = dict( anat = [['subject_id',anat_name]] , rest = [['subject_id',rest_name]])
	datasource.iterables = ('subject_id', sublist)


def getFields():

	fields = ['subject_id','recon_subjects']
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
	infosource.inputs.recon_subjects = os.path.abspath(recon_subjects)
	infosource.inputs.subject_id = sublist




def getSink(analysisdirectory):

	sink = pe.Node(nio.DataSink(), name = 'sinker')
	sink.inputs.base_directory = os.path.abspath(analysisdirectory)
	sink.inputs.substitutions = [('_subject_id_','results_'),('_func_','')]
	return sink


def makeOutputConnections():

	global workflow
	global datasink
	global infosource

	workflow.connect(jhu_nodes.anat_reorient,('out_file',extract_subjectID), datasink, 'container')
	
	## Connect anatpreproc nodes to datasink
	workflow.connect(jhu_nodes.anat_refit,'out_file',datasink,'anat-results.@refit')
	workflow.connect(jhu_nodes.anat_reorient,'out_file',datasink,'anat-results.@reorient')
	workflow.connect(jhu_nodes.anat_skullstrip,'out_file', datasink,'anat-results.@skullstrip')
	workflow.connect(jhu_nodes.anat_calc,'out_file',datasink,'anat-results.@brain')

	## Connect funcpreproc nodes to datasink

	workflow.connect( jhu_nodes.func_calc , 'out_file', datasink , 'rest-results.@rest_dr' )
	workflow.connect( jhu_nodes.func_refit , 'out_file' , datasink, 'rest-results.@rest_dr_1' )
	workflow.connect( jhu_nodes.func_reorient , 'out_file' , datasink, 'rest-results.@rest_ro' )
	workflow.connect( jhu_nodes.func_tstat , 'out_file' , datasink, 'rest-results.@rest_ro_mean' )
	workflow.connect( jhu_nodes.func_volreg , 'oned_file', datasink, 'rest-results.@rest_mc_1D')
	workflow.connect( jhu_nodes.func_volreg , 'out_file', datasink, 'rest-results.@rest_mc')
	workflow.connect( jhu_nodes.func_bbreg, 'out_reg_file', datasink, 'rest-results.@out_reg_file')
	workflow.connect( jhu_nodes.func_sampler_lh, 'out_file', datasink, 'rest-results.@lh_out_file')
	workflow.connect( jhu_nodes.func_sampler_lh, 'out_file', datasink, 'rest-results.@rh_out_file')


def processS(sublist, analysisdirectory):


	global infosource
	global datasource
	global workflow
	global datasink

	numCores = int(sys.argv[1])
	
	getInfoSource(sublist, analysisdirectory)
	gatherData(sublist, analysisdirectory)
	
	wfname =  'SurfaceRegistration'
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
	#uncomment once eric sets up FSLDIR
	FSLDIR = parsermap['fsldir']
	scripts_dir = parsermap['scripts_dir']
	working_dir = parsermap['working_dir']
	batch_list = scripts_dir + '/' + parsermap['batch_file']
	rest_name = parsermap['rest_name']
	anat_name = parsermap['anat_name']
	#uncomment once eric sets up FSLDIR
	recon_subjects = parsermap['recon_subjects']

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
