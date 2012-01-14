#!/usr/local/epd-7.0-2-rh5-x86_64/bin/python
import e_afni
import sys
import os
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util
from nipype.interfaces.utility import Rename
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
from jhu_pipeline import working_dir

FWHM = 6.0
sigma = 2.5479870901
hp = 0.005
lp = 0.1
HP = 0.01
LP = 0.1

def getDIR():

	import os
	from jhu_pipeline import working_dir

	print '>>>>>>>> ' + working_dir + ' <<<<<< ' +'fsub_recon'
	return os.path.abspath(working_dir + '/fsub_recon')


def getImgTR(in_files):

	out = []
	from nibabel import load
	if(isinstance(in_files,list)):
		for in_file in in_files:
			img = load(in_file)
			hdr = img.get_header()
			tr = int(hdr.get_zooms()[3])
			if tr > 10:
				tr = float(float(tr)/1000.0)
			out.append(tr)
		return out
	else:
		img = load(in_files)
		hdr = img.get_header()
		tr = int(hdr.get_zooms()[3])
		if tr > 10:
			tr = float(float(tr)/1000.0)
		return [tr]

def getImgNVols(in_files):

	out = []
	from nibabel import load
	if(isinstance(in_files,list)):
		for in_file in in_files:
			img = load(in_file)
			hdr = img.get_header()
			out.append(int(hdr.get_data_shape()[3]))
		return out
	
	else:
		img = load(in_files)
		hdr = img.get_header()
		n_vol = int(hdr.get_data_shape()[3])
		return [n_vol]





#TR and n_vols nodes

TR = pe.Node(util.Function(input_names = ['in_files'], output_names = ['TR'], function = getImgTR), name='TR')
NVOLS = pe.Node(util.Function(input_names = ['in_files'], output_names = ['nvols'], function = getImgNVols), name='NVOLS')

DIR = pe.Node(util.Function(input_names = [], output_names = ['DIR'], function = getDIR), name='DIR')

#anatomical nodes


anat_refit = pe.Node(interface=e_afni.Threedrefit(), name='anat_refit')
anat_refit.inputs.deoblique = True

anat_reorient = pe.Node(interface=e_afni.Threedresample(), name='anat_reorient')
anat_reorient.inputs.orientation = 'RPI'

anat_skullstrip = pe.Node(interface= e_afni.ThreedSkullStrip(), name='anat_skullstrip')
anat_skullstrip.inputs.options = '-no_use_edge'


anat_calc = pe.Node(interface=e_afni.Threedcalc(), name='anat_calc')
anat_calc.inputs.expr = '\'a*step(b)\''
anat_calc.inputs.out_file = 'mprage_brain.nii.gz'


anat_reconall = pe.Node(interface=fs.ReconAll(), name="anat_reconall")
anat_reconall.inputs.directive = 'all'
#anat_reconall.inputs.flags = '-use-gpu'

#funcpreproc Nodes

func_calc = pe.MapNode(interface = e_afni.Threedcalc(), name='func_calc', iterfield = ["infile_a","stop_idx"])
func_calc.inputs.start_idx = 0
func_calc.inputs.expr = '\'a\''
func_calc.inputs.out_file = 'rest_dr.nii.gz'

func_refit = pe.MapNode(interface=e_afni.Threedrefit(), name='func_refit', iterfield = ["in_file"])
func_refit.inputs.deoblique = True

func_slice_time_correction = pe.MapNode(interface=fsl.SliceTimer(), name = 'func_slice_time_correction', iterfield=["in_file","time_repetition"])

func_slice_time_correction.inputs.interleaved = True 

func_reorient = pe.MapNode(interface=e_afni.Threedresample(), name='func_reorient', iterfield = ["in_file"])
func_reorient.inputs.orientation = 'RPI'

func_tstat = pe.MapNode(interface=e_afni.ThreedTstat(), name='func_tstat', iterfield = ["in_file"])
func_tstat.inputs.args = "-mean"
func_tstat.inputs.out_file = 'rest_ro_mean.nii.gz'

func_volreg = pe.MapNode(interface=e_afni.Threedvolreg(), name='func_volreg', iterfield = ["in_file","basefile"])
func_volreg.inputs.other = '-Fourier -twopass'
func_volreg.inputs.zpad = '4'
func_volreg.inputs.oned_file = 'rest_mc.1D'
func_volreg.inputs.out_file = 'rest_mc.nii.gz'

func_bbreg = pe.MapNode(interface = fs.BBRegister(init='fsl', contrast_type='t2',registered_file = True, out_fsl_file = True), name = 'func_bbreg', iterfield = ["source_file"] ) 

func_sampler_lh = pe.MapNode(interface = fs.SampleToSurface(hemi = "lh"),name ='func_sampler_lh', iterfield = ["source_file","reg_file"])
func_sampler_lh.inputs.no_reshape = True
func_sampler_lh.inputs.interp_method = 'trilinear'
func_sampler_lh.inputs.sampling_method = "point"
func_sampler_lh.inputs.sampling_range = 0.5
func_sampler_lh.inputs.sampling_units = "frac"

func_sampler_rh = pe.MapNode(interface = fs.SampleToSurface(hemi = "rh"),name ='func_sampler_rh', iterfield = ["source_file","reg_file"])
func_sampler_rh.inputs.no_reshape = True
func_sampler_rh.inputs.interp_method = 'trilinear'
func_sampler_rh.inputs.sampling_method = "point"
func_sampler_rh.inputs.sampling_range = 0.5
func_sampler_rh.inputs.sampling_units = "frac"
