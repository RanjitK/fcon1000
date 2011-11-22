#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
import e_afni
import sys
import os
import commands
import nipype.pipeline.engine as pe
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
import nipype.interfaces.fsl as fsl
import nipype.interfaces.io as nio
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

FWHM = 6.0
sigma = 2.5479870901
#n_lp = float(HP) * float(int(n_vols)) * float (TR)
#n1 = float("%1.0f" %(float(n_lp - 1.0)))
#n_hp = float(LP) * float(int(n_vols)) * float (TR)
#n2 = float("%1.0f" %(float(n_hp - n_lp + 1.0)))
seed_list = []

def setSeedList(seed_file):

	global seed_list

	f= open(seed_file,'r')
	lines = f.readlines()
	f.close()

	for line in lines:
		line = line.rstrip('\r\n')
		seed_list.append(line)
	
	print str(seed_list)


def getN1(TR,nvols,HP):

	n_lp = float(HP) * float(int(nvols)) * float (TR)
	n1 = int("%1.0f" %(float(n_lp - 1.0)))

	return n1

def getN2(TR,nvols,LP,HP):

	n_lp = float(HP) * float(int(nvols)) * float (TR)
	n_hp = float(LP) * float(int(nvols)) * float (TR)
	n2 = int("%1.0f" %(float(n_hp - n_lp + 1.0)))

	return n2


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




def createMC(in_file,pp):

	import os
	import sys
	import commands
	EV_Lists = []
	idx = 0
	for file in in_file:
		print file
		parent_path = os.path.dirname(pp[idx])

		EV_List = []
		cmd = "awk '{print $1}' %s > %s/mc1.1D" %(file, parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd = "awk '{print $2}' %s > %s/mc2.1D" %(file, parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd = "awk '{print $3}' %s > %s/mc3.1D" %(file, parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd = "awk '{print $4}' %s > %s/mc4.1D" %(file, parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd = "awk '{print $5}' %s > %s/mc5.1D" %(file, parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd = "awk '{print $6}' %s > %s/mc6.1D" %(file, parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		EV_List.append(os.path.abspath(parent_path+'/mc1.1D'))
		EV_List.append(os.path.abspath(parent_path+'/mc2.1D'))
		EV_List.append(os.path.abspath(parent_path+'/mc3.1D'))
		EV_List.append(os.path.abspath(parent_path+'/mc4.1D'))
		EV_List.append(os.path.abspath(parent_path+'/mc5.1D'))
		EV_List.append(os.path.abspath(parent_path+'/mc6.1D'))
		
		EV_Lists.append(EV_List)
		idx += 1

        return EV_Lists

def createFSF(nuisance_template,rest_pp,TR,n_vols):

	import os
	import sys
	import commands
	nuisance_files = []
	idx = 0
	for rest_p in rest_pp:
		parent_path = os.path.dirname(rest_p)

		cmd = "sed -e s:nuisance_dir:\"%s\":g <%s >%s/temp1" %(parent_path,nuisance_template,parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))


		cmd = "sed -e s:nuisance_model_outputdir:\"%s/residuals.feat\":g <%s/temp1 >%s/temp2" %(parent_path,parent_path,parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd = "sed -e s:nuisance_model_TR:\"%s\":g <%s/temp2 >%s/temp3" %(TR[idx],parent_path,parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd = "sed -e s:nuisance_model_numTRs:\"%s\":g <%s/temp3 >%s/temp4" %(n_vols[idx],parent_path,parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd = "sed -e s:nuisance_model_input_data:\"%s\":g <%s/temp4 >%s/nuisance.fsf" %(rest_p,parent_path,parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd = "rm %s/temp*" %(parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))
		nuisance_files.append(os.path.abspath("%s/nuisance.fsf"%(parent_path)))
		
		idx += 1


	return nuisance_files

def get_EV_LIST(EV_Lists, global1Ds, csf1Ds, wm1Ds):

	idx = 0
	EV_final_lists_set = []
	for EV_List in EV_Lists:
		print '>>>>><<<<', EV_List
		EV_List.append(global1Ds[idx])
		EV_List.append(csf1Ds[idx])
		EV_List.append(wm1Ds[idx])
		EV_final_lists_set.append(EV_List)
		
		idx += 1

	return EV_final_lists_set

def copyStuff(EV_Lists, nuisance_files, global1Ds,csf1Ds,wm1Ds):

	idx = 0
	import os
	import sys
	import commands

	global1ds = []
	csf1ds = []
	wm1ds = []
	EV_final_lists_set = []
	for nuisance_file in nuisance_files:

		parent_path = os.path.dirname(nuisance_file)

		cmd  = "cp %s %s/global.1D" %(global1Ds[idx],parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))
		
		cmd  = "cp %s %s/csf.1D" %(csf1Ds[idx],parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		cmd  = "cp %s %s/wm.1D" %(wm1Ds[idx],parent_path)
		print cmd
		sys.stderr.write(commands.getoutput(cmd))

		EV_Lists[idx].append("%s/global.1D" %(parent_path))
		EV_Lists[idx].append("%s/csf.1D" %(parent_path))
		EV_Lists[idx].append("%s/wm.1D" %(parent_path))
		EV_final_lists_set.append(EV_Lists[idx])
		idx += 1


	return EV_final_lists_set

def getStatsDir(in_files):

	import os

	stats_dir = []

	for in_file in in_files:

		parent_path = os.path.dirname(in_file)

		stats_dir.append(parent_path + '/stats')

	return stats_dir


def getOpString(mean, std_dev):

	str1 = "-sub %f -div %f" %(float(mean),float(std_dev))

	op_string = str1 + " -mas %s"

	return op_string

def pToFile(time_series,seed,out_dir):

	import os

	ts_oneD =  os.path.dirname(out_dir) + '/%s.1D'%( ( (os.path.basename(seed)).split('.') )[0] )
	print time_series
	f = open(ts_oneD,'w')

	for tr in time_series:
		print >>f,"%f" %(tr)
	
	f.close()
	
	return ts_oneD


RSFC_printToFile = pe.MapNode(util.Function(input_names = ['time_series','seed','out_dir'], output_names = ['ts_oneD'], function = pToFile), name='RSFC_printToFile', iterfield = ["time_series","out_dir"])

RSFC_printToFile.iterables = ("seed",seed_list) 

#TR and n_vols nodes

TR = pe.Node(util.Function(input_names = ['in_files'], output_names = ['TR'], function = getImgTR), name='TR')
NVOLS = pe.Node(util.Function(input_names = ['in_files'], output_names = ['nvols'], function = getImgNVols), name='NVOLS')

calcN1 = pe.MapNode(util.Function(input_names = ['nvols','TR','HP'], output_names = ['n1'], function = getN1), name='calcN1', iterfield = ["nvols","TR"])

calcN2 = pe.MapNode(util.Function(input_names = ['nvols','TR','LP','HP'], output_names = ['n2'], function = getN2), name='calcN2', iterfield = ["nvols","TR"])

alff_op_string = pe.MapNode(util.Function(input_names = ['mean','std_dev'], output_names = ['op_string'], function = getOpString), name='alff_op_string', iterfield = ["mean","std_dev"])

alff_op_string1 = pe.MapNode(util.Function(input_names = ['mean','std_dev'], output_names = ['op_string'], function = getOpString), name='alff_op_string1', iterfield = ["mean","std_dev"])


#create Motion Correction ID files for all six motion parameters
MC = pe.Node(util.Function(input_names = ['in_file','pp'], output_names = ['EV_Lists'], function = createMC), name='MC')


FSF = pe.Node(util.Function(input_names = ['nuisance_template','rest_pp','TR','n_vols'], output_names = ['nuisance_files'], function = createFSF), name='FSF')

copyS = pe.Node(util.Function(input_names = ['EV_Lists','nuisance_files','global1Ds','csf1Ds','wm1Ds'], output_names = ['EV_final_lists_set'], function = copyStuff), name='copyS')

EV_List = pe.Node(util.Function(input_names = ['EV_Lists','global1Ds','csf1Ds','wm1Ds'], output_names = ['EV_final_lists_set'], function = get_EV_LIST), name = 'EV_List')

getStats = pe.Node(util.Function(input_names = ['in_files'], output_names = ['stats_dir'], function = getStatsDir), name='getStats')


lifeSaver = pe.MapNode(interface = util.Rename(), name = 'lifeSaver', iterfield = ["in_file","format_string"])

#anatomical nodes


anat_refit = pe.Node(interface=e_afni.Threedrefit(), name='anat_refit')
anat_refit.inputs.deoblique = True

anat_reorient = pe.Node(interface=e_afni.Threedresample(), name='anat_reorient')
anat_reorient.inputs.orientation = 'RPI'

anat_skullstrip = pe.Node(interface= e_afni.ThreedSkullStrip(), name='anat_skullstrip')
anat_skullstrip.inputs.options = '-o_ply'

anat_calc = pe.Node(interface=e_afni.Threedcalc(), name='anat_calc')
anat_calc.inputs.expr = '\'a*step(b)\''
anat_calc.inputs.out_file = 'mprage_brain.nii.gz'

#funcpreproc Nodes

func_calc = pe.MapNode(interface = e_afni.Threedcalc(), name='func_calc', iterfield = ["infile_a","stop_idx"])
func_calc.inputs.start_idx = 0
#calc.inputs.infile_a = analysisdirectory + '/' + subject + '/' + 'rest_1/rest.nii.gz[%d..%d]' %(int(first_vol),int(last_vol))
func_calc.inputs.expr = '\'a\''
func_calc.inputs.out_file = 'rest_dr.nii.gz'

func_refit = pe.MapNode(interface=e_afni.Threedrefit(), name='func_refit', iterfield = ["in_file"])
func_refit.inputs.deoblique = True

func_reorient = pe.MapNode(interface=e_afni.Threedresample(), name='func_reorient', iterfield = ["in_file"])
func_reorient.inputs.orientation = 'RPI'
#reorient.inputs.infile = anat_name+'.nii.gz'
#reorient.inputs.out_file = os.path.abspath( analysisdirectory + '/' + subject + '/' + 'rest_1/rest_ro.nii.gz')

func_tstat = pe.MapNode(interface=e_afni.ThreedTstat(), name='func_tstat', iterfield = ["in_file"])
func_tstat.inputs.args = "-mean"
func_tstat.inputs.out_file = 'rest_ro_mean.nii.gz'

func_volreg = pe.MapNode(interface=e_afni.Threedvolreg(), name='func_volreg', iterfield = ["in_file","basefile"])
func_volreg.inputs.other = '-Fourier -twopass'
func_volreg.inputs.zpad = '4'
func_volreg.inputs.oned_file = 'rest_mc.1D'
func_volreg.inputs.out_file = 'rest_mc.nii.gz'

func_automask = pe.MapNode(interface=e_afni.ThreedAutomask(), name = 'func_automask', iterfield = ["in_file"])
func_automask.inputs.dilate = 1
#automask.inputs.out_file = os.path.abspath(analysisdirectory + '/' + subject + '/' + 'rest_1/rest_mask.nii.gz')


func_calcR = pe.MapNode(interface=e_afni.Threedcalc(), name='func_calcR', iterfield = ["infile_a","infile_b"])
#func_calcR.inputs.infile_a = analysisdirectory + '/' + subject + '/' + 'rest_1/rest_mc.nii.gz'
#func_calcR.inputs.infile_b = analysisdirectory + '/' + subject + '/' + 'rest_1/rest_mask.nii.gz'
func_calcR.inputs.expr = '\'a*b\''
func_calcR.inputs.out_file =  'rest_ss.nii.gz'

func_calcI = pe.MapNode(interface = e_afni.Threedcalc(), name='func_calcI', iterfield = ["infile_a"])
#func_calcI.inputs.infile_a = analysisdirectory + '/' + subject + '/' + 'rest_1/rest_ss.nii.gz[7]'
func_calcI.inputs.single_idx = 7
func_calcI.inputs.expr = '\'a\''
func_calcI.inputs.out_file = 'example_func.nii.gz'

func_despike = pe.MapNode(interface = e_afni.ThreedDespike(), name = 'func_despike', iterfield = ["in_file"])
func_despike.inputs.out_file = 'rest_ds.nii.gz'

func_smooth = pe.MapNode(interface = MultiImageMaths(), name = 'func_smooth', iterfield = ["in_file","operand_files"])
func_str1 = "-kernel gauss %f -fmean -mas" %(sigma)
func_smooth.inputs.op_string = func_str1 + " %s"
#func_smooth.inputs.out_file = os.path.abspath(analysisdirectory + '/' + subject + '/' + 'rest_1/rest_sm.nii.gz')

func_scale = pe.MapNode(interface = fsl.ImageMaths(), name = 'func_scale', iterfield = ["in_file"])
func_scale.inputs.op_string = '-ing 10000'
func_scale.inputs.out_data_type = 'float'
#func_scale.inputs.out_file = os.path.abspath(analysisdirectory + '/' + subject + '/' + 'rest_1/rest_gms.nii.gz')

func_filter = pe.MapNode(interface = e_afni.ThreedFourier(), name = 'func_filter', iterfield = ["in_file"])
#func_filter.inputs.highpass = str(hp)
#func_filter.inputs.lowpass = str(lp)
func_filter.inputs.other = '-retrend'
func_filter.inputs.out_file = 'rest_filt.nii.gz'

func_detrenda = pe.MapNode(interface = e_afni.ThreedTstat(), name = 'func_detrenda', iterfield = ["in_file"])
func_detrenda.inputs.options = '-mean'
func_detrenda.inputs.out_file = 'rest_filt_mean.nii.gz'

func_detrendb = pe.MapNode(interface = e_afni.ThreedDetrend(), name = 'func_detrendb', iterfield = ["in_file"])
func_detrendb.inputs.options = '-polort 2'
func_detrendb.inputs.out_file = 'rest_dt.nii.gz'


func_detrendc = pe.MapNode(interface = e_afni.Threedcalc(), name = 'func_detrendc', iterfield = ["infile_a","infile_b"])
func_detrendc.inputs.expr = '\'a+b\''
func_detrendc.inputs.out_file = 'rest_pp.nii.gz'

func_mask = pe.MapNode(interface = fsl.ImageMaths(), name = 'func_mask', iterfield = ["in_file"])
func_mask.inputs.op_string = '-Tmin -bin'
func_mask.inputs.out_data_type = 'char'
#func_mask.inputs.out_file = os.path.abspath(analysisdirectory + '/' + subject + '/' + 'rest_1/rest_pp_mask.nii.gz')



##registration

reg_flirt = pe.MapNode(interface= fsl.FLIRT(), name='reg_flirt', iterfield = ["in_file"])
#reg_flirt.inputs.in_file  = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#reg_flirt.inputs.reference = os.path.abspath(anat_reg_dir+'/highres.nii.gz')
#reg_flirt.inputs.out_file = 'example_func2highres.nii.gz'
#reg_flirt.inputs.out_matrix_file = os.path.abspath('example_func2highres.mat')
reg_flirt.inputs.cost = 'corratio'
reg_flirt.inputs.dof = 6
reg_flirt.inputs.interp = 'trilinear'

# Create mat file for conversion from subject's anatomical to functional
reg_xfm1 = pe.MapNode(interface= fsl.ConvertXFM(), name='reg_xfm1', iterfield = ["in_file"])
#reg_xfm1.inputs.out_file = 'highres2example_func.mat'
reg_xfm1.inputs.invert_xfm = True

## 4. T1->STANDARD
## NOTE THAT THIS IS Linear registration, you may want to use FNIRT (non-linear)
reg_flirt1 = pe.Node(interface=fsl.FLIRT(), name='reg_flirt1')
#reg_flirt1.inputs.in_file  = os.path.abspath(anat_reg_dir+'/highres.nii.gz')
#reg_flirt1.inputs.reference = os.path.abspath(anat_reg_dir+'/standard.nii.gz')
#reg_flirt1.inputs.out_file = 'highres2standard.nii.gz'
#reg_flirt1.inputs.out_matrix_file ='highres2standard.mat'
reg_flirt1.inputs.cost = 'corratio'
reg_flirt1.inputs.cost_func = 'corratio'
reg_flirt1.inputs.dof = 12
reg_flirt1.inputs.interp = 'trilinear'

## Create mat file for conversion from standard to high res
reg_xfm2 = pe.Node(interface= fsl.ConvertXFM(), name='reg_xfm2')
#reg_xfm2.inputs.in_file = os.path.abspath(reg_dir+'/highres2standard.mat')
#reg_xfm2.inputs.out_file = 'standard2highres.mat'
reg_xfm2.inputs.invert_xfm = True

## 5. FUNC->STANDARD
## Create mat file for registration of functional to standard
reg_xfm3 = pe.MapNode(interface= fsl.ConvertXFM(), name='reg_xfm3', iterfield = ["in_file"])
#reg_xfm3.inputs.in_file = os.path.abspath(reg_dir+'/example_func2highres.mat')
#reg_xfm3.inputs.in_file2 = os.path.abspath(reg_dir+'/highres2standard.mat')
#reg_xfm3.inputs.out_file = 'example_func2standard.mat'
reg_xfm3.inputs.concat_xfm = True

## apply registration
reg_flirt2 = pe.MapNode(interface=fsl.FLIRT(), name='reg_flirt2', iterfield = ["in_file","in_matrix_file"])
#reg_flirt2.inputs.in_file  = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#reg_flirt2.inputs.reference = os.path.abspath(anat_reg_dir+'/standard.nii.gz')
#reg_flirt2.inputs.out_file = 'example_func2standard.nii.gz'
#reg_flirt2.inputs.in_matrix_file = os.path.abspath(reg_dir+'/example_func2standard.mat')
reg_flirt2.inputs.apply_xfm = True
reg_flirt2.inputs.interp = 'trilinear'

## Create inverse mat file for registration of standard to functional
reg_xfm4 = pe.MapNode(interface= fsl.ConvertXFM(), name='reg_xfm4', iterfield = ["in_file"])
#reg_xfm4.inputs.in_file = os.path.abspath(reg_dir+'/example_func2standard.mat')
#reg_xfm4.inputs.out_file = 'standard2example_func.mat'
reg_xfm4.inputs.invert_xfm = True

## 6. T1->STANDARD NONLINEAR
# Perform nonlinear registration (higres to standard)
reg_fnt = pe.Node(interface = fsl.FNIRT(), name = 'reg_fnt')
#reg_fnt.inputs.in_file = os.path.abspath(anat_reg_dir+'/highres_head.nii.gz')
#reg_fnt.inputs.affine_file = os.path.abspath(anat_reg_dir+'/highres2standard.mat')
reg_fnt.inputs.fieldcoeff_file = True
#reg_fnt.inputs.warped_file = os.path.abspath(anat_reg_dir+'/highres2standard_NL.nii.gz')
reg_fnt.inputs.jacobian_file = True
reg_fnt.inputs.fieldcoeff_file = True
#reg_fnt.inputs.config_file = os.path.abspath(FSLDIR+'/etc/flirtsch/T1_2_MNI152_%s.cnf'%(standard_res) )
#reg_fnt.inputs.ref_file = os.path.abspath(FSLDIR+'/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#reg_fnt.inputs.refmask_file = os.path.abspath(FSLDIR+'/data/standard/MNI152_T1_%s_brain_mask_dil.nii.gz' %(standard_res))
reg_fnt.inputs.warp_resolution = (10,10,10)

reg_warp = pe.MapNode(interface = fsl.ApplyWarp(), name = 'reg_warp', iterfield = ["in_file","premat"])
#reg_warp.inputs.in_file = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#reg_warp.inputs.ref_file = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#reg_warp.inputs.field_file = 'highres2standard_warp.nii.gz'
#reg_warp.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#reg_warp.inputs.out_file = os.path.abspath(func_reg_dir + '/example_func2standard_NL.nii.gz')


#segmentation

seg_segment = pe.Node(interface= fsl.FAST(), name = 'seg_segment')
#seg_segment.inputs.in_files = os.path.abspath(anat_dir+'/mprage_brain.nii.gz')
seg_segment.inputs.img_type = 1
seg_segment.inputs.segments = True
seg_segment.inputs.probability_maps = True
seg_segment.inputs.out_basename = 'segment'

## 4. Copy functional mask from FSLpreproc step 5 - this is the global signal mask
seg_copy = pe.MapNode(interface= e_afni.Threedcopy(), name='seg_copy', iterfield = ["in_file"])
#seg_copy.inputs.in_file = os.path.abspath(func_dir+'/rest_pp_mask.nii.gz')
#seg_copy.inputs.out_file = 'global_mask.nii.gz'

## CSF
## 5. Register csf to native space
seg_flirt = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt', iterfield = ["reference","in_matrix_file"])
#seg_flirt.inputs.in_file  = os.path.abspath(anat_dir+'/segment_prob_0.nii.gz')
#seg_flirt.inputs.reference = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#seg_flirt.inputs.out_file = os.path.abspath(segment_dir+'/csf2func.nii.gz')
#seg_flirt.inputs.in_matrix_file = os.path.abspath(func_reg_dir+'/highres2example_func.mat')
seg_flirt.inputs.apply_xfm = True

## 6. Smooth image to match smoothing on functional
seg_smooth = pe.MapNode(interface = fsl.ImageMaths(), name = 'seg_smooth', iterfield = ["in_file"])
seg_str1 = "-kernel gauss %f -fmean " %(sigma)
#seg_smooth.inputs.in_file = os.path.abspath(segment_dir+'/csf2func.nii.gz')
seg_smooth.inputs.op_string = seg_str1
#seg_smooth.inputs.out_file = os.path.abspath(segment_dir+'/csf_sm.nii.gz')

## 7. register to standard
seg_flirt1 = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt1', iterfield = ["in_file","in_matrix_file"])
#seg_flirt1.inputs.in_file  = os.path.abspath(segment_dir+'/csf_sm.nii.gz')
#seg_flirt1.inputs.reference = os.path.abspath(anat_reg_dir+'/standard.nii.gz')
#seg_flirt1.inputs.out_file = os.path.abspath(segment_dir+'/csf2standard.nii.gz')
#seg_flirt1.inputs.in_matrix_file = os.path.abspath(func_reg_dir+'/example_func2standard.mat')
seg_flirt1.inputs.apply_xfm = True

## 8. find overlap with prior
seg_smooth1 = pe.MapNode(interface = MultiImageMaths(), name = 'seg_smooth1', iterfield = ["in_file"])
seg_str1 = "-mas %s "
#seg_smooth1.inputs.in_file = os.path.abspath(segment_dir+'/csf2standard.nii.gz')
seg_smooth1.inputs.op_string = seg_str1
#seg_smooth1.inputs.operand_files = os.path.abspath(PRIOR_CSF)
#seg_smooth1.inputs.out_file = os.path.abspath(segment_dir+'/csf_masked.nii.gz')

## 9. revert back to functional space
seg_flirt2 = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt2', iterfield = ["in_file", "reference","in_matrix_file"])
#seg_flirt2.inputs.in_file  = os.path.abspath(segment_dir+'/csf_masked.nii.gz')
#seg_flirt2.inputs.reference = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#seg_flirt2.inputs.out_file = os.path.abspath(segment_dir+'/csf_native.nii.gz')
#seg_flirt2.inputs.in_matrix_file = os.path.abspath(func_reg_dir+'/standard2example_func.mat')
seg_flirt2.inputs.apply_xfm = True

## 10. Threshold and binarize probability map of csf
seg_thresh = pe.MapNode(interface= fsl.ImageMaths(), name = 'seg_thresh', iterfield = ["in_file"])
seg_str1 = "-thr 0.4 -bin "
#seg_thresh.inputs.in_file = os.path.abspath(segment_dir+'/csf_native.nii.gz')
seg_thresh.inputs.op_string = seg_str1
#seg_thresh.inputs.out_file = os.path.abspath(segment_dir+'/csf_bin.nii.gz')

## 11. Mask again by the subject's functional
seg_mask = pe.MapNode(interface = fsl.MultiImageMaths(), name = 'seg_mask', iterfield = ["in_file","operand_files"])
seg_str1 = "-mas %s "
#seg_mask.inputs.in_file = os.path.abspath(segment_dir+'/csf_bin.nii.gz')
seg_mask.inputs.op_string = seg_str1
#seg_mask.inputs.operand_files = os.path.abspath(segment_dir+'global_mask.nii.gz')
#seg_mask.inputs.out_file = os.path.abspath(segment_dir+'/csf_mask.nii.gz')

## WM
## 12. Register wm to native space
seg_flirt3 = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt3', iterfield = ["reference","in_matrix_file"])
#seg_flirt3.inputs.in_file  = os.path.abspath(anat_dir+'/segment_prob_2.nii.gz')
#seg_flirt3.inputs.reference = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#seg_flirt3.inputs.out_file = os.path.abspath(segment_dir+'/wm2func.nii.gz')
#seg_flirt3.inputs.in_matrix_file = os.path.abspath(func_reg_dir+'/highres2example_func.mat')
seg_flirt3.inputs.apply_xfm = True

## 13. Smooth image to match smoothing on functional
seg_smooth2 = pe.MapNode(interface= fsl.ImageMaths(), name = 'seg_smooth2', iterfield = ["in_file"])
seg_str1 = "-kernel gauss %f -fmean " %(sigma)
#seg_smooth2.inputs.in_file = os.path.abspath(segment_dir+'/wm2func.nii.gz')
seg_smooth2.inputs.op_string = seg_str1
#seg_smooth2.inputs.out_file = os.path.abspath(segment_dir+'/wm_sm.nii.gz')

## 14. register to standard
seg_flirt4 = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt4', iterfield = ["in_file","in_matrix_file"])
#seg_flirt4.inputs.in_file  = os.path.abspath(segment_dir+'/wm_sm.nii.gz')
#seg_flirt4.inputs.reference = os.path.abspath(anat_reg_dir+'/standard.nii.gz')
#seg_flirt4.inputs.out_file = os.path.abspath(segment_dir+'/wm2standard.nii.gz')
#seg_flirt4.inputs.in_matrix_file = os.path.abspath(func_reg_dir+'/example_func2standard.mat')
seg_flirt4.inputs.apply_xfm = True

## 15. find overlap with prior
seg_prior1 = pe.MapNode(interface = fsl.MultiImageMaths(), name = 'seg_prior1', iterfield = ["in_file"])
seg_str1 = "-mas %s "
#seg_prior1.inputs.in_file = os.path.abspath(segment_dir+'/wm2standard.nii.gz')
seg_prior1.inputs.op_string = seg_str1
#seg_prior1.inputs.operand_files = os.path.abspath(PRIOR_WHITE)
#seg_prior1.inputs.out_file = os.path.abspath(segment_dir+'/wm_masked.nii.gz')

## 16. revert back to functional space
seg_flirt5 = pe.MapNode(interface=fsl.FLIRT(), name='seg_flirt5', iterfield = ["in_file","reference","in_matrix_file"])
#seg_flirt5.inputs.in_file  = os.path.abspath(segment_dir+'/wm_masked.nii.gz')
#seg_flirt5.inputs.reference = os.path.abspath(func_reg_dir+'/example_func.nii.gz')
#seg_flirt5.inputs.out_file = os.path.abspath(segment_dir+'/wm_native.nii.gz')
#seg_flirt5.inputs.in_matrix_file = os.path.abspath(func_reg_dir+'/standard2example_func.mat')
seg_flirt5.inputs.apply_xfm = True

## 17. Threshold and binarize probability map of wm
seg_thresh1 = pe.MapNode(interface= fsl.ImageMaths(), name = 'seg_thresh1', iterfield = ["in_file"])
seg_str1 = "-thr 0.66 -bin "
#seg_thresh1.inputs.in_file = os.path.abspath(segment_dir+'/wm_native.nii.gz')
seg_thresh1.inputs.op_string = seg_str1
#seg_thresh1.inputs.out_file = os.path.abspath(segment_dir+'/wm_bin.nii.gz')

## 18. Mask again by the subject's functional
seg_mask1 = pe.MapNode(interface = fsl.MultiImageMaths(), name = 'seg_mask1', iterfield = ["in_file","operand_files"])
seg_str1 = "-mas %s "
#seg_mask1.inputs.in_file = os.path.abspath(segment_dir+'/wm_bin.nii.gz')
seg_mask1.inputs.op_string = seg_str1
#seg_mask1.inputs.operand_files = os.path.abspath(segment_dir+'global_mask.nii.gz')
#seg_mask1.inputs.out_file = os.path.abspath(segment_dir+'/wm_mask.nii.gz')


nuisance_globalE = pe.MapNode(interface=e_afni.ThreedMaskave(), name = 'nuisance_globalE', iterfield = ["in_file","mask"])
#nuisance_globalE.inputs.in_file = os.path.abspath(func_dir + '/' + 'rest_pp.nii.gz')
#nuisance_globalE.inputs.mask = os.path.abspath(segment_dir + '/' + 'global_mask.nii.gz')
nuisance_globalE.inputs.quiet = True
nuisance_globalE.inputs.out_file = 'global.1D'

## 4. csf
nuisance_csf = pe.MapNode(interface=e_afni.ThreedMaskave(), name = 'nuisance_csf', iterfield = ["in_file","mask"])
#nuisance_csf.inputs.in_file = os.path.abspath(func_dir + '/' + 'rest_pp.nii.gz')
#nuisance_csf.inputs.mask = os.path.abspath(segment_dir + '/' + 'csf_mask.nii.gz')
nuisance_csf.inputs.quiet = True
nuisance_csf.inputs.out_file =  'csf.1D'

## 5. wm
nuisance_wm = pe.MapNode(interface=e_afni.ThreedMaskave(), name = 'nuisance_wm', iterfield = ["in_file","mask"])
#nuisance_wm.inputs.in_file = os.path.abspath(func_dir + '/' + 'rest_pp.nii.gz')
#nuisance_wm.inputs.mask = os.path.abspath(segment_dir + '/' + 'wm_mask.nii.gz')
nuisance_wm.inputs.quiet = True
nuisance_wm.inputs.out_file = 'wm.1D'

## feat model
nuisance_featM = pe.MapNode(interface = fsl.FEATModel(), name = 'nuisance_featM', iterfield = ["fsf_file","ev_files"])
#nuisance_featM.inputs.fsf_file = os.path.abspath(nuisance_dir + '/nuisance.fsf')
#nuisance_featM.inputs.ev_files = EV_List
#nuisance_featM.inputs.design_file = os.path.abspath(nuisance_dir + '/nuisance.mat')

nuisance_brick = pe.MapNode(interface = e_afni.ThreedBrickStat(), name = 'nuisance_brick', iterfield = ["in_file","mask"])
#nuisance_brick.inputs.in_file = os.path.abspath(func_dir + '/rest_pp.nii.gz')
#nuisance_brick.inputs.mask = os.path.abspath(func_dir + '/rest_pp_mask.nii.gz')
nuisance_brick.inputs.min = True

## 7. Get residuals
nuisance_fgls = pe.MapNode(interface = fsl.FILMGLS() , name = 'nuisance_fgls', iterfield = ["in_file","results_dir","design_file","threshold"])
#nuisance_fgls.inputs.in_file = os.path.abspath(func_dir + '/rest_pp.nii.gz')
#nuisance_fgls.inputs.results_dir = nuisance_dir + '/stats'
nuisance_fgls.inputs.mask_size = 5
nuisance_fgls.inputs.smooth_autocorr = True
nuisance_fgls.inputs.autocorr_noestimate = True
#nuisance_fgls.inputs.design_file = os.path.abspath(nuisance_dir + '/nuisance.mat')

## 8. Demeaning residuals and ADDING 100
nuisance_stat = pe.MapNode(interface=e_afni.ThreedTstat(), name = 'nuisance_stat', iterfield = ["in_file"])
#nuisance_stat.inputs.in_file = os.path.abspath(nuisance_dir + '/stats/res4d.nii.gz')
nuisance_stat.inputs.options = ' -mean '
nuisance_stat.inputs.out_file = 'res4d_mean.nii.gz'

nuisance_calc = pe.MapNode(interface = e_afni.Threedcalc(), name = 'nuisance_calc', iterfield = ["infile_a","infile_b"])
#nuisance_calc.inputs.infile_a = os.path.abspath(nuisance_dir + '/stats/res4d.nii.gz')
#nuisance_calc.inputs.infile_b = os.path.abspath(nuisance_dir + '/stats/res4d_mean.nii.gz')
nuisance_calc.inputs.expr = '\'(a-b)+100\''
nuisance_calc.inputs.out_file = 'rest_res.nii.gz'

## 9. Resampling residuals to MNI space
nuisance_warp = pe.MapNode(interface = fsl.ApplyWarp(), name = 'nuisance_warp', iterfield = ["in_file","premat"])
#nuisance_warp.inputs.in_file = os.path.abspath(func_dir+'/rest_res.nii.gz')
#nuisance_warp.inputs.ref_file = os.path.abspath( FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res) )
#nuisance_warp.inputs.field_file = os.path.abspath(anat_reg_dir + '/highres2standard_warp.nii.gz')
#nuisance_warp.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#nuisance_warp.inputs.out_file = os.path.abspath(func_dir + '/rest_res2standard.nii.gz')


alff_detrend = pe.MapNode(interface = e_afni.ThreedTcat(), name = 'alff_detrend', iterfield = ["in_file"])
#alff_detrend.inputs.in_file = os.path.abspath(func_dir + '/rest_ds.nii.gz')
alff_detrend.inputs.out_file = 'prefiltered_func_data_rlt.nii.gz'
alff_detrend.inputs.rlt = '+'

## 3. Spatial Smoothing
alff_smooth = pe.MapNode(interface = MultiImageMaths(), name = 'alff_smooth', iterfield = ["in_file","operand_files"])
alff_str1 = "-kernel gauss %f -fmean -mas " %(sigma)
alff_smooth.inputs.op_string = alff_str1 + " %s"
#alff_smooth.inputs.in_file = os.path.abspath(alff_dir+'/prefiltered_func_data_rlt.nii.gz')
#alff_smooth.inputs.operand_files = os.path.abspath(func_dir + '/rest_mask.nii.gz')
#alff_smooth.inputs.out_file = os.path.abspath(alff_dir+'/prefiltered_func_data_smooth.nii.gz')

## 4. Grand mean Scaling
alff_scale = pe.MapNode(interface = fsl.ImageMaths(), name = 'alff_scale', iterfield = ["in_file"])
alff_scale.inputs.op_string = '-ing 10000'
#alff_scale.inputs.in_file = os.path.abspath(alff_dir+'/prefiltered_func_data_smooth.nii.gz')
#alff_scale.inputs.out_file = os.path.abspath(alff_dir+'/prefiltered_func_data_inn_volsorm.nii.gz')

alff_cp = pe.MapNode(interface = fsl.ImageMaths(), name = 'alff_cp', iterfield = ["in_file"])
#alff_cp.inputs.in_file = os.path.abspath(alff_dir+'/prefiltered_func_data_inn_volsorm.nii.gz')
#alff_cp.inputs.out_file = os.path.abspath(alff_dir+'/rest_alff_pp.nii.gz')

alff_mean = pe.MapNode(interface = fsl.ImageMaths(), name = 'alff_mean', iterfield = ["in_file"])
#alff_mean.inputs.in_file = os.path.abspath(alff_dir+'/rest_alff_pp.nii.gz')
alff_mean.inputs.op_string = '-Tmean'
#alff_mean.inputs.out_file = os.path.abspath(alff_dir+'/mean_rest_alff_pp.nii.gz')

## 4. Calculate fALFF
alff_falff = pe.MapNode(interface = fsl.ImageMaths(), name = 'alff_falff', iterfield = ["in_file","op_string"])
#alff_falff.inputs.op_string = '-Tmean -mul %s -div 2' %(n_vols)
#alff_falff.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sqrt.nii.gz')
#alff_falff.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sum.nii.gz')

alff_falff1 = pe.MapNode(interface = MultiImageMaths(), name = 'falff1', iterfield = ["in_file","operand_files"])
alff_falff1.inputs.op_string = "-div %s"
#alff_falff1.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_ps_alff4slow.nii.gz')
#alff_falff1.inputs.operand_files = os.path.abspath(alff_dir+'/prealff_func_data_ps_sum.nii.gz')
#alff_falff1.inputs.out_file = os.path.abspath(alff_dir+'/fALFF.nii.gz')

## 5. Z-normalisation across whole brain
alff_normM = pe.MapNode(interface = fsl.ImageStats(), name = 'alff_normM', iterfield = ["in_file","mask_file"])
#alff_normM.inputs.in_file = os.path.abspath(alff_dir+'/ALFF.nii.gz')
#alff_normM.inputs.mask_file = os.path.abspath(func_dir+'/rest_mask.nii.gz')
alff_normM.inputs.op_string = "-k %s -m"

alff_normS = pe.MapNode(interface = fsl.ImageStats(), name = 'alff_normS', iterfield = ["in_file","mask_file"])
#alff_normS.inputs.in_file = os.path.abspath(alff_dir+'/ALFF.nii.gz')
#alff_normS.inputs.mask_file = os.path.abspath(func_dir+'/rest_mask.nii.gz')
alff_normS.inputs.op_string = "-k %s -s"

alff_normM1 = pe.MapNode(interface = fsl.ImageStats(), name = 'alff_normM1', iterfield = ["in_file","mask_file"])
#alff_normM1.inputs.in_file = os.path.abspath(alff_dir+'/fALFF.nii.gz')
#alff_normM1.inputs.mask_file = os.path.abspath(func_dir+'/rest_mask.nii.gz')
alff_normM1.inputs.op_string = "-k %s -m"

alff_normS1 = pe.MapNode(interface = fsl.ImageStats(), name = 'alff_normS1', iterfield = ["in_file","mask_file"])
#alff_normS1.inputs.in_file = os.path.abspath(alff_dir+'/fALFF.nii.gz')
#alff_normS1.inputs.mask_file = os.path.abspath(func_dir+'/rest_mask.nii.gz')
alff_normS1.inputs.op_string = "-k %s -s"

alff_Z_alff = pe.MapNode(interface = MultiImageMaths(), name = 'alff_Z_alff', iterfield = ["in_file","operand_files","op_string"])
#alff_str = "-sub %f -div %f" %(mean,std_dev)
#alff_str1 = str + " -mas %s "
#alff_Z_alff.inputs.in_file = os.path.abspath(alff_dir+'/ALFF.nii.gz')
#alff_Z_alff.inputs.op_string = alff_str1
#alff_Z_alff.inputs.operand_files = os.path.abspath(func_dir+'/rest_mask.nii.gz')
#alff_Z_alff.inputs.out_file = os.path.abspath(alff_dir+'/ALFF_Z.nii.gz')

alff_Z_falff = pe.MapNode(interface = MultiImageMaths(), name = 'alff_Z_falff', iterfield = ["in_file","operand_files","op_string"])
#alff_str = "-sub %f -div %f" %(mean1,std_dev1)
#alff_str1 = str + " -mas %s "
#alff_Z_falff.inputs.in_file = os.path.abspath(alff_dir+'/fALFF.nii.gz')
#alff_Z_falff.inputs.op_string = alff_str1

#alff_Z_falff.inputs.operand_files = os.path.abspath(func_dir+'/rest_mask.nii.gz')
#alff_Z_falff.inputs.out_file = os.path.abspath(alff_dir+'/fALFF_Z.nii.gz')


#Registering Z-transformed ALFF to standard space
alff_warp_alff = pe.MapNode(interface = fsl.ApplyWarp(), name = 'alff_warp_alff', iterfield = ["in_file","premat"])
#alff_warp_alff.inputs.in_file = os.path.abspath(alff_dir+'/ALFF_Z.nii.gz')
#alff_warp_alff.inputs.ref_file = os.path.abspath( FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#alff_warp_alff.inputs.field_file = os.path.abspath(anat_reg_dir + '/highres2standard_warp.nii.gz')
#alff_warp_alff.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#alff_warp_alff.inputs.out_file = os.path.abspath(alff_dir + '/ALFF_Z_2standard.nii.gz')

#Registering Z-transformed fALFF to standard space
alff_warp_falff = pe.MapNode(interface = fsl.ApplyWarp(), name = 'alff_warp_falff', iterfield = ["in_file","premat"])
#alff_warp_falff.inputs.in_file = os.path.abspath(alff_dir+'/fALFF_Z.nii.gz')
#alff_warp_falff.inputs.ref_file = os.path.abspath( FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#alff_warp_falff.inputs.field_file = os.path.abspath(anat_reg_dir + '/highres2standard_warp.nii.gz')
#alff_warp_falff.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#alff_warp_falff.inputs.out_file = os.path.abspath(alff_dir + '/fALFF_Z_2standard.nii.gz')

alff_roi = pe.MapNode(interface = fsl.ExtractROI(), name = "alff_roi", iterfield = ["in_file","t_size"])
#alff_roi.inputs.in_file = os.path.abspath(alff_dir+'/rest_alff_pp.nii.gz')
#alff_roi.inputs.roi_file = os.path.abspath(alff_dir+'/prealff_func_data.nii.gz')
alff_roi.inputs.t_min = 1
#alff_roi.inputs.t_size = int(n_vols)

alff_cp1 = pe.MapNode(interface = fsl.ImageMaths(), name = 'alff_cp1', iterfield = ["in_file"])
#alff_cp1.inputs.in_file = os.path.abspath(alff_dir+'/rest_alff_pp.nii.gz')
#alff_cp1.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_data.nii.gz')

alff_concatnode = pe.MapNode(interface=util.Merge(2), name = 'alff_concatnode', iterfield = ["in1","in2"])

alff_selectnode = pe.MapNode(interface=util.Select(),name='alff_selectnode', iterfield = ["inlist","index"])

alff_pspec = pe.MapNode(interface = fsl.PowerSpectrum(), name = 'alff_pspec', iterfield = ["in_file"])
#alff_pspec.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_data.nii.gz')
#alff_pspec.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_data_ps.nii.gz')

##compute sqrt of power spectrum
alff_sqrt = pe.MapNode(interface = fsl.ImageMaths(), name = 'alff_sqrt', iterfield = ["in_file"])
alff_sqrt.inputs.op_string = '-sqrt'
#alff_sqrt.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_data_ps.nii.gz')
#alff_sqrt.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sqrt.nii.gz')

alff_roi1 = pe.MapNode(interface = fsl.ExtractROI(), name = "alff_roi1", iterfield = ["in_file","t_min","t_size"])
#alff_roi1.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_data_ps_sqrt.nii.gz')
#alff_roi1.inputs.roi_file = os.path.abspath(alff_dir+'/prealff_func_ps_slow.nii.gz')
#alff_roi1.inputs.t_min = int(n1)
#alff_roi1.inputs.t_size = int(n2)

## calculate ALFF as the sum of the amplitudes in the low frequency band
alff_sum = pe.MapNode(interface = fsl.ImageMaths(), name = 'alff_sum',iterfield = ["in_file","op_string"])
#alff_sum.inputs.op_string = '-Tmean -mul %f' %(n2)
#alff_sum.inputs.in_file = os.path.abspath(alff_dir+'/prealff_func_ps_slow.nii.gz')
#alff_sum.inputs.out_file = os.path.abspath(alff_dir+'/prealff_func_ps_alff4slow.nii.gz')


## 1. Extract Timeseries
RSFC_time_series = pe.MapNode(interface = e_afni.ThreedROIstats(), name = 'RSFC_time_series', iterfield = ["in_file"])
#RSFC_time_series.inputs.in_file = os.path.abspath(func_dir+'/rest_res2standard.nii.gz')
#RSFC_time_series.inputs.mask = os.path.abspath(seed)
RSFC_time_series.inputs.quiet = True
RSFC_time_series.inputs.mask_f2short = True
RSFC_time_series.iterables = ("mask",seed_list)
#RSFC_(time_series.run()).outputs.stats

## 2. Compute voxel-wise correlation with Seed Timeseries
RSFC_corr = pe.MapNode(interface = e_afni.Threedfim(), name = 'RSFC_corr', iterfield = ["in_file","ideal_file"])
#RSFC_corr.inputs.in_file = os.path.abspath(func_dir + '/rest_res.nii.gz')
#RSFC_corr.inputs.ideal_file = os.path.abspath(seed_ts_dir + '/%s.1D'%(seed_name))
RSFC_corr.inputs.fim_thr = 0.0009
RSFC_corr.inputs.out = 'Correlation'
RSFC_corr.inputs.out_file = 's_corr.nii.gz'

## 3. Z-transform correlations
RSFC_z_trans = pe.MapNode(interface = e_afni.Threedcalc(), name = 'RSFC_z_trans', iterfield = ["infile_a"])
#RSFC_z_trans.inputs.infile_a = os.path.abspath(RSFC_dir+'/%s_corr.nii.gz'%(seed_name))
RSFC_z_trans.inputs.expr = '\'log((a+1)/(a-1))/2\''
RSFC_z_trans.inputs.out_file = 's_Z.nii.gz'

## 4. Register Z-transformed correlations to standard space
RSFC_register = pe.MapNode(interface = fsl.ApplyWarp(), name = 'RSFC_register', iterfield = ["premat","in_file"])
#RSFC_register.inputs.in_file = os.path.abspath(RSFC_dir+'/%s_Z.nii.gz'%(seed_name))
#RSFC_register.inputs.ref_file = os.path.abspath(FSLDIR + '/data/standard/MNI152_T1_%s.nii.gz' %(standard_res))
#RSFC_register.inputs.field_file = os.path.abspath(anat_reg_dir + '/highres2standard_warp.nii.gz')
#RSFC_register.inputs.premat = os.path.abspath(func_reg_dir + '/example_func2highres.mat')
#RSFC_register.inputs.out_file = os.path.abspath(RSFC_dir + '/%s_Z_2standard.nii.gz' %(seed_name))

## Linear registration of T1 --> symmetric standard
VMHC_flirt = pe.Node(interface= fsl.FLIRT(), name='VMHC_flirt')
#VMHC_flirt.inputs.in_file  = os.path.abspath(${anat_dir}/brain.nii.gz)
#VMHC_flirt.inputs.reference = os.path.abspath(${symm_standard_brain})
#VMHC_flirt.inputs.out_file = '${anat_reg_dir}/highres2symmstandard.nii.gz'
#VMHC_flirt.inputs.out_matrix_file = os.path.abspath('${anat_reg_dir}/highres2symmstandard.mat')
VMHC_flirt.inputs.cost = 'corratio'
VMHC_flirt.inputs.cost_func = 'corratio'
VMHC_flirt.inputs.dof = 12
VMHC_flirt.inputs.interp = 'trilinear'

## Perform nonlinear registration (higres to standard) to symmetric standard brain
VMHC_fnt = pe.Node(interface = fsl.FNIRT(), name = 'VMHC_fnt')
#VMHC_fnt.inputs.in_file = os.path.abspath('${anat_dir}/head.nii.gz')
#VMHC_fnt.inputs.affine_file = os.path.abspath('${anat_reg_dir}/highres2symmstandard.mat')
VMHC_fnt.inputs.fieldcoeff_file = True
#VMHC_fnt.inputs.warped_file = os.path.abspath(${anat_reg_dir}/highres2symmstandard_warp.nii.gz)
#VMHC_fnt.inputs.jacobian_file = os.path.abspath(${anat_reg_dir}/highres2symmstandard_jac)
#VMHC_fnt.inputs.config_file = os.path.abspath(${T1_2_MNI152_2mm_symmetric} )
#VMHC_fnt.inputs.ref_file = os.path.abspath(${symm_standard})
#VMHC_fnt.inputs.refmask_file = os.path.abspath(${FSLDIR}/data/standard/MNI152_T1_2mm_brain_mask_symmetric_dil.nii.gz)
VMHC_fnt.inputs.warp_resolution = (10,10,10)

## Apply nonlinear registration (func to standard)
VMHC_warp = pe.MapNode(interface = fsl.ApplyWarp(), name = 'VMHC_warp', iterfield = ["in_file","premat"])
#VMHC_warp.inputs.in_file = os.path.abspath(${func_dir}/${image2use}.nii.gz)
#VMHC_warp.inputs.ref_file = os.path.abspath(${symm_standard})
#VMHC_warp.inputs.field_file = '${anat_reg_dir}/highres2symmstandard_warp.nii.gz'
#VMHC_warp.inputs.premat = os.path.abspath(${reg_dir}/example_func2highres.mat)
#VMHC_warp.inputs.out_file = os.path.abspath('${vmhc_dir}/${image2use}2symmstandard.nii.gz')


## copy and L/R swap file
VMHC_swap = pe.MapNode(interface = fsl.SwapDimensions(), name = 'VMHC_swap', iterfield = ["in_file"])
#VMHC_swap.inputs.in_file = '${vmhc_dir}/${image2use}.nii.gz'
VMHC_swap.inputs.new_dims = ('-x','y','z')
VMHC_swap.inputs.out_file = 'tmp_LRflipped.nii.gz'

## caculate vmhc
VMHC_corr = pe.MapNode(interface = e_afni.ThreedTcorrelate(), name = 'VMHC_corr', iterfield = ["xset"])
VMHC_corr.inputs.pearson = True
VMHC_corr.inputs.polort = -1
VMHC_corr.inputs.out_file = 'VMHC.nii.gz'
#VMHC_corr.inputs.xset = '${vmhc_dir}/${image2use}.nii.gz'
#VMHC_corr.inputs.yset = 'tmp_LRflipped.nii.gz'

## Fisher Z transform map

VMHC_z_trans = pe.MapNode(interface = e_afni.Threedcalc(), name = 'VMHC_z_trans', iterfield = ["infile_a"])
#VMHC_z_trans.inputs.infile_a = os.path.abspath('VMHC.nii.gz')
VMHC_z_trans.inputs.expr = '\'log((a+1)/(a-1))/2\''
VMHC_z_trans.inputs.out_file = 'VMHC_Z.nii.gz'

VMHC_z_stat = pe.MapNode(interface = e_afni.Threedcalc(), name = 'VMHC_z_stat', iterfield = ["infile_a","expr"])
VMHC_z_stat.inputs.out_file = 'VMHC_z_stat.nii.gz'


