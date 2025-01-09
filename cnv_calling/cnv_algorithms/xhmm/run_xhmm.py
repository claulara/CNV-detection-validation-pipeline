'''
Created on October 2020

@author: adelpozo
'''

#!/usr/bin/python

import os,sys,gzip,re,shutil

import ConfigParser, optparse

from subprocess import Popen , PIPE

import logging

from operator import itemgetter

from itertools import groupby

from collections import OrderedDict

from pybedtools import BedTool

import numpy as np

import pandas as pd

import glob

######################################################################

class OptionParser(optparse.OptionParser):

	def check_required (self, opt):

		option = self.get_option(opt)

		atrib = getattr(self.values, option.dest)
		
		if atrib is None:
#			self.error("%s option not supplied" % option)
			return False
		else:
			return True

##########################################################################

def read_cfg_file(cfg_filename):
	
	hash_cfg = {}
	
	fi = open(cfg_filename,'r')

	config = ConfigParser.ConfigParser()
	config.readfp(fi)

	for field in config.options('MARKERS'):
		hash_cfg[field] = config.get('MARKERS',field)

	for field in config.options('SOFTWARE'):
		hash_cfg[field] = config.get('SOFTWARE',field)

	fi.close()

	return hash_cfg 

##########################################################################

def configure_process(xhmm_path):
	
	hash_parameters = {}
	
	version = ""
	
	try:
		xhmm_sal = Popen([xhmm_path,'--version'],stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(output,logdata) = xhmm_sal.communicate()
		xhmm_sal.wait()
		
		version = output
	except:
		pass
		
	hash_parameters['version'] = version
	hash_parameters['Q_threshold']   = 15
	hash_parameters['minTargetSize'] = 10
	hash_parameters['maxTargetSize'] = 10000
	hash_parameters['minMeanTargetRD'] = 10
	hash_parameters['maxMeanTargetRD'] = 2000
	hash_parameters['minMeanSampleRD'] = 10
	hash_parameters['maxMeanSampleRD'] = 2000
	hash_parameters['minSdSampleRD'] = 0
	hash_parameters['maxSdSampleRD'] = 500
	
	return hash_parameters

##########################################################################

def write_configuration(hash_param,output_path):
	
	output_filename = os.path.join(output_path,'xhmm_param.cfg')
	
	f = open(output_filename,'w')
	
	f.write('xhmm version:\t%s\n' % (hash_param['version']))
	f.write('mosdepth quality threshold:\t%d\n' % (hash_param['Q_threshold']))
	f.write('xhmm minTargetSize:\t%d\n' % (hash_param['minTargetSize']))
	f.write('xhmm maxTargetSize:\t%d\n' % (hash_param['maxTargetSize']))
	f.write('xhmm minMeanTargetRD:\t%d\n' % (hash_param['minMeanTargetRD']))
	f.write('xhmm maxMeanTargetRD:\t%d\n' % (hash_param['maxMeanTargetRD']))
	f.write('xhmm minMeanSampleRD:\t%d\n' % (hash_param['minMeanSampleRD']))
	f.write('xhmm maxMeanSampleRD:\t%d\n' % (hash_param['maxMeanSampleRD']))
	f.write('xhmm minSdSampleRD:\t%d\n' % (hash_param['minSdSampleRD']))
	f.write('xhmm maxSdSampleRD:\t%d\n\n' % (hash_param['maxSdSampleRD']))
	
	f.close()
	
	return output_filename

##########################################################################

def write_excluded_samples(hash_excluded_samples,output_path):
	
	output_filename = os.path.join(output_path,'xhmm_excluded_samples.tsv')
	
	fo = open(output_filename,"w")
	
	for sample in sorted(hash_excluded_samples.keys()):
		fo.write("%s\t%s\n" % (sample,";".join(sorted(hash_excluded_samples[sample]))))
	
	fo.close()
	
	return output_filename

##########################################################################

def run_mosdepth(bam_file,bed_analysis,mosdepth_path,output_path,Q_threshold=20,**kwargs):
	
	correct_overlap = kwargs.get('overlap',False)
	Q_threshold     = kwargs.get('base_qual',15)
		
	fileName,extName = os.path.splitext(os.path.basename(bam_file))
	fileName = os.path.join(output_path,fileName+"_Mosdepthcoverage")
		
	if extName <> ".cram":
		args = [mosdepth_path,fileName,'--no-per-base','-t','12','-b',bed_analysis,'-F','4','-T','0,20,50','-Q',"%d" % (Q_threshold),bam_file]
	else:
		ref_fasta = kwargs.get('ref_fasta',None)
		
		args = [mosdepth_path,fileName,'--no-per-base','-t','12','-b',bed_analysis,'-f',ref_fasta,'-F','4','-T','0,20,50','-Q',"%d" % (Q_threshold),bam_file]
			
	if correct_overlap:
		args = args + ['--fast-mode']
		
	try:
		mosdepth_sal = Popen(args,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(output,logdata) = mosdepth_sal.communicate()
		mosdepth_sal.wait()
	except:
		raise RuntimeError("run_coverage: Error in mosdepth with bam file: %s" % (bam_file))
	
	if logdata.lower().find('error') <> -1:
		raise RuntimeError("run_coverage: Error in mosdepth with bam file: %s\n%s" % (bam_file,logdata))
			
	covFileName_i = fileName+'.regions.bed.gz'
	covFileName_o = fileName+".regions.bed"
		
	f = gzip.open(covFileName_i, 'rb')
	file_content = f.read()
	f.close()
			
	fo = open(covFileName_o,'w')
	fo.write(file_content)
	fo.close()
		
	for f in glob.glob(fileName+'.*gz.csi'):
		os.remove(f) 
	for f in glob.glob(fileName+'*thresholds*'):
		os.remove(f)
	for f in glob.glob(fileName+'.*mosdepth*txt'):
		os.remove(f)
	for f in glob.glob(fileName+'.*regions.bed.gz'):
		os.remove(f)
	
	return covFileName_o

##########################################################################

def __sort_key(x):
	
	chrom = x[0]

	if chrom[3:].isdigit():
		return int(chrom[3:]),int(x[1]),int(x[2])
	else:
		return x[0],int(x[1]),int(x[2])
	
def __sort_key2(x):
	
	chrom = x[0]

	if chrom[3:].isdigit():
		return int(chrom[3:]),int(x[1]),int(x[2]),x[3]
	else:
		return x[0],int(x[1]),int(x[2]),x[3]
	
def __sort_key3(x):
	
	chrom = x[0]

	if chrom[3:].isdigit():
		return int(chrom[3:]),int(x[1]),int(x[2]),x[3:]
	else:
		return x[0],int(x[1]),int(x[2]),x[3:]
	

##########################################################################

def __parse_cov_file(covFileName):
	
	bed_cov = BedTool(covFileName)

	hash_cov_parsed = dict(map(lambda intv: ((intv.chrom,intv.end,intv.start),float(intv.score)), bed_cov))
	
	return hash_cov_parsed

##########################################################################

def merge_coverage_into_matrix(hash_cov):
	
	l_samples = hash_cov.keys() 
	
	hash_cov_parsed = {}
	
	for s in l_samples:
		
		hash_cov_parsed_s = __parse_cov_file(hash_cov[s])
		
		map(lambda intv: hash_cov_parsed.setdefault(intv,[]).append(hash_cov_parsed_s[intv]), hash_cov_parsed_s.keys())
	
	l_intv = sorted(hash_cov_parsed.keys(),key=__sort_key)
	
	cov_mat = np.zeros((len(l_samples),len(l_intv)),float)
	
	for (i,intv) in enumerate(l_intv):
		
		l_cov = hash_cov_parsed[intv]

		cov_mat[:,i] = l_cov
		
	l_intv_columns = map(lambda (chrm,s,e): "%s:%d-%d" % (chrm,min(s,e),max(s,e)),l_intv)
		
	cov_mat_df = pd.DataFrame(cov_mat,columns=l_intv_columns,index=l_samples)
	
	return cov_mat_df 

##########################################################################

def xhmm_center_and_filter(cov_mat_filename,hash_parameters,xhmm_path,label,output_path):
	
	"""
		### 1. Filters samples and targets and then mean-centers the targets
		
		### Input: 
		###	coverage_mat and configuration parameters
		
		### Output:
		### coverage_mat_centered: Centered coverage matrix
		### excluded_targets and excluded_samples: Files with the excluded intervals and samples due to not supply parameters above explained
	
		### Parameters:
		
		--minTargetSize:    Minimum size of target (in bp) to process (10)
		--maxTargetSize:    Maximum size of target (in bp) to process (1000)
		--minMeanTargetRD:  Minimum per-target mean RD to require for target to be processed (10)
		--maxMeanTargetRD:  Maximum per-target mean RD to require for target to be processed (2000)
		--minMeanSampleRD:  Minimum per-sample mean RD to require for sample to be processed (10)
		--maxMeanSampleRD:  Maximum per-sample mean RD to require for sample to be processed (2000)
		--minSdSampleRD:    Minimum per-sample standard deviation of RD to require for sample to be processed (0)
		--maxSdSampleRD:    Maximum per-sample standard deviation of RD to require for sample to be processed (500)
            
        ### Call:                      
		$xhmm_path --matrix --centerData --centerType target -r $coverage_mat -o $coverage_mat_centered  \
		--outputExcludedTargets $excluded_targets \
		--outputExcludedSamples $excluded_samples \
		--minTargetSize 10 --maxTargetSize 10000 \
		--minMeanTargetRD 10 --maxMeanTargetRD 2000 \
		--minMeanSampleRD 10 --maxMeanSampleRD 2000 \
		--minSdSampleRD 0 --maxSdSampleRD 500 
	"""

	cov_mat_center_filename = os.path.join(output_path,"xhmm.filtered_centered.%s.mat" % (label))
	
	filename_excluded_intv    = os.path.join(output_path,'xhmm.excluded_intervals_by_center.%s.tsv' % (label))
	filename_excluded_samples = os.path.join(output_path,'xhmm.excluded_samples_by_center.%s.tsv' % (label))
	
	args_param = ['--minTargetSize',str(hash_parameters['minTargetSize']),'--maxTargetSize',str(hash_parameters['maxTargetSize']),'--minMeanTargetRD',str(hash_parameters['minMeanTargetRD']),'--maxMeanTargetRD',str(hash_parameters['maxMeanTargetRD']),'--minMeanSampleRD',str(hash_parameters['minMeanSampleRD']),'--maxMeanSampleRD',str(hash_parameters['maxMeanSampleRD']),'--minSdSampleRD',str(hash_parameters['minSdSampleRD']),'--maxSdSampleRD',str(hash_parameters['maxSdSampleRD'])]
	args = [xhmm_path,'--matrix','--centerData','--centerType','target','-r',cov_mat_filename,'-o',cov_mat_center_filename,'--outputExcludedTargets',filename_excluded_intv,'--outputExcludedSamples',filename_excluded_samples]
	
	try:
		xhmm_sal = Popen(args+args_param,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(output,logdata) = xhmm_sal.communicate()
		xhmm_sal.wait()
	except:
		err = sys.exc_info()[0]
		raise RuntimeError("run_xhmm.xhmm_filter_and_center: Error in xhmm center data: %s" % (err))
	
	if logdata.lower().find('error')<>-1:
		raise RuntimeError("run_xhmm.xhmm_filter_and_center: Error in xhmm center data: %s" % (logdata))
	
	return (cov_mat_center_filename,filename_excluded_intv,filename_excluded_samples)

##########################################################################

def xhmm_PCA_normalizate_and_filter(cov_mat_center_filename,xhmm_path,label,output_path):
	
	"""
		The function runs PCA with centered and filtered coverage-matrix. 
		Then, it normalizes mean-centered data using PCA information
		Finally, it filters and z-score centers (by sample) the PCA-normalized data
	"""
	
	### 2.a Runs PCA on mean-centered data:
	
	### Call:
	# $xhmm_path --PCA -r $coverage_mat_centered --PCAfiles $pca_file
	
	### Input:
	### cov_mat_center_filename
	
	### Output:
	### pca_filenames: *PC_LOADINGS.txt, *PC_SD.txt, *PC.txt (this last one is PCA data)
	
	pca_filename               = os.path.join(output_path,'xhmm.%s.pca' % (label))
	pca_filename_norm          = os.path.join(output_path,'xhmm.%s.pca.norm.txt' % (label))
	pca_filename_norm_filtered = os.path.join(output_path,'xhmm.%s.pca.norm.filter.txt' % (label))
	
	excluded_intervals = os.path.join(output_path,'xhmm.excluded_intervals_by_zscore.%s.tsv' % (label))
	excluded_samples   = os.path.join(output_path,'xhmm.excluded_samples_by_zscore.%s.tsv' % (label))
		
	args1 = [xhmm_path,'--PCA','-r',cov_mat_center_filename,'--PCAfiles',pca_filename]
	
	try:
		xhmm_sal = Popen(args1,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(output,logdata1) = xhmm_sal.communicate()
		xhmm_sal.wait()
	except:
		err = sys.exc_info()[0]
		raise RuntimeError("run_xhmm.xhmm_filter_and_center: Error in xhmm center data: %s" % (err))
	
	if logdata1.lower().find('error')<>-1:
		raise RuntimeError("run_xhmm.xhmm_PCA_normalizate_filter: Error in xhmm: %s" % (logdata1))
	
	### 2.b Normalizes mean-centered data using PCA information:
	# --PCnormalizeMethod:  Method to choose which prinicipal components are removed for data normalization  (possible values="numPCtoRemove", "PVE_mean", "PVE_contrib" default=`PVE_mean')
	# --numPCtoRemove:      Number of highest principal components to filter out  (default=`20')
	# --PVE_mean_factor:    Remove all principal components that individually explain more variance than this factor times the average (in the original PCA-ed data)  (default=`0.7')
	# --PVE_contrib:        Remove the smallest number of principal components that explain this percent of the variance (in the original PCA-ed data) (default=`50')
	# $xhmm_path --normalize -r $coverage_mat_centered --PCAfiles $pca_file --normalizeOutput $pca_file_norm --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
	
	args2 = [xhmm_path,'--normalize','-r',cov_mat_center_filename,'--PCAfiles',pca_filename,'--normalizeOutput',pca_filename_norm,'--PCnormalizeMethod','PVE_mean','--PVE_mean_factor','0.7']
	
	try:
		xhmm_sal = Popen(args2,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(output,logdata2) = xhmm_sal.communicate()
		xhmm_sal.wait()
	except:
		err = sys.exc_info()[0]
		raise RuntimeError("run_xhmm.xhmm_filter_and_center: Error in xhmm center data: %s" % (err))
	
	if logdata2.lower().find('error')<>-1:
		raise RuntimeError("run_xhmm.xhmm_PCA_normalizate_filter: Error in xhmm: %s" % (logdata2))
		
	### 2.c. Filters and z-score centers (by sample) the PCA-normalized data:
	# $xhmm_path --matrix --centerData --centerType sample --zScoreData -r $pca_file_norm  -o $pca_file_norm_filtered \
	# --outputExcludedTargets $excluded_targets_filter \
	# --outputExcludedSamples $excluded_samples_filter --maxSdTargetRD 30 
	
	args3 = [xhmm_path,'--matrix','--centerData','--centerType','sample','--zScoreData','-r',pca_filename_norm,'-o',pca_filename_norm_filtered,'--outputExcludedTargets',excluded_intervals,'--outputExcludedSamples',excluded_samples,'--maxSdTargetRD','30']

	try:
		xhmm_sal = Popen(args3,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(output,logdata3) = xhmm_sal.communicate()
		xhmm_sal.wait()
	except:
		err = sys.exc_info()[0]
		raise RuntimeError("run_xhmm.xhmm_filter_and_center: Error in xhmm center data: %s" % (err))
	
	if logdata3.lower().find('error')<>-1:
		raise RuntimeError("run_xhmm.xhmm_PCA_normalizate_filter: Error in xhmm: %s" % (logdata3))
	
	return pca_filename_norm_filtered,excluded_intervals,excluded_samples

##########################################################################

def xhmm_create_final_matrix(cov_mat_filename,excluded_intervals1,excluded_samples1,excluded_intervals2,excluded_samples2,xhmm_path,output_path):
	
	"""
		### Filters original read-depth data to be the same as filtered, normalized data:
		
		$xhmm_path --matrix -r $coverage_mat \
		--excludeTargets $excluded_targets_filter --excludeTargets $excluded_targets \
		--excludeSamples $excluded_samples --excludeSamples $excluded_samples_filter \
		-o $coverage_mat_filtered
	"""
	
	cov_mat_filename_filtered = os.path.splitext(cov_mat_filename)[0]+'.final.mat'
	
	args = [xhmm_path,'--matrix','-r',cov_mat_filename,'--excludeTargets',excluded_intervals1,'--excludeTargets',excluded_intervals2,'--excludeSamples',excluded_samples1,'--excludeSamples',excluded_samples2,'-o',cov_mat_filename_filtered]
	
	try:
		xhmm_sal = Popen(args,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(output,logdata) = xhmm_sal.communicate()
		xhmm_sal.wait()
	except:
		err = sys.exc_info()[0]
		raise RuntimeError("run_xhmm.xhmm_filter_and_center: Error in xhmm center data: %s" % (err))
	
	if logdata.lower().find('error')<>-1:
		raise RuntimeError("run_xhmm.xhmm_filter_and_center: Error in xhmm center data: %s" % (logdata))
	
	return cov_mat_filename_filtered

##########################################################################

def xhmm_discover_cnvs(cov_mat_filtered,pca_data,label,xhmm_path,output_path):
	
	"""
		### It discovers CNVs in normalized data
	 
		### Call:
			$xhmm_path --discover -p $xhmm_params_file \
			-r $pca_file_norm_filtered -R $coverage_mat_filtered \
			-c $cnv_output -a $cnv_output_aux -s $posterior_path
		
		### Input:
			cov_mat_filtered: Coverage matrix with intervals and samples removed
			pca_data: PCA components
			param_filename: file with parameters. Copy from default distribution
			output_path_posterior: path of directory with posterior prob
			
		### Output:
			cnv_output_filename: This file is the final result
			cnv_output_filename_aux
	"""
	
	param_filename = os.path.join(output_path,"params.txt")
	f = open(param_filename,"w")
	f.write("1e-8	6	70	-3	1.00	0	1.00	3	1.00\n")
	f.close()
	
	output_path_posterior = os.path.join(output_path,'posterior')
	
	cnv_output_filename     = os.path.join(output_path,'xhmm.%s.xcnv' % (label))
	cnv_output_filename_aux = os.path.join(output_path,'xhmm.%s.aux.xcnv' % (label))
	
	args = [xhmm_path,'--discover','-p',param_filename,'-R',cov_mat_filtered,'-r',pca_data,'-s',output_path_posterior,'-c',cnv_output_filename,'-a',cnv_output_filename_aux]
	
	try:
		xhmm_sal = Popen(args,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(output,logdata) = xhmm_sal.communicate()
		xhmm_sal.wait()
	except:
		err = sys.exc_info()[0]
		raise RuntimeError("run_xhmm.xhmm_filter_and_center: Error in xhmm center data: %s" % (err))
	
	if logdata.lower().find('error')<>-1:
		raise RuntimeError("run_xhmm.xhmm_filter_and_center: Error in xhmm center data: %s" % (logdata))
	
	return cnv_output_filename

##########################################################################

def split_into_samples_and_save(l_bams,cnv_output_filename,bed_analysis,output_path):

	"""
		SAMPLE: sample name
		CNV: type of copy number variation (DEL or DUP)
		INTERVAL: genomic range of the called CNV
		KB: length in kilobases of called CNV
		CHR: chromosome name on which CNV falls
		MID_BP: the midpoint of the CNV (to have one genomic number for plotting a single point, if desired)
		TARGETS: the range of the target indices over which the CNV is called (NOTE: considering only the FINAL set of post-filtering targets)
		NUM_TARG: # of exome targets of the CNV
		Q_EXACT: Phred-scaled quality of the exact CNV event along the entire interval- Identical to EQ in .vcf output from genotyping
		Q_SOME: Phred-scaled quality of some CNV event in the interval - Identical to SQ in .vcf output from genotyping
		Q_NON_DIPLOID:  Phred-scaled quality of not being diploid, i.e., DEL or DUP event in the interval - Identical to NDQ in .vcf output from genotyping
		Q_START: Phred-scaled quality of "left" breakpoint of CNV - Identical to LQ in .vcf output from genotyping
		Q_STOP: Phred-scaled quality of "right" breakpoint of CNV - Identical to RQ in .vcf output from genotyping
		MEAN_RD: Mean normalized read depth (z-score) over interval - Identical to RD in .vcf output from genotyping
		MEAN_ORIG_RD: Mean read depth (# of reads) over interval - Identical to ORD in .vcf output from genotyping
	"""
	
	fo = open(cnv_output_filename,'r')
	l_lines = map(lambda x: x.strip().split('\t'), fo.readlines())
	fo.close()
	
	hash_cnvs = {}

	for line in l_lines[1:]:
				
		chrom	 = line[4]
		interval = line[2]
		(start,end) = interval.replace(chrom+":",'').split('-')
		
		sample = line[0]
		
		cnv_type   = line[1]
		phred_qual = line[8]  #Q_EXACT
		dosis	   = line[13] #MEAN_RD
		
		#hash_cnvs.setdefault(sample,[]).append("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom,start,end,sample,cnv_type,dosis,phred_qual))
		hash_cnvs.setdefault(sample,[]).append((chrom,start,end,sample,cnv_type,dosis,phred_qual))
	
	### Annotation with reference bed 
	l_bed_cnv = []
		
	bed_file = BedTool(bed_analysis)
	
	for bam_file in l_bams:
		
		sample = os.path.splitext(os.path.basename(bam_file))[0]
		
		cnv_bed_filename = os.path.join(output_path,"%s_xhmm_cnv.bed" % (sample))
		
		if hash_cnvs.has_key(sample):
			
			l_intv = []
			
			for intv in BedTool(sorted(hash_cnvs[sample],key=__sort_key3)).intersect(bed_file,wb=True):
				
				chrom  = intv[0]
				start  = intv[1]
				end    = intv[2]
				sample = intv[3]
				cnv_type = intv[4]
				dosis  = intv[5]
				phred  = intv[6]
			
				if len(intv) == 11:
					gene  = intv[10]
				else:
					gene = '.'
					
				if len(intv) == 13:
					nm    = intv[11]
					exon  = intv[12]
				else:
					nm    = '.'
					exon  = '.'
				
				l_intv.append((chrom,start,end,sample,cnv_type,dosis,phred,gene,nm,exon))

			BedTool(l_intv).saveas(cnv_bed_filename)

		else:
			open(cnv_bed_filename,'w').close()
			
		l_bed_cnv.append(cnv_bed_filename)
	
	return l_bed_cnv

##########################################################################

def __parse_coord(coord):
		
	try:
		(chom,_,start,_,end) = re.split(r'(-|:|)?',coord)
	except:
		return None
	
	return (chom,start,end)

##########################################################################

def split_matrix_in_chromosomes(mat_cov_filtered_filename,hash_gender,tmp_path):
	
	root_name = os.path.splitext(os.path.basename(mat_cov_filtered_filename))[0]
	
	mat_cov = pd.read_csv(mat_cov_filtered_filename,sep="\t",header=0,index_col=0)
	
	l_coord = map(lambda coord: (coord,__parse_coord(coord)[0]), list(mat_cov.columns))
	
	l_coord_chrX    = map(itemgetter(0),filter(lambda (c,chrom): chrom.lower().find("x")<>-1, l_coord))
	l_coord_chrY    = map(itemgetter(0),filter(lambda (c,chrom): chrom.lower().find("y")<>-1 , l_coord))
	l_coord_no_sex  = map(itemgetter(0),filter(lambda (c,chrom): chrom.lower().find("x")==-1 and chrom.lower().find("y")==-1 , l_coord))
	
	l_samples_XY = filter(lambda s: hash_gender[s]=="XY", hash_gender.keys())
	l_samples_XX = filter(lambda s: hash_gender[s]=="XX", hash_gender.keys()) 
			
	mat_cov_autosomes = mat_cov.loc[:,l_coord_no_sex]
	mat_cov_autosomes_filename = os.path.join(tmp_path,root_name+'.autosomes.mat')
	mat_cov_autosomes.to_csv(mat_cov_autosomes_filename,sep='\t')
	
	mat_cov_chrX_f = mat_cov.loc[l_samples_XX,l_coord_chrX]
	mat_cov_chrX_f_filename = os.path.join(tmp_path,root_name+'.chrX.F.mat')
	mat_cov_chrX_f.to_csv(mat_cov_chrX_f_filename,sep='\t')
	
	mat_cov_chrX_m = mat_cov.loc[l_samples_XY,l_coord_chrX]
	mat_cov_chrX_m_filename = os.path.join(tmp_path,root_name+'.chrX.M.mat')
	mat_cov_chrX_m.to_csv(mat_cov_chrX_m_filename,sep='\t')
	
	mat_cov_chrY = mat_cov.loc[l_samples_XY,l_coord_chrY]
	mat_cov_chrY_filename = os.path.join(tmp_path,root_name+'.chrY.M.mat')
	mat_cov_chrY.to_csv(mat_cov_chrY_filename,sep='\t')
	
	return mat_cov_autosomes_filename,mat_cov_chrX_f_filename,mat_cov_chrX_m_filename,mat_cov_chrY_filename

##########################################################################

def join_cnvs(cnv_output_filename_autosomes,cnv_output_filename_chrX_F,cnv_output_filename_chrX_M,cnv_output_filename_chrY_M,output_path):
	
	cnv_joined_filename = os.path.join(output_path,'xhmm.joined.xcnv')
	
	f = open(cnv_joined_filename,'w')
	
	if cnv_output_filename_autosomes <> None:
		f.write(open(cnv_output_filename_autosomes,'r').read())
	else:
		f.write("SAMPLE\tCNV\tINTERVAL\tKB\tCHR\tMID_BP\tTARGETS\tNUM_TARG\tQ_EXACT\tQ_SOME\tQ_NON_DIPLOID\tQ_START\tQ_STOP\tMEAN_RD\tMEAN_ORIG_RD\n")
	if cnv_output_filename_chrX_F <> None:
		f.write("".join(open(cnv_output_filename_chrX_F,'r').readlines()[1:]))
	if cnv_output_filename_chrX_M <> None:
		f.write("".join(open(cnv_output_filename_chrX_M,'r').readlines()[1:]))
	if cnv_output_filename_chrY_M <> None:
		f.write("".join(open(cnv_output_filename_chrY_M,'r').readlines()[1:]))
	
	f.close()
	
	return cnv_joined_filename

##########################################################################

def annotate_excluded_intervals(l_excluded_intervals_filename,bed_analysis,output_path):
	
	l_excluded_intervals = []
	
	for filename in l_excluded_intervals_filename:
		
		if filename == None:
			continue
		
		f = open(filename,'r')
		l_excluded_intervals.extend(map(lambda x: x.strip().split('\t')+[os.path.splitext(os.path.basename(filename))[0]], f.readlines()))
		f.close()
	
	hash_excluded_coord = dict(map(lambda (coord,step): (__parse_coord(coord),step) , l_excluded_intervals))
	
	l_tagged_intv = []
	
	for intv in BedTool(bed_analysis):
		
		key_intv = (intv[0],intv[1],intv[2])
		
		num_feat = len(intv.fields)
		
		gene = nm = exon = '.' 
		
		if num_feat == 6:
			gene  = intv[3]
			nm    = intv[4]
			exon  = intv[5]
		elif num_feat == 4:
			gene = intv[3]
		elif num_feat == 5:
			gene = intv[3]
			nm = intv[4]
		
		if hash_excluded_coord.has_key(key_intv):
			l_tagged_intv.append((key_intv+(gene,nm,exon) + ('excluded',hash_excluded_coord[key_intv])))
		else:
			l_tagged_intv.append((key_intv+(gene,nm,exon) + ('.','.')))
	
	filename_tagged_bed = os.path.join(output_path,'xhmm.tagged_excluded_intervals.bed')
	
	BedTool(l_tagged_intv).saveas(filename_tagged_bed)
	
	return filename_tagged_bed
	
##########################################################################

def get_excluded_samples(l_bamfiles,excluded_samples1,excluded_samples2,label):
	
	l_samples_excluded1 = []
	l_samples_excluded2 = []
	
	if excluded_samples1 == None:
		l_samples_excluded1 = map(lambda s: (os.path.splitext(os.path.basename(s))[0],label), l_bamfiles)
	elif os.stat(excluded_samples1).st_size == 0:
		l_samples_excluded1 = []
	else:
		l_samples_excluded1 = map(lambda x: (x.strip(),label), open(excluded_samples1,'r').readlines())
	
	if excluded_samples2 == None:
		l_samples_excluded2 = map(lambda s: (os.path.splitext(os.path.basename(s))[0],label), l_bamfiles)
	elif os.stat(excluded_samples2).st_size == 0:
		l_samples_excluded2 = []
	else:
		l_samples_excluded2 = map(lambda x: (x.strip(),label), open(excluded_samples2,'r').readlines())
		
	return list(set(l_samples_excluded1+l_samples_excluded2))

##########################################################################

def write_xlsx(output_xlsx_file,param_config_file,excluded_samples_filename,bed_tagged_file,cnv_file,**kwargs):

	writer = pd.ExcelWriter(output_xlsx_file, engine='xlsxwriter')
	
	df_param = pd.read_csv(param_config_file, delimiter='\t', na_values=['.'], header=None, index_col=None,names=['Field','value'])
	df_param.to_excel(writer, index=False, sheet_name='XHMM param')
	
	df_excl = pd.read_csv(excluded_samples_filename, delimiter='\t', na_values=['.'], header=None, index_col=None,names=['sample','set excluded'])
	df_excl.to_excel(writer, index=False, sheet_name='samples excluded')
	
	df_intv = pd.read_csv(bed_tagged_file, delimiter='\t', na_values=['.'], header=None, index_col=None, names = ['Chrom','start','end','gene','exon','transcript','is_excluded','filter'])
	df_intv.to_excel(writer, index=False, sheet_name='Intervals analyzed')
	
	if kwargs.get('is_excluded',False):
		df_bed = pd.DataFrame(columns=['sample excluded by the XHMM tool for cnv analysis'])
		df_bed.to_excel(writer, index=False, sheet_name='XHMM cnvs')
	else:
		if os.stat(cnv_file).st_size == 0:
			df_bed = pd.DataFrame(columns=['There are not CNVs determined in this sample'])
		else:
			df_bed = pd.read_csv(cnv_file, delimiter='\t', na_values=['.'], header=None, index_col=None,names=['Chrom','start','end','sample','cnv type','dosis','phred score qual','gene','exon','transcript'])
				
		df_bed.to_excel(writer, index=False, sheet_name='XHMM cnvs')
	
	writer.save()
	writer.close()

##########################################################################

def get_chromosome_sex(sample,bam_file,bed_file_chrom,samtools_path,output_path,ref_fasta=None):

	fi = open(bed_file_chrom,'r')
	l_regions = map(lambda r: r.strip().split('\t'), fi.readlines())
	fi.close()

	chrY_regions = []
	chrX_regions = []
	
	l_pos_chr = groupby(map(lambda (i,r): (i,r[0]), enumerate(l_regions)),lambda x: x[1])
	
	for chrom,l_pos in l_pos_chr:
		if chrom == "chrX" or chrom == "X":
			chrX_regions = map(lambda (pos,c): l_regions[pos], l_pos)
		elif chrom == "chrY" or chrom == "Y":
			chrY_regions = map(lambda (pos,c): l_regions[pos], l_pos)
	
	mean_depth_chrX = 0
	
	sex_output_file = os.path.join(output_path,sample+'.gender')
	
	fo_g = open(sex_output_file,'w')
	fo_g.write("Interval\tTotal reads\tMean depth\n")
	
	for region in chrX_regions:
		
		chrm  = region[0]
		start = int(region[1])
		end   = int(region[2])

		if ref_fasta <> None:
			args_samtools = [samtools_path,'depth',bam_file,'--reference',ref_fasta,'-r',"%s:%d-%d" % (chrm,start,end)]
		else:
			args_samtools = [samtools_path,'depth',bam_file,'-r',"%s:%d-%d" % (chrm,start,end)]
		
		samtools_sal = Popen(args_samtools,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(str_stdout,logdata) = samtools_sal.communicate()
		samtools_sal.wait()
		
		if logdata <> "":
			if logdata.lower().find('fail')<>-1 or logdata.lower().find('error')<>-1:
				raise IOError('CNV_analysis.get_chromosome_sex: error in samtools depth: %s' % (logdata))
		
		total_depth = mean_depth_intv = 0
		
		if str_stdout <> "":
			
			total_depth = sum(map(lambda (chrm,pos,depth): int(depth), map(lambda l: l.split('\t'),str_stdout.strip().split('\n'))))
			
			mean_depth_intv = float(total_depth)/(end-start+1)
		
			mean_depth_chrX += mean_depth_intv
			
		fo_g.write("%s:%d-%d\t%d\t%1.2f\n" % (chrm,start,end,total_depth,mean_depth_intv))
	
	mean_depth_chrX = mean_depth_chrX/len(chrX_regions)
	
	mean_depth_chrY = 0
	
	for region in chrY_regions:
		
		chrm  = region[0]
		start = int(region[1])
		end   = int(region[2])
		
		if ref_fasta <> None:
			args_samtools = [samtools_path,'depth',bam_file,'--reference',ref_fasta,'-r',"%s:%d-%d" % (chrm,start,end)]
		else:
			args_samtools = [samtools_path,'depth',bam_file,'-r',"%s:%d-%d" % (chrm,start,end)]
	
		samtools_sal = Popen(args_samtools,stdin=PIPE, stdout=PIPE, stderr=PIPE,close_fds=True,bufsize=1)
		(str_stdout,logdata) = samtools_sal.communicate()
		samtools_sal.wait()
		
		total_depth = mean_depth_intv = 0
		
		if str_stdout <> "":
			
			total_depth = sum(map(lambda (chrm,pos,depth): int(depth), map(lambda l: l.split('\t'),str_stdout.strip().split('\n'))))
			
			mean_depth_intv = float(total_depth)/(end-start+1)
		
			mean_depth_chrY += mean_depth_intv
			
		fo_g.write("%s:%d-%d\t%d\t%1.2f\n" % (chrm,start,end,total_depth,mean_depth_intv))
			
	fo_g.close()
		
	mean_depth_chrY = mean_depth_chrY/len(chrY_regions)
	
	#print mean_depth_chrX,mean_depth_chrY
	
	if mean_depth_chrY >= 6 and mean_depth_chrX > 0:
		return 'XY'
	elif mean_depth_chrX > 2:
		return 'XX'
	else:
		return "UN"
	
#######################################################################

def run(argv=None):

	if argv is None: argv = sys.argv    

	parser = OptionParser(add_help_option=True,description="The script launches XHMM tool that determines cnvs")
		
	parser.add_option("--i",default=None,help="Input file with the list of bam files to calculate CNVs with xhmm",dest="f_listbam")
	parser.add_option("--f",default=None,action="append",help="Bam file",dest="f_bam")
	parser.add_option("--b",default=None,help="Bed file with the region of interest",dest="bed_file")
	parser.add_option("--o",default=None,help="Path of output xhmm files",dest="output_path")
	parser.add_option("--cfg",default=None,help="Configuration file cfg",dest="cfg_file")
	parser.add_option("--no-cov",action="store_true",default=False,help="Flag to disconnect coverage calculus due to the files are already generated",dest="no_cov")
	parser.add_option("--clean",action="store_true",default=False,help="Flag to remove the temporal directory",dest="clean")
	
	# Se leen las opciones aportadas por el usuario
	(options, args) = parser.parse_args(argv[1:])

	if len(argv) == 1:
		sys.exit(0)

	if not parser.check_required("--b"):
		raise IOError('run_xhmm: The bed file with intervals has not been provided')
	if not parser.check_required("--o"):
		raise IOError('run_xhmm: The output path has not been provided')
	if not parser.check_required("--cfg"):
		raise IOError('run_xhmm: The cfg file has not been provided')
	
	l_bams_raw = []
	
	if options.f_listbam <> None:
		
		if not os.path.exists(options.f_listbam):
			raise IOError('run_xhmm: The file with the list of bam files does not exist. %s' % (options.f_listbam))
		
		fi = open(options.f_listbam,'r')
		l_bams_raw = map(lambda l: l.strip(), fi.readlines())
		fi.close()
	
	elif options.f_bam <> None:
		l_bams_raw = options.f_bam 
	else:
		raise IOError("run_xhmm: No bam file has been provided")
	
	l_bamfiles = []
		
	for bam_file in l_bams_raw:
			
		if bam_file[0] == '#':
			continue
			
		if not os.path.exists(bam_file) and not options.no_cov:
			raise IOError("run_xhmm: The bam file does not exist:\n%s" % (bam_file))
			
		l_bamfiles.append(bam_file)
			
		
	bed_file = options.bed_file
	
	if not os.path.exists(bed_file):
		raise IOError('run_xhmm: The bed file does not exist: %s' % (bed_file))
	
	output_path = options.output_path
	
	xhmm_output_path = os.path.join(output_path,"xhmm")
	
	if not os.path.exists(output_path):
		raise IOError('run_xhmm: The output file does not exist: %s' % (output_path))
	
	if not os.path.exists(xhmm_output_path):
		raise IOError('run_xhmm: The xhmm output file does not exist: %s' % (xhmm_output_path))
			
	cov_output_path = os.path.join(output_path,"coverage")
			
	if not os.path.exists(cov_output_path):
		os.mkdir(cov_output_path)
		
	tmp_path = os.path.join(xhmm_output_path,'tmp')
	
	if not os.path.exists(tmp_path):
		os.mkdir(tmp_path)
		
	cfg_filename = options.cfg_file
	
	if not os.path.exists(cfg_filename):
		raise IOError('run_xhmm: The cfg file does not exist: %s' % (cfg_filename))
	
	hash_cfg = read_cfg_file(cfg_filename)
	
	xhmm_path = hash_cfg.get("xhmm_path",'None')
	
	if not os.path.exists(xhmm_path):
		raise IOError('run_xhmm: The path of xhmm tool does not exist: %s' % (xhmm_path))
	
	samtools_path = hash_cfg.get("samtools_path",'None')
	
	if not os.path.exists(samtools_path):
		raise IOError('run_xhmm: The path of samtools tool does not exist: %s' % (samtools_path))
	
	sex_markers_bed = hash_cfg.get("chrom_sexual",'None')
	
	if not os.path.exists(sex_markers_bed):
		raise IOError('run_xhmm: The bed of sexual markers does not exist: %s' % (sex_markers_bed))
		
	hash_param = configure_process(xhmm_path)
	
	
	#Configure logger
	formatter = logging.Formatter('%(asctime)s - %(module)s - %(levelname)s - %(message)s')
	console = logging.StreamHandler()
	console.setFormatter(formatter)
	console.setLevel(logging.INFO)

	log_file = os.path.join(xhmm_output_path,"xhmm.log")

	if os.path.exists(log_file):
		os.remove(log_file)

	logging.basicConfig(filename=log_file,level=logging.INFO,format='%(asctime)s,%(msecs)d %(levelname)-8s %(message)s',datefmt='%Y-%m-%d:%H:%M:%S')

	logger = logging.getLogger("logging")
	
	if ( len( logger.handlers ) == 0 ):
		logger.addHandler(console)	
	
	logger.info("Checking the chromosomal sex of input samples:")
		
	### 0. It is checked the chromosomal sex of the samples
	hash_gender = OrderedDict()
		
	for bam_file in l_bamfiles:
		
		sample = os.path.splitext(os.path.basename(bam_file))[0]
		
		chrom_sex = get_chromosome_sex(sample,bam_file,sex_markers_bed,samtools_path,tmp_path,ref_fasta=None)
				
		hash_gender[sample] = chrom_sex
		
		logger.info("%s --> %s" % (sample,chrom_sex))
		
	logger.info("done!!\n") 
		
	### 1a. It is run the coverage and it is created a final joined matrix
	hash_cov = OrderedDict()
	
	if options.no_cov:
		
		logger.info("WARNING: It was disabled the coverage option\n")
		
		for bam_file in l_bamfiles:
		
			fileName,extName = os.path.splitext(os.path.basename(bam_file))
			covFileName_o = os.path.join(cov_output_path,fileName+"_Mosdepthcoverage.regions.bed")
		
			if not os.path.exists(covFileName_o):
				raise IOError("run_xhmm: It has been selected no_cov option but coverage files does not exist.\n%s" % (covFileName_o))
			
			hash_cov[fileName] = covFileName_o 
	else:
		
		Q_THRES = hash_param['Q_threshold']
		
		mosdepth_path = hash_cfg.get("mosdepth_path",'None')
		
		if not os.path.exists(mosdepth_path):
			raise IOError('run_xhmm: The path of mosdepth tool does not exist: %s' % (mosdepth_path))
		
		logger.info("Coverage calculus:")
		
		for bam_file in l_bamfiles:
			
			logger.info("Performing coverage of file: %s" % (bam_file))
			
			fileName,extName = os.path.splitext(os.path.basename(bam_file))
			covFileName_o = run_mosdepth(bam_file,bed_file,mosdepth_path,cov_output_path,Q_THRES) ## If cram file, check that fasta file is provided
			hash_cov[fileName] = covFileName_o
					
			
		logger.info("done!!\n")
	
	### The joined matrix is then, splitted into sexual/non-sexual chromosmes
	logger.info("Building coverage matrix")
	
	mat_cov_filename = os.path.join(xhmm_output_path,"xhmm.mat")
	
	cov_mat_df = merge_coverage_into_matrix(hash_cov)
	
	cov_mat_df.to_csv(mat_cov_filename,sep='\t')
	
	mat_cov_autosomes_filename,mat_cov_chrX_f_filename,mat_cov_chrX_m_filename,mat_cov_chrY_filename = split_matrix_in_chromosomes(mat_cov_filename,hash_gender,tmp_path)
	
	logger.info("done!!\n")
		
	### Autosome calling
	logger.info("CNV calling of autosome regions")
	
	try:
		label = "autosomes"
		
		### 1. Filters samples and targets and then mean-centers the targets
		(mat_cov_autosomes_filename_center,excluded_intervals_autosomes1,excluded_samples_autosomes1) = xhmm_center_and_filter(mat_cov_autosomes_filename,hash_param,xhmm_path,label,tmp_path)
		### 2. PCA step
		(pca_norm_autosomes_filtered,excluded_intervals_autosomes2,excluded_samples_autosomes2) = xhmm_PCA_normalizate_and_filter(mat_cov_autosomes_filename_center,xhmm_path,label,tmp_path)
		### 3. Final filtering
		mat_cov_autosomes_filtered_filename = xhmm_create_final_matrix(mat_cov_autosomes_filename,excluded_intervals_autosomes1,excluded_samples_autosomes1,excluded_intervals_autosomes2,excluded_samples_autosomes2,xhmm_path,tmp_path)
		### 4. CNV calling
		cnv_output_filename_autosomes = xhmm_discover_cnvs(mat_cov_autosomes_filtered_filename,pca_norm_autosomes_filtered,label,xhmm_path,tmp_path)
	except Exception as e: 
		raise RuntimeError("%s" % (e))
		
		#excluded_intervals_autosomes1 = excluded_intervals_autosomes2 = None
		#excluded_samples_autosomes1   = excluded_samples_autosomes2 = None
		#cnv_output_filename_autosomes = None
	
	### ChrX in females calling 
	logger.info("CNV calling of chrX regions in female individuals")
	try:
		label = "chrX_F"
		
		### 1. Filters samples and targets and then mean-centers the targets
		(mat_cov_chrX_f_filename_center,excluded_intervals_chrX_f1,excluded_samples_chrX_f1) = xhmm_center_and_filter(mat_cov_chrX_f_filename,hash_param,xhmm_path,label,tmp_path)
		### 2. PCA step
		(pca_norm_chrX_f_filtered,excluded_intervals_chrX_f2,excluded_samples_chrX_f2) = xhmm_PCA_normalizate_and_filter(mat_cov_chrX_f_filename_center,xhmm_path,label,tmp_path)
		### 3. Final filtering
		mat_cov_chrX_f_filtered_filename = xhmm_create_final_matrix(mat_cov_chrX_f_filename,excluded_intervals_chrX_f1,excluded_samples_chrX_f1,excluded_intervals_chrX_f2,excluded_samples_chrX_f2,xhmm_path,tmp_path)
		### 4. CNV calling
		cnv_output_filename_chrX_F = xhmm_discover_cnvs(mat_cov_chrX_f_filtered_filename,pca_norm_chrX_f_filtered,label,xhmm_path,tmp_path)
	except:
		excluded_intervals_chrX_f1 = excluded_intervals_chrX_f2 = None 
		excluded_samples_chrX_f1   = excluded_samples_chrX_f2 = None
		cnv_output_filename_chrX_F = None
	
	### ChrX in males calling
	logger.info("CNV calling of chrX regions in male individuals")
	try:
		label = "chrX_M"
		
		### 1. Filters samples and targets and then mean-centers the targets
		(mat_cov_chrX_m_filename_center,excluded_intervals_chrX_m1,excluded_samples_chrX_m1) = xhmm_center_and_filter(mat_cov_chrX_m_filename,hash_param,xhmm_path,label,tmp_path)
		### 2. PCA step
		(pca_norm_chrX_m_filtered,excluded_intervals_chrX_m2,excluded_samples_chrX_m2) = xhmm_PCA_normalizate_and_filter(mat_cov_chrX_m_filename_center,xhmm_path,label,tmp_path)
		### 3. Final filtering
		mat_cov_chrX_m_filtered_filename = xhmm_create_final_matrix(mat_cov_chrX_m_filename,excluded_intervals_chrX_m1,excluded_samples_chrX_m1,excluded_intervals_chrX_m2,excluded_samples_chrX_m2,xhmm_path,tmp_path)
		### 4. CNV calling
		cnv_output_filename_chrX_M = xhmm_discover_cnvs(mat_cov_chrX_m_filtered_filename,pca_norm_chrX_m_filtered,label,xhmm_path,tmp_path)
	except:
		excluded_intervals_chrX_m1 = excluded_intervals_chrX_m2 = None 
		excluded_samples_chrX_m1   = excluded_samples_chrX_m2 = None
		cnv_output_filename_chrX_M = None
	
	### ChrY in males calling
	logger.info("CNV calling of chrY regions in male individuals")
	try:
		label = "chrY_M"
		
		### 1. Filters samples and targets and then mean-centers the targets
		(mat_cov_chrY_filename_center,excluded_intervals_chrY1,excluded_samples_chrY1) = xhmm_center_and_filter(mat_cov_chrY_filename,hash_param,xhmm_path,label,tmp_path)
		### 2. PCA step
		(pca_norm_chrY_filtered,excluded_intervals_chrY2,excluded_samples_chrY2) = xhmm_PCA_normalizate_and_filter(mat_cov_chrY_filename_center,xhmm_path,label,tmp_path)
		### 3. Final filtering
		mat_cov_chrY_filtered_filename = xhmm_create_final_matrix(mat_cov_chrY_filename,excluded_intervals_chrY1,excluded_samples_chrY1,excluded_intervals_chrY2,excluded_samples_chrY2,xhmm_path,tmp_path)
		### 4. CNV calling
		cnv_output_filename_chrY_M = xhmm_discover_cnvs(mat_cov_chrY_filtered_filename,pca_norm_chrY_filtered,label,xhmm_path,tmp_path)
	except:
		excluded_intervals_chrY1   = excluded_intervals_chrY2 = None 
		excluded_samples_chrY1     = excluded_samples_chrY2 = None
		cnv_output_filename_chrY_M = None
	
	cnv_output_filename = join_cnvs(cnv_output_filename_autosomes,cnv_output_filename_chrX_F,cnv_output_filename_chrX_M,cnv_output_filename_chrY_M,xhmm_output_path)
	
	logger.info("done!!\n")
			
	### 5. Saving results
	### 5.1 Get excluded samples into a list
	l_excluded_samples_autosomes = get_excluded_samples(l_bamfiles,excluded_samples_autosomes1,excluded_samples_autosomes2,"autosomes")
	l_excluded_samples_chrX_f    = get_excluded_samples(l_bamfiles,excluded_samples_chrX_f1,excluded_samples_chrX_f2,"chrX_F")
	l_excluded_samples_chrX_m    = get_excluded_samples(l_bamfiles,excluded_samples_chrX_m1,excluded_samples_chrX_m2,"chrX_M")
	l_excluded_samples_chrY      = get_excluded_samples(l_bamfiles,excluded_samples_chrY1,excluded_samples_chrY2,"chrY")
	
	hash_excluded_samples = {} 
	map(lambda (s,l): hash_excluded_samples.setdefault(s,[]).append(l), sorted(set(l_excluded_samples_autosomes+l_excluded_samples_chrX_f+l_excluded_samples_chrX_m+l_excluded_samples_chrY)))
	excl_filename = write_excluded_samples(hash_excluded_samples,xhmm_output_path)
		
	### 5.2 Save bed files with cnvs detected
	l_bed_cnvs = split_into_samples_and_save(l_bamfiles,cnv_output_filename,bed_file,xhmm_output_path)
		
	### 5.3 Save a bed with excluded intervals tagged  
	bed_analysis_tagged = annotate_excluded_intervals([excluded_intervals_autosomes1,excluded_intervals_autosomes2,excluded_intervals_chrX_f1,excluded_intervals_chrX_f2,excluded_intervals_chrX_m1,excluded_intervals_chrX_m2,excluded_intervals_chrY1,excluded_intervals_chrY2],bed_file,xhmm_output_path)
	
	### 5.4 Write configuration param file
	config_filename = write_configuration(hash_param,xhmm_output_path)

	### 5.5 Write the xlsx files
	for i,bam_file in enumerate(l_bamfiles):
		
		sample = os.path.splitext(os.path.basename(bam_file))[0]
		
		logger.info("Writing the results of sample: %s" % (sample))
						
		is_excluded = False
		
		if len(hash_excluded_samples.get(sample,[])) == 4:
			is_excluded = True
		
		output_xlsx_file = os.path.join(xhmm_output_path,"%s_xhmm_cnv.xlsx" % (sample))
		
		write_xlsx(output_xlsx_file,config_filename,excl_filename,bed_analysis_tagged,l_bed_cnvs[i],is_excluded=is_excluded)
		
	if options.clean:
		shutil.rmtree(tmp_path)
	
	logger.info("XHMM finished!!")
	
###########################################################################

if __name__=='__main__':
	
	try:
		run()
	except Exception as e: 
		print(e)
