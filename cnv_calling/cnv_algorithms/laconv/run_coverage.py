'''
Created on August 2020

@author: adelpozo
'''

#!/usr/bin/python


import sys, re, shlex , os, string, urllib, time, math, random, subprocess, shutil

import ConfigParser, optparse

from subprocess import Popen , PIPE

import gzip

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

def run_mosdepth(bam_list,output_path,bed_analysis,mosdepth_path,ref_fasta,**kwargs):

	correct_overlap = kwargs.get('overlap',False)
	Q_threshold     = kwargs.get('base_qual',15)
		
	for bam_file in bam_list:
		
		fileName,extName = os.path.splitext(os.path.basename(bam_file))
		fileName = os.path.join(output_path,fileName+"_Mosdepthcoverage")
		
		sys.stdout.write("%s\n" % (bam_file))
		sys.stdout.flush()
		
		if extName <> ".cram":
			args = [mosdepth_path,fileName,'--no-per-base','-t','12','-b',bed_analysis,'-F','4','-T','0,20,50','-Q',"%d" % (Q_threshold),bam_file]
		else:
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
			
	return

##########################################################################

def run(argv=None):

	if argv is None: argv = sys.argv    

	parser = OptionParser(add_help_option=True,description="The script runs mosdepth tool to get coverage files for CNV analysis")
		
	parser.add_option("--i",default=None,help="File with the list of bam files",dest="f_listbam")
	parser.add_option("--f",default=None,action="append",help="Bam file",dest="f_bam")
	parser.add_option("--b",default=None,help="Bed file with the list of intervals",dest="bed_filename")
	parser.add_option("--r",default=None,help="Fasta reference file (only for .cram files)",dest="fasta_filename")
	parser.add_option("--m",default=None,help="mosdepth path",dest="mosdepth_path")
	parser.add_option("--o",default=None,help="Output path of coverage directory",dest="output_path")
	parser.add_option("--no-overlap",action="store_true",default=False,help="Flag to correct mate overlaps in coverage counting",dest="no_overlap")
		
	# Se leen las opciones aportadas por el usuario
	(options, args) = parser.parse_args(argv[1:])

	if len(argv) == 1:
		sys.exit(0)

	if not parser.check_required("--b"):
		raise IOError('run_coverage: The bed file has not been provided')
	if not parser.check_required("--m"):
		raise IOError('run_coverage: The path of mosdepth tool has not been provided')
	if not parser.check_required("--o"):
		raise IOError('run_coverage: The output path has not been provided')
	
	bed_filename = options.bed_filename
	
	if not os.path.exists(bed_filename):
		raise IOError("run_coverage: The bed does not exist %s" % (bed_filename))
		
	mosdepth_path = options.mosdepth_path
	
	if not os.path.exists(mosdepth_path):
		raise IOError("run_coverage: The path of mosdepth tool does not exist %s" % (mosdepth_path))
			
	output_path = options.output_path
	
	if not os.path.exists(output_path):
		raise IOError("run_coverage: The output path provided does not exist %s" % (output_path))
		
	if options.f_listbam <> None:
		
		if not os.path.exists(options.f_listbam):
			raise IOError('run_coverage: The file with the list of bam files does not exist')
		
		fi = open(options.f_listbam,'r')
		l_bams_raw = map(lambda l: l.strip(), fi.readlines())
		fi.close()
	
	elif options.f_bam <> None:
		l_bams_raw = options.f_bam 
	else:
		raise IOError("run_coverage: No bam file has been provided")
			
	l_bams = []
	
	is_cram = False
		
	for bam_file in l_bams_raw:
			
		if bam_file[0] == '#':
			continue
			
		if not os.path.exists(bam_file) and not options.no_cov:
			raise IOError("run_coverage: The bam file does not exist:\n%s" % (bam_file))
		
		if os.path.splitext(bam_file)[1] == '.cram':
			is_cram = True
			
		l_bams.append(bam_file)
	
	ref_fasta = options.fasta_filename
	
	if is_cram:
		if not os.path.exists(ref_fasta):
			raise IOError("run_coverage: The fasta file does not exist:\n%s" % (ref_fasta))
	
	sys.stdout.write("Running mosdepth...\n\noutput path: %s\n\n" % (output_path))
	sys.stdout.flush()
	
	run_mosdepth(l_bams,output_path,bed_filename,mosdepth_path,ref_fasta,overlap=options.no_overlap)
		
	sys.stdout.write("\nFinished!!!")
	sys.stdout.flush()
	

############################################################################

if __name__=='__main__':
	
	run()
