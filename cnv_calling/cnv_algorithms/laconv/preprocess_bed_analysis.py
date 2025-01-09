'''
Created on July 2020

@author: adelpozo
'''

#!/usr/bin/python

import sys, os

import optparse

from pybedtools import BedTool

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

#########################################################################

def __sort_by_chrom(x):
	
	chrom = x[0]

	if chrom[3:].isdigit():
		return int(chrom[3:]),int(x[1]),int(x[2])
	else:
		return x[0],int(x[1]),int(x[2])

#########################################################################

def run(argv=None):

	if argv is None: argv = sys.argv    

	parser = OptionParser(add_help_option=True,description="The script performs CNV estimation within the regions of interest")
	
	parser.add_option("--b",default=None,help="Bed file with the intervals to preprocess",dest="f_bed")
	parser.add_option("--o",default=None,help="Output path where to write the preprocessed bed",dest="output_path")
	parser.add_option("--no-flank",action="store_true",default=False,help="Flag to disconnect the inclusion of 10bp as flanking region",dest="no_flank")
	parser.add_option("--s",default=10,help="lenght of the slop to increase each interval. Default: 10 bp",dest="slop")
	
	# Se leen las opciones aportadas por el usuario
	(options, args) = parser.parse_args(argv[1:])

	if len(argv) == 1:
		sys.exit(0)

	if not parser.check_required("--b"):
		raise IOError('The bed file has not been provided')

	if not parser.check_required("--o"):
		raise IOError('The output path has not been provided')

	bed_analysis = options.f_bed
	
	if not os.path.exists(bed_analysis):
		raise IOError('preprocess_bed_analysis: The bed file provided does not exist: %s' % (bed_analysis))
	
	output_path = options.output_path
	
	if not os.path.exists(output_path):
		raise IOError('preprocess_bed_analysis: The output path must be created before running the script')
	
	slop_length = int(options.slop)

	try:
		
		fi = open(bed_analysis,'r')
		l_intervals_raw = map(lambda l: l.strip().split('\t'), fi.readlines())
		fi.close()
	
		l_intervals = []
	
		for i,l_intv in enumerate(l_intervals_raw):
		
			chrom = l_intv[0]
			try:
				start = int(l_intv[1])
				end   = int(l_intv[2])
			except:
				raise RuntimeError('preprocess_bed_analysis: Line %i of bed file is malformed')
		
			gene = '.'
			exon = '.'
			nm   = '.'
		
			if len(l_intv) >= 4:
				gene = l_intv[3]
			if len(l_intv) >= 5:
				exon = l_intv[4]
				if exon.find('exon')==-1:
					exon = '.'
			if len(l_intv) >= 7:
				nm = l_intv[6]
				if nm.find("NM_") == -1 and nm.find("NR_") == -1:
					nm = '.'
			
			start_new = min(start,end)
			end_new   = max(start,end)
			
			if not options.no_flank:
				if start_new == 0:
					start_new = start
				else:
					start_new = start_new-slop_length
				end_new = end_new+slop_length
			
			if start_new < 0 or end_new < 0:
				print "The line %i of the bed is malformed" % (i)
				continue
			
			l_intervals.append((chrom,start_new,end_new,gene,exon,nm))
			
		if l_intervals == []:
			raise RuntimeError('preprocess_bed_analysis: The bed file is empty or there are too many malformed lines')
	
		l_intervals_sorted = sorted(l_intervals,key=__sort_by_chrom)
		
		preproc_bed_analysis = os.path.join(output_path,os.path.splitext(os.path.basename(bed_analysis))[0]+'.preprocessed.bed')
		
		BedTool(l_intervals_sorted).merge(c='4,5,6',o='distinct',delim='|').saveas(preproc_bed_analysis)
		
		sys.stdout.write('Created the bed file: %s\n' % (preproc_bed_analysis))
		
	except:
		print >> sys.stderr , '\n%s\t%s' % (sys.exc_info()[0],sys.exc_info()[1])
		sys.exit(2)

############################################################################

if __name__=='__main__':
	
	run()
