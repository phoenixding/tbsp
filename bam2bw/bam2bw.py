#!/usr/bin/env python
"""
bam2bw.py <bam_file> <chromesome size file>

<bam_file>: bam file
<chromesome size file>: chromsome size files downloaded form UCSC genome browser (e.g. mm10.chrom.szies) 

function:
convert bam files to bigWig files

prerequisite:
samtools
bedtools
bedGraphtoBigWig 
"""

import pdb,sys,os,subprocess
from subprocess import Popen,PIPE
import argparse
	
def requirements_check():
    """
    Ensure we have programs needed to download/manipulate the data
    """
    required_programs = [
        ('samtools',
         'http://samtools.sourceforge.net/'),
        ('bedtools',
         'http://bedtools.readthedocs.org/en/latest/'),
        ('bedGraphToBigWig',
         'http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/'),
    ]
    for req, url in required_programs:
        try:
            p = subprocess.Popen(
                [req], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stdout, stderr = p.communicate()
        except OSError:
            raise ValueError("Please install %s (%s)" % (req, url))



def main():
	requirements_check()
	print("all needed programs were installed! Start processing bam files ...")
	if len(sys.argv[1:])!=2:
		print(__doc__)
		sys.exit(0)
		
	bamFile_pre=sys.argv[1]
	genomeSize=sys.argv[2]
	# sort bam files
	subpx=subprocess.Popen(["samtools","sort",bamFile_pre,bamFile_pre+'.sort'],stdout=PIPE,stderr=PIPE)
	stdout,stderr=subpx.communicate()
	#pdb.set_trace()
	if (len(stderr)>0):
		subpx=subprocess.Popen(["samtools","sort","-o",bamFile_pre+".sort.bam", bamFile_pre],stdout=PIPE,stderr=PIPE)
		stdout,stderr=subpx.communicate()
	
	bamFile=bamFile_pre+".sort.bam"
	outfile=open(bamFile+'.bedGraph','w')
	subp=subprocess.Popen(["bedtools","genomecov","-ibam",bamFile,"-bga","-g",genomeSize],stdout=outfile,stderr=PIPE)
	stdout,stderr=subp.communicate()
	outfile.close()

	# sort bedGraph
	sortoutfile=open(bamFile+'.bedGraph.sort','w')
	subp2=subprocess.Popen(['sort',"-k1,1","-k2,2n",bamFile+'.bedGraph'],stdout=sortoutfile,stderr=PIPE)
	stdout,stderr=subp2.communicate()
	sortoutfile.close()
	
	# bedGraph to bigwig
	subp3=subprocess.Popen(["bedGraphToBigWig",bamFile+".bedGraph.sort",genomeSize,bamFile+'.bedGraph.sort.bw'],stdout=PIPE,stderr=PIPE)
	stdout,stderr=subp3.communicate()
	
	# remove temp files
	
	try:
		os.remove(bamFile)
		os.remove(bamFile+'.bedGraph')
		os.remove(bamFile+'.bedGraph.sort')
	
	except:
		pass 
	
if __name__=="__main__":
	main()
	

