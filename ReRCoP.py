import ReRCoP_checkPrerequisite 
import ReRCoP_preprocessing
import ReRCoP_matrix
import ReRCoP_outlierDetection
import ReRCoP_postprocessing
import copy
import sys
import os
from optparse import OptionParser
from optparse import OptionGroup


#################################################################
##### Input options
#################################################################
# Parse input options
usage = "usage: %prog [options] Genomes.fasta"
parser = OptionParser(usage=usage, version="%prog 1.0")

group = OptionGroup(parser, "Input Options")
group.add_option("-a", "--aligned", action="store_true", dest="aligned", help="Set this if genome sequences in the input file are aligned.")
group.add_option("--gbk",action="store",type="string",dest="gbk",help="Input GenBank file of the reference genome.")
group.add_option("-w", "--window", action="store_true", dest="window", help="Set this if sliding windows instead of genes are to be considered.")
group.add_option("--fSize",action="store",type="int",dest="fSize",default=1000,help="Fragment size if using sliding window. [Default: 1000]")
group.add_option("--sSize",action="store",type="int",dest="sSize",default=500,help="Step size if using sliding window. [Default: 500]")
group.add_option("--cds",action="store",type="string",dest="cds",help="Input coding sequences in fasta format. Used to determine the core genome when input genomes are not aligned.")
parser.add_option_group(group)

group = OptionGroup(parser, "Core Gene Identification Options")
group.add_option("--cov",action="store",type="float",dest="cov",default=0.7,help="Minimum sequence coverage to regard genes as present [Default: 0.7]")
group.add_option("--sim",action="store",type="int",dest="sim",default=70,help="Minimum sequence similarity to regard genes as present [Default: 70]")
parser.add_option_group(group)

group = OptionGroup(parser, "Outlier Removal Options")
group.add_option("-m", "--method",action="store",type="string",dest="method",help="Outlier removal method. Can be 'Grubbs', 'kNN', or 'DBSCAN', or can be multiple methods separated by ','")
group.add_option("--alpha", action="store",type="float",dest="alpha",default=0.05,help="For 'Grubbs' method: Significance level in Grubbs test. [Default: 0.05]")
group.add_option("--radius", action="store",type="float",dest="radius",default=1.5,help="For 'kNN' method: Maximum number of differences for a point to be considered as a neighbor (in the unit of standard deviation of all pair-wise nubmer of differences). [Default: 1.5)]")
group.add_option("--k",action="store",type="float",dest="k", default=0.2, help="For 'kNN' method: Minimum number of neighbors for a non-outlier point (in the unit of total number of points). [Default: 0.2]")
group.add_option("--eps",action="store",type="float",dest="eps",default=1, help="For 'DBSCAN' method: Maximum number of differences between two points for them to be considered as in the same neighborhood (in the unit of standard deviation of all pair-wise nubmer of differences). [Default: 1]")
group.add_option("--minP",action="store",type="float",dest="minP",default=0.2, help="For 'DBSCAN' method: Minimum number of points required to form a dense region (in the unit of total number of points). [Default: 0.2]")
parser.add_option_group(group)

group = OptionGroup(parser, "Output Options")
group.add_option("-o", "--outdir",action="store",type="string",dest="outdir",default=".", help="Output directory. [Default: running directory]")
group.add_option("-p", "--prefix",action="store",type="string",dest="prefix",default="ReRCoP", help="Output prefix. [Default: ReRCoP]")
parser.add_option_group(group)

(options, args) = parser.parse_args()


# Check input options
if len(args) == 0:
    print usage
    sys.exit()

if options.aligned:
    if options.cds:
        parser.error("Options --aligned and --cds are mutually exclusive.")
    if options.gbk and options.window:
        parser.error("Options --gbk and --window are mutually exclusive.")
    if not options.gbk and not options.window:
        parser.error("Either --gbk or --window should be set while using --aligned.")
    if options.window:
        if options.sSize <=0 or options.sSize>options.fSize:
            parser.error("Option --sSize should be larger than 0 and no larger than --fSize")
else:
    if not options.cds:
        parser.error("Option --cds should be indicated while not --aligned.")

if options.cov<0 or options.cov>1:
    parser.error("Option --cov should be within the range of 0-1.")
if options.sim<0 or options.sim>100:
    parser.error("Option --sim should be within the range of 0-100")

if not options.method:
    parser.error("Option --method is required.")

if options.alpha<0 or options.alpha>1:
    parser.error("Option --alpha should be within the range of 0-1.")
if options.radius <0:
    parser.error("Option --radius should be no less than 0.")
if options.k<0 or options.k>1:
    parser.error("Option --k should be within the range of 0-1.")
if options.eps <0:
    parser.error("Option --eps should be no less than 0.")
if options.minP <0 or options.minP >1:
    parser.error("Option --minP should be within the range of 0-1.")


# Get the values of the options
inputGenome = args[0]		# Input complete genome sequences

aligned = options.aligned	# Input genome sequences are aligned or not
inputGbk = options.gbk		# Input genebank file
window = options.window		# Using sliding windows
fragSize = int(options.fSize)	# Fragment size of sliding windows
stepSize = int(options.sSize)	# Step size of sliding windows
inputGene = options.cds		# input coding sequences

covCut = float(options.cov)	# 0-1, coverage cutoff in core genome identification
simCut = float(options.sim)	# 0-100, similarity cutoff in core genome identification

outlierMethod = options.method.split(',')	# Outlier removal method
for item in outlierMethod:
    if item not in ['Grubbs', 'kNN', 'DBSCAN']:
        parser.error("Option --method should be 'Grubbs', 'kNN', or 'DBSCAN'")

alpha = float(options.alpha)	# Grubbs: Significance level
radius = float(options.radius)	# kNN: parameter
k = options.k			# kNN: parameter
eps = options.eps		# DBSCAN: radius
minP = options.minP		# minP: minPts

outdir = options.outdir		# output directory
prefix = options.prefix		# output prefix


#################################################################
##### Check prerequisite
#################################################################

# Check python module
ReRCoP_checkPrerequisite.checkModule('scipy')

# Check the necessary software
if inputGene:
    ReRCoP_checkPrerequisite.checkCommand('which makeblastdb')
    ReRCoP_checkPrerequisite.checkCommand('which blastn')

# Check input files and duplicate names in the fasta files
ReRCoP_checkPrerequisite.checkFile(inputGenome)
ReRCoP_checkPrerequisite.checkName(inputGenome)
if inputGbk:
    ReRCoP_checkPrerequisite.checkFile(inputGbk)
if inputGene:
    ReRCoP_checkPrerequisite.checkFile(inputGene)
    ReRCoP_checkPrerequisite.checkName(inputGene)

# Check the output directory
ReRCoP_checkPrerequisite.checkDir(outdir)


#################################################################
##### Form concatenate core genome and generate concatenation log
#################################################################

seqConcat = {}	# A fasta object of the concatenated genomes
logConcat = []	# A list of list as the concatenation log
inGenome = ReRCoP_preprocessing.readFasta(inputGenome)
if aligned:
    ReRCoP_checkPrerequisite.checkLen(inGenome)

# Input: complete genome + coding sequences, require identification
if not aligned and inputGene:
    inGene = ReRCoP_preprocessing.readFasta(inputGene)
    tmpFile = outdir + '/' + prefix + ".ReRCoP.tmp"
    [seqConcat, logConcat] = ReRCoP_preprocessing.parseRaw(inGene, inGenome, simCut, covCut, tmpFile+".1", tmpFile+".2", tmpFile+".3")

# Input: sequence alignment + gbk file, parse based on gbk
if aligned and inputGbk:
    inGbk = ReRCoP_preprocessing.readGbk(inputGbk)
    [seqConcat, logConcat] = ReRCoP_preprocessing.parseGbk(inGbk, inGenome, covCut)

# Input: sequence alignmet, keep non-coding regions in a sliding-window manner.
if aligned and window:
    fullLen = ReRCoP_preprocessing.fastaLen(inGenome)
    logConcat = ReRCoP_preprocessing.slidingWindow(fullLen, fragSize, stepSize)
    seqConcat = inGenome


################################################################
##### Generate matrix of SNP number for each gene
################################################################
SNPmat = ReRCoP_matrix.GeneDiff(logConcat, seqConcat)

################################################################
##### Outlier detection
################################################################
if "Grubbs" in outlierMethod:
    Outliermat_Grubbs = ReRCoP_matrix.methodMat(SNPmat)
    for i in range(1,len(SNPmat)):
        SNPs = SNPmat[i][3:]
        outlierPos = ReRCoP_outlierDetection.Grubbs(SNPs, alpha)
        for j in outlierPos:
            Outliermat_Grubbs[i][j+3] = 1
    finalConcat_Grubbs = ReRCoP_postprocessing.removeOutlier(seqConcat, Outliermat_Grubbs)
    ReRCoP_postprocessing.writeMat(Outliermat_Grubbs, outdir+'/'+prefix+".Grubbs.outliermat")
    ReRCoP_postprocessing.writeFasta(finalConcat_Grubbs, outdir+'/'+prefix+".Grubbs.removal.fasta")


if "kNN" in outlierMethod:
    Outliermat_kNN = ReRCoP_matrix.methodMat(SNPmat)
    for i in range(1,len(SNPmat)):
        SNPs = SNPmat[i][3:]
        outlierPos = ReRCoP_outlierDetection.kNN(SNPs, k, radius)
        for j in outlierPos:
            Outliermat_kNN[i][j+3] = 1
    finalConcat_kNN = ReRCoP_postprocessing.removeOutlier(seqConcat, Outliermat_kNN)
    ReRCoP_postprocessing.writeMat(Outliermat_kNN, outdir+'/'+prefix+".kNN.outliermat")
    ReRCoP_postprocessing.writeFasta(finalConcat_kNN, outdir+'/'+prefix+".kNN.removal.fasta")

if 'DBSCAN' in outlierMethod:
    Outliermat_DBSCAN = ReRCoP_matrix.methodMat(SNPmat)
    if not minP:
            minP = len(SNPs)*0.3
    for i in range(1,len(SNPmat)):
        SNPs = SNPmat[i][3:]
        outlierPos = ReRCoP_outlierDetection.DBSCAN(SNPs, eps, minP)
        for j in outlierPos:
            Outliermat_DBSCAN[i][j+3] = 1

    finalConcat_DBSCAN = ReRCoP_postprocessing.removeOutlier(seqConcat, Outliermat_DBSCAN)
    ReRCoP_postprocessing.writeMat(Outliermat_DBSCAN, outdir+'/'+prefix+".DBSCAN.outliermat")
    ReRCoP_postprocessing.writeFasta(finalConcat_DBSCAN, outdir+'/'+prefix+".DBSCAN.removal.fasta")


################################################################
##### Write to output
################################################################

## write the SNP mat
ReRCoP_postprocessing.writeMat(SNPmat, outdir+'/'+prefix+".snpmat")

## write the Concatinated sequences
ReRCoP_postprocessing.writeFasta(seqConcat, outdir+'/'+prefix+".core.fasta")

## write the log file
ReRCoP_postprocessing.writeMat(logConcat, outdir+'/'+prefix+".concatenation.log")
