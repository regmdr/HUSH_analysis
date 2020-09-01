import sys
import pybedtools
import pysam
import pyBigWig
import numpy as np

##############################################



def genome_file_as_dict(fileName):
    """ 
    Function to read genome file. 
    It will also provide the list of chromosomes in the correct order.
    """
    chrom_sizes = {}
    chrom_names = []
    with open(fileName, "r") as f:
        for line in f:
            values = line.strip().split("\t")
            if len(values) != 2:
                raise ValueError("Wrong format for genome file")
            else:
                chrom_sizes[values[0]] = values[1]
                chrom_names.append(values[0])
    return chrom_sizes, chrom_names


def get_chrom_windows_dict(chromosomeDict, windowSize):
    bedtoolDict = {}
    for chrom in chromosomeDict:
        if int(chromosomeDict[chrom]) < windowSize:
            raise ValueError("Error, chromosome " 
                             + chrom + " smaller than window")
        chromString = chrom + '\t1\t' + chromosomeDict[chrom]
        chrom_bt = pybedtools.BedTool(chromString, from_string=True)
        bedtoolDict[chrom] = pybedtools.BedTool().window_maker(b=chrom_bt, 
                                                               w=windowSize)
    return bedtoolDict


def get_covArray(bamAlignFile, chrWindowBed):
    """ 
    Function to get read count for a chromosome, 
    using pysam. Using an indexed bam file 
    """
    chromArray = np.array([])
    for w in chrWindowBed:
        chromArray = np.append(chromArray, 
                               bamAlignFile.count(w.chrom, w.start, w.end))
        # with this method, PE reads are counted as singles
    return chromArray


def feat_count(featureBed, chrWindowBed):
    """ 
    Function to get feature counts, 
    using pybedtools. Use bed file 
    """
    # could use number of features or percent of the window covered?
    chromArray = np.array([])
    repChrom = chrWindowBed.intersect(featureBed, c=True)
    for w in repChrom:
        chromArray = np.append(chromArray, w.count)
    return chromArray


def bigW_count(bigWigSignal, chrWindowBed):
    """ 
    Function to get bigwig signal, 
    using pybigwig. Use bigwig
    """
    # could use number of features or percent of the window covered?
    chromArray = np.array([])
    for w in chrWindowBed:
        regionStats=bigWigSignal.stats(w.chrom, w.start, w.end, type="mean")[0]
        if regionStats is None:
            chromArray = np.append(chromArray, 0)
        else:
            chromArray = np.append(chromArray, regionStats)
    return chromArray


def get_values_All_tracks(files, formatList, chromList, chromWindict, windowSize):
    numDataSets = len(files)
    chrNum = len(chromList)
    data_allChrom = []
    for dataCount in range(numDataSets):
        dataTracks = []
        if formatList[dataCount] == "bed":
            features = pybedtools.BedTool(files[dataCount])
            if len(features) == 0:
                raise RuntimeError("Please provide at least one region for plotting!")
            for chrom in chromList:
                chromWindow = chromWindict[chrom]
                fArray = feat_count(features, chromWindow)
                dataTracks.append(fArray)
        elif formatList[dataCount] == "bam":
            bam = pysam.AlignmentFile(files[dataCount], "rb")
            totMapReads = bam.mapped
            for chrom in chromList:
                chromWindow = chromWindict[chrom]
                tempArray = get_covArray(bam, chromWindow)
                rpkm = (10e9*tempArray) / (totMapReads*windowSize)
                dataTracks.append(rpkm)
            bam.close()
        elif formatList[dataCount] == "bigWig":
            bigWigO = pyBigWig.open(files[dataCount])
            for chrom in chromList:
                chromWindow = chromWindict[chrom]
                bwArray = bigW_count(bigWigO, chromWindow)
                dataTracks.append(bwArray)
            bigWigO.close()
        else:
            raise ValueError("Wrong format")
        if dataCount == 0:
            data_allChrom = dataTracks
        else:
            for cNum in range(chrNum):
                data_allChrom[cNum] = np.vstack((data_allChrom[cNum], 
                                                dataTracks[cNum]))
    return data_allChrom




##############################################

# General arguments:      
genome_file = "hg38.chromosomes.genome"
windowSize = 500
percentile_t = 85

# files = ["siren_high.merged.fwd.bw", "siren_high.merged.rev.bw"]
fw_file = sys.argv[1]
rev_file = sys.argv[2]
outprfx = sys.argv[3]
files = [fw_file, rev_file]
fileF = ["bigWig", "bigWig"]

# Window name
wname = windowSize

# Getting ordered list of chromosome and dict with lengths
genomeDict, chrom_names = genome_file_as_dict(genome_file)
  
# Dictionary of bedtools with window-divided chromosomes
chromWdict = get_chrom_windows_dict(genomeDict, windowSize)
data_allChrom= get_values_All_tracks(files, fileF, chrom_names, chromWdict, windowSize)

# Flattening data:
array1 = np.array([])
array2 = np.array([])
for chrom in data_allChrom:
    array1 = np.append(array1, chrom[0])
    array2 = np.append(array2, chrom[1])

# Combining data excluding bins with expression of 0.
non_zero_array = all_array[all_array != 0]
# Getting quantile threshold
exp_thr2 = np.percentile(non_zero_array,percentile_t)


# Which bins have expression above the threshold in both strands?
chrom_list = []
position_list = []
c = 0 
for chrom in data_allChrom:
    for i in range(chrom.shape[1]):
        if chrom[0][i] > exp_thr2 and chrom[1][i] > exp_thr2:
            chrom_list.append(chrom_names[c])
            start = (i*windowSize) + 1
            end = (i*windowSize) + 1 + windowSize
            position_list.append([start, end])
    c += 1


# Writing a bed file with bins expressed from both strands:
out_file_name = outprfx + "_expressed_both_strands" + str(percentile_t) + "_non-zero." + wname + ".bed"
oFILE = open(out_file_name, "w")
for i in range(len(chrom_list)):
    oFILE.write( "%s\t%s\t%s\n" % (chrom_list[i], position_list[i][0], position_list[i][1]))
oFILE.close()
