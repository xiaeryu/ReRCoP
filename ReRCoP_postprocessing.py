import copy
import ReRCoP_preprocessing


def mergeInterval(intervals):
    '''
    This function merges intercals in the sliding window approach.

    @param intervals: list of list [[start.1, end.1], [start.2, end.2]]
    @return list of list after merging
    '''
    if len(intervals) == 0:
        return []
    output = [intervals[0]]
    for pair in intervals[1:]:
        if pair[0] <= output[-1][1]:
            output[-1][1] = max(pair[1], output[-1][1])
        else:
            output.append(pair)
    return output


def removeOutlier(inFasta, inMat):
    '''
    This function removes the outlier genes from the aligned fasta file generated by sliding window.

    @param inFasta: a fasta object of the input genome
    @param inMat: the matrix (Outliermat) recording the outlier genes
    @return a fasta object with outlier genes removed.
    '''
    fasta = copy.deepcopy(inFasta)

    for i in range(3,len(inMat[0])):
        oInterval = []
        for j in range(1,len(inMat)):
            if inMat[j][i] == 1:
                oInterval.append([inMat[j][1], inMat[j][2]])
        merged = mergeInterval(oInterval)

        for pair in merged:
            pair = map(int, pair)
            fasta[inMat[0][i]] = fasta[inMat[0][i]][:pair[0]-1] + '-'*(pair[1]-pair[0]+1) + fasta[inMat[0][i]][pair[1]:]

    return fasta


def writeFasta(inFasta, outfile):
    '''
    This function writes a fasta object into a fasta file

    @param inFasta: input fasta object
    @param outfile: output fasta file
    '''
    outH = open(outfile, 'w')
    for key in inFasta:
        outH.write(">%s\n" % key)
        outH.write("%s\n" % inFasta[key])
    outH.close()


def writeMat(inMat, outfile):
    '''
    This function writes a matrix to the output file

    @param inMat: input matrix
    @param outfile: output file
    '''
    outH = open(outfile, 'w')
    for i in range(len(inMat)):
        outH.write(str(inMat[i][0]))
        for j in range(1, len(inMat[i])):
             outH.write("\t%s" % str(inMat[i][j]))
        outH.write('\n')
    outH.close()