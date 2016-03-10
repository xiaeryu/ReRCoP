import math
import copy
import ReRCoP_preprocessing


def mutCount(ref, query, log):
    '''
    This function compares the reference sequence and the query sequence
    and returns the number of differences in each gene as defined in the
    log file.

    @param ref: a string of the reference sequence
    @param query: a string of the query sequence
    @param log: the concatenation log
    @return [total, storage] with total the total nubmer of differences
    and storage a list with the number of differences in each region 
    defined in the log
    '''
    storage = []
    total = 0
    for i in range(len(log)):
        diff = 0
        for j in range(int(log[i][1])-1,int(log[i][2])):
            if query[j] != ref[j]:
                diff += 1
                total += 1
        storage.append(diff)
    return [total, storage] # [int, list]


def median(inarr):
    '''
    This function calculates the median of a list.

    @param inarr: a numerical list
    @return the median
    '''
    tmp = inarr[:]
    sorts = sorted(tmp)
    l = len(sorts)
    if l % 2 == 0:
        return (int(sorts[l/2]) + int(sorts[l/2-1])) / 2
    return sorts[l / 2]


def average(inarr):
    '''
    This function calculates the mean of a list.
    
    @parm inarr: a numerical list
    @return the mean
    '''
    return sum(inarr) * 1.0 / len(inarr)


def stdDeviation(inarr):
    '''
    This function calculates the standard deviation of a list.

    @param inarr: a numerical list
    @return the standard deviation
    '''
    avg = average(inarr)
    variance = map(lambda x: (x - avg)**2, inarr)
    return math.sqrt(average(variance))


def GeneDiff(record, inFasta):
    '''
    This function calculates the relative number of SNPs

    @param record: the log file returned by concatenation
    @param inFasta: fasta object of the sequence concatenations
    @return list of list of relative number of SNPs
    '''
    ref = ReRCoP_preprocessing.consensus(inFasta)

    output = [[0 for i in xrange(len(inFasta)+3)] for j in xrange(len(record)+1)]
    allSum = []

    output[0][0] = "Name"
    output[0][1] = "From"
    output[0][2] = "To"

    storage = []
    for header in inFasta:
        [tmp1, tmp2] = mutCount(ref, inFasta[header], record)
        storage.append([header, tmp1, tmp2])
        allSum.append(tmp1)
    
    totalMedian = median(allSum)

    for m in range(len(storage)):
        output[0][m+3] = storage[m][0]

    for i in range(len(record)):
        output[i+1][0] = record[i][0]
        output[i+1][1] = record[i][1]
        output[i+1][2] = record[i][2]
        for j in range(len(storage)):
            output[i+1][j+3] = int(storage[j][2][i]) * int(totalMedian) * 1.0 / storage[j][1]

    return output


def methodMat(source):
    '''
    This function makes a copy of the SNPmat, maintains the descriptive features
    and set all other values to 0.

    @param source: the SNPmat to copy
    @return a matrix with the values set to 0
    '''
    mat = copy.deepcopy(source)
    for i in range(1,len(mat)):
        for j in range(3, len(mat[0])):
            mat[i][j] = 0
    return mat


def sdSelection(inarr, perc):
    '''
    This function adjust the standard deviation to reduce the effect of outliers.

    @param inarr: input numerical array
    @param perc: the percentage deviation for define a point as not to be included
    toward the standard deviation calculation
    @return the adjusted standard deviation.
    '''
    arr = inarr[:]

    while True:
        if len(arr) <= 2:
            return stdDeviation(inarr)    # Values are too diverse

        sd = stdDeviation(arr)
        if sd == 0:
            return 0.000001                             # Values are the same

        deviation = [stdDeviation(arr[:i]+arr[i+1:])/sd for i in range(len(arr))]
        removal = []
        for i in range(len(deviation)):
            if (deviation[i] > 1+perc) or (deviation[i]<1-perc):
                removal.append(i)

        if len(removal) == 0:
            return sd
        else:
            for index in sorted(removal, reverse=True):
                del arr[index]

    return stdDeviation(arr)

