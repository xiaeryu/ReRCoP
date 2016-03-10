import math
import random
from scipy import stats
import ReRCoP_matrix


def kNN(inarr, Pk, Pthreshold):
    '''
    This function implements the kNN outlier detection.

    @param inarr: input numerical array
    @param Pk: k value in the kNN method input as a percentage of the number 
    of points
    @param Pthreshold: threshold used to define a neighbor input as a factor
    of the standard deviation
    @return a list of indices of the outliers
    '''
    total = len(inarr)		# Total number of points
    sd = ReRCoP_matrix.stdDeviation(inarr)
    threshold = Pthreshold * sd	# Threshold used to define a neighbor
    k = int(total * Pk)		# k value in the kNN method

    outlier = []
    for i in range(total):
        count = 0
        for j in range(total):
            if abs(inarr[i] - inarr[j]) <= threshold:
                count += 1
                if count >=k:
                    break
        if count < k:
            outlier.append(i)
    return outlier


def Grubbs(inarr, alpha):
    '''
    This function implements the Grubb's outlier detectoin.

    @param inarr: input numerical array
    @param alpha: significance level for the statistical test
    @return a list of indices of the outliers
    '''
    avg = ReRCoP_matrix.average(inarr)
    std = ReRCoP_matrix.stdDeviation(inarr)
    if std == 0:
        return []
    G = [abs(item - avg)/std for item in inarr]
    N = len(inarr)
    t = stats.t.isf(1-alpha/(2*N), N-2)
    Gtest = (N-1)/math.sqrt(N) * math.sqrt(t**2 / (N-2+t**2))
    outlier = []
    for i in range(len(G)):
        if G[i] > Gtest:
            outlier.append(i)
    return outlier


def region_query(inarr, point, eps):
    '''
    This is a helper function for 'DBSCAN'
    '''
    nP = len(inarr)
    seeds = []
    for i in range(nP):
        if abs(inarr[point] - inarr[i]) < eps:
            seeds.append(i)
    return seeds


def expand_cluster(inarr, classification, point, cluster, eps, minP):
    '''
    This is a helper function for 'DBSCAN'
    '''
    seeds = region_query(inarr, point, eps)
    if len(seeds) < minP:
        classification[point] = -1
        return False
    else:
        classification[point] = cluster
        for seed in seeds:
            classification[seed] = cluster

        while len(seeds) > 0:
            currentP = seeds[0]
            results = region_query(inarr, currentP, eps)
            if len(results) >= minP:
                for result in results:
                    if classification[result] in [0,-1]:
                        if classification[result] == 0:
                            seeds.append(result)
                    classification[result] = cluster
            seeds = seeds[1:]
        return True


def DBSCAN(inarr, Peps, PminP):
    '''
    This function implements the DBSCAN outlier detection.

    @param inarr: input numerical array
    @param Peps: eps value in the DBSCAN method input as a factor of the
    standard deviation
    @param PminP: minP value in the DBSCAN method input as a percentage of
    the total number of points
    @return a list of indices of the outliers
    '''
    cluster = 1
    totalP = len(inarr)
    minP = int(totalP * PminP)
    sd = ReRCoP_matrix.sdSelection(inarr, 0.15)       # Estimated sd of the non-outlier points
    eps = Peps * sd
    classification = [0 for i in xrange(totalP)]
    for point in range(totalP):
        if classification[point] == 0:
            if expand_cluster(inarr, classification, point, cluster, eps, minP):
                cluster += 1
    outlier = []
    for i in range(len(classification)):
        if classification[i] == -1:
            outlier.append(i)
    return outlier

