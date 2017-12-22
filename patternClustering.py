#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/5/24
Description     :

"""

import math
import os
import BioinfoComm

def GetCovSimi(pattern1, pattern2, userKmer2seqIdSet):
    seqIdSet1 = userKmer2seqIdSet[pattern1]
    seqIdSet2 = userKmer2seqIdSet[pattern2]
    simi = len(seqIdSet1.intersection(seqIdSet2)) * 1.0 / len(seqIdSet1 | seqIdSet2)
    return simi
#end_func

def outputPatternCov(patternInfoLis, pattern2covFn):
    pattern2cov = {}
    with open(pattern2covFn, 'w') as pattern2covObj:
        for patternItems in patternInfoLis:
            pattern = patternItems[0]
            cov = patternItems[2]
            pattern2cov[pattern] = cov
            pattern2covOutputLine = '%s\t%s' % (pattern, cov)
            pattern2covObj.write('%s\n' % pattern2covOutputLine)
    return pattern2cov
#end_func

def CalcPatternCovSimi(patternInfoLis, userKmer2seqIdSet, pattern2simiFn):
    pattern2covSimi = {}
    with open(pattern2simiFn, 'w') as pattern2simiObj:
        for patternItems1 in patternInfoLis:
            for patternItems2 in patternInfoLis:
                pattern1 = patternItems1[0]
                pattern2 = patternItems2[0]
                if pattern1 <= pattern2: continue
                simi = GetCovSimi(pattern1, pattern2, userKmer2seqIdSet)
                pattern2covSimi[(pattern1, pattern2)] = simi
                pattern2simiOutputLine = '%s\t%s\t%s' % (pattern1, pattern2, simi)
                pattern2simiObj.write('%s\n' % pattern2simiOutputLine)
    return pattern2covSimi
#end_func

def CalcPatternPwmSimi(patternInfoLis, pattern2PPM, patternPairInfo, pattern2simiFn, alphabet, doNormalization=True):
    pattern2pwmSimi = CalcPwmSimi(pattern2PPM, patternPairInfo)
    if doNormalization: pattern2pwmSimi = NormalizeSimi(pattern2pwmSimi)
    with open(pattern2simiFn, 'w') as pattern2simiObj:
        for (pattern1, pattern2), simi in pattern2pwmSimi.iteritems():
            pattern2simiOutputLine = '%s\t%s\t%s' % (pattern1, pattern2, simi)
            pattern2simiObj.write('%s\n' % pattern2simiOutputLine)
    return pattern2pwmSimi
#end_func

def loadSinglePWM(pwmFn):
    """
    load a single PWM as a matrix
    """
    pwmKmerMatrix = []
    with open(pwmFn) as pwmF:
        for line in pwmF:
            line = line.strip()
            if not line or line[0] == '>': continue
            pwmKmerMatrix.append(line)
    return pwmKmerMatrix
#end_func

def loadAllPMM(patternInfoLis, PMMDir, alphabet):
    """
    load all the PWM file in the dir
    colCnt[0] = [0, 312, 0, 138] means that in pos 0, there are 0 A, 312 C, 0 G, 138 T
    """
    pattern2pwm = {}
    pattern2colCnt = {}
    for patternItems in patternInfoLis:
        pattern = patternItems[0]
        pwmFn = os.path.join(PMMDir, 'weighted_kmer-%s.fa' % pattern)
        pwmKmerMatrix = loadSinglePWM(pwmFn)
        pattern2colCnt[pattern] = []
        for colId in range(len(pwmKmerMatrix[0])):
            strInCol = map(lambda x:x[colId], pwmKmerMatrix)
            colCntLis = [strInCol.count(curChar) for curChar in alphabet]
            pattern2colCnt[pattern].append(colCntLis)
        pattern2pwm[pattern] = pwmKmerMatrix
    return pattern2pwm, pattern2colCnt
#end_func

def CalcPwmSimi(pattern2PPM, patternPairInfo):
    """
    get the similarity of each pattern pair
    """
    pattern2Simi = {}
    for patternX, PPMx in pattern2PPM.iteritems():
        for patternY, PPMy in pattern2PPM.iteritems():
            if patternX <= patternY:continue

            offset, pvalue, evalue, qvalue = patternPairInfo.get((patternX, patternY), (0, 1, 1, 1)) # fetch tomtom result from dict
            (shiftedPatternPPM, refPatternPPM, offset) = (PPMx, PPMy, offset) if offset > 0 else (PPMy, PPMx, -offset)

            # add empty slot for offset
            emptyHeader = []
            for _ in range(offset):
                emptyHeader.append({})
            shiftedPatternPPM = emptyHeader + shiftedPatternPPM
            nonEmptyLen = len(shiftedPatternPPM) - offset

            # simi = 0 if pvalue > 0.05 else PCC(shiftedPatternPPM, refPatternPPM)
            # simi = 0 if qvalue > 0.05 else PCC(shiftedPatternPPM, refPatternPPM)
            simi = PCC(shiftedPatternPPM, refPatternPPM) * 1.0 / nonEmptyLen # simi is the avg of PCC distance
            pattern2Simi[(patternX, patternY)] = simi
    return pattern2Simi
#end_func

def PCC(PPMx, PPMy):
    """
    Pearson correlation coefficient
    """
    colPCCLis = []

    for colId, posDicX in enumerate(PPMx):
        if colId >= len(PPMy): continue # PPMx is longer because of offset
        posDicY = PPMy[colId]
        if not posDicX or not posDicY: continue # simi=0 when a certain position is empty because of offset

        colProbX = [v * 1.0 / sum(posDicX.values()) for k, v in posDicX.iteritems()]
        colProbY = [v * 1.0 / sum(posDicY.values()) for k, v in posDicY.iteritems()]

        avgX = 1.0 / len(colProbX)
        avgY = 1.0 / len(colProbY)

        sdXLis = [(Xa - avgX) ** 2 for Xa in colProbX]
        sdX = math.sqrt(sum(sdXLis))
        sdYLis = [(Ya - avgY) ** 2 for Ya in colProbY]
        sdY = math.sqrt(sum(sdYLis))

        denominator = sdX * sdY
        numerator = sum([(Xa - avgX) * (colProbY[idx] - avgY) for idx, Xa in enumerate(colProbX)])

        colPCC = numerator * 1.0 / denominator
        colPCCLis.append(colPCC)
    #end_func

    simi = sum(colPCCLis)
    return simi
#end_func

def NormalizeSimi(pattern2Simi, sudoCount=0.0001):
    patternSimiLis = [(pattern1, pattern2, simi) for (pattern1, pattern2), simi in pattern2Simi.iteritems()]
    simiLis = map(lambda x:x[2], patternSimiLis)
    minSimi = min(simiLis)
    maxSimi = max(simiLis)
    patternSimiLis = [(pattern1, pattern2, (simi - minSimi) * 1.0 / (maxSimi - minSimi)) for pattern1, pattern2, simi in patternSimiLis]
    for idx, (pattern1, pattern2, simi) in enumerate(patternSimiLis):
        if simi == 0: patternSimiLis[idx] = (pattern1, pattern2, simi + sudoCount)
        if simi == 1: patternSimiLis[idx] = (pattern1, pattern2, simi - sudoCount)

    pattern2Simi = {(pattern1, pattern2): simi for pattern1, pattern2, simi in patternSimiLis}
    return pattern2Simi
#end_func

def patternClustering(pattern2cov, pattern2pwmSimi, sudoCount=0.0001):
    # print pattern2cov
    # print '--'
    # print pattern2pwmSimi

    pattern1 = None

    # use mediate as a filter
    simiLis = sorted(pattern2pwmSimi.values())
    med = max(simiLis[len(simiLis) / 2], sudoCount) # simi must be positive

    # build graph, no more filter since simi is already filtered by tomtom significance
    g = igraph.Graph(directed=False)
    for pattern, cov in pattern2cov.iteritems():
        g.add_vertex(pattern, weight=cov)
    for (pattern1, pattern2), simi in pattern2pwmSimi.iteritems():
        # if simi < med: continue
        if simi < 0: continue
        g.add_edge(pattern1, pattern2, weight=simi)

    # print g.get_edgelist()

    # do clustering
    clusterLis = []
    if len(g.es) == 0: return []
    clusterRlt = g.community_multilevel(g.es['weight'], return_levels=False)
    # print clusterRlt
    clusterInfoStr = str(clusterRlt).replace(',\n', ', ')
    allCluster = clusterInfoStr.split('\n')[1:]
    for singleClusterStr in allCluster:
        clusterItems = singleClusterStr.partition('] ')[2].split(', ')
        if clusterItems == ['None']: continue
        clusterItems = map(lambda x:x.strip(), clusterItems)
        clusterLis.append(set(clusterItems))

    # if BioinfoComm.getPatternLength(pattern1) == 4:
    #     print len(clusterLis)
    #     exit()
    return clusterLis
#end_func


def loadClusterRlt(clusterRltFn):
    """
    load the cluster result from txt file processed by R script
    """
    clusterLis = []
    with open(clusterRltFn) as clusterRltFileobj:
        for line in clusterRltFileobj:
            elementsInOneCluster = set(line.strip().split('\t'))
            clusterLis.append(elementsInOneCluster)
    return clusterLis
#end_func


def findCore(clusterLis, userKmer2covIdxInt, seqCnt, pattern2IC, RNAPattern2pvalue, miRPattern2pvalue, maxNodeCnt=5, fetchGroupMode=0, covPrecThres=0.1):
    """
    :param fetchGroupMode:
        0 - original method, only consider coverage in each cluster. Each cluster has a motif
        1 - only consider pvalue in each cluster. Each cluster has a motif
        2 - merge all one-node cluster and handle it using method 0. Then only consider multi-node cluster
    """
    clusterCoreLis = []
    if fetchGroupMode == 0: # based on coverage
        for cluster in clusterLis:
            # looking for the 1st node, which has the largest coverage
            clusterCoreLis.append([])
            clusterSeqIdSet = set()
            maxCovPattern = ''
            maxCovValue = 0
            for pattern in cluster: # loop to find the pattern with largest coverage
                covIdxInt = userKmer2covIdxInt[pattern]
                covIdxSet = BioinfoComm.covIdxInt2covIdxSet(covIdxInt)
                if (len(covIdxSet) > maxCovValue) or (len(covIdxSet) == maxCovValue and pattern2IC[pattern] > pattern2IC[maxCovPattern]):
                    maxCovValue = len(covIdxSet)
                    maxCovPattern = pattern
                clusterSeqIdSet |= covIdxSet
            cluster.remove(maxCovPattern)
            restSeqIdSet = clusterSeqIdSet - BioinfoComm.covIdxInt2covIdxSet(userKmer2covIdxInt[maxCovPattern])
            clusterCoreLis[-1].append([maxCovPattern, len(clusterSeqIdSet), len(restSeqIdSet)])

            # looking for the nodes that cover most sequence in the rest set
            while cluster and restSeqIdSet and len(clusterCoreLis[-1]) < maxNodeCnt:
                nextPattern = cluster.pop()
                cluster.add(nextPattern)
                nextRestCov = len(restSeqIdSet - BioinfoComm.covIdxInt2covIdxSet(userKmer2covIdxInt[nextPattern]))
                for pattern in cluster:
                    restCov = len(restSeqIdSet - BioinfoComm.covIdxInt2covIdxSet(userKmer2covIdxInt[pattern]))
                    if restCov < nextRestCov:  # cover more sequences in the rest
                        nextPattern = pattern
                        nextRestCov = restCov
                    elif restCov == nextRestCov and BioinfoComm.covIdxInt2covCnt(userKmer2covIdxInt[pattern]) > BioinfoComm.covIdxInt2covCnt(userKmer2covIdxInt[nextPattern]):
                        # cover equal sequences in the rest but have better overall coverage
                        nextPattern = pattern
                if len(restSeqIdSet) - nextRestCov < seqCnt * 0.1: break # if coverage 'step' is tiny, stop the loop
                clusterCoreLis[-1].append([nextPattern, len(restSeqIdSet), nextRestCov])
                cluster.remove(nextPattern)
                restSeqIdSet = restSeqIdSet - BioinfoComm.covIdxInt2covIdxSet(userKmer2covIdxInt[nextPattern])
        return clusterCoreLis
    elif fetchGroupMode == 1: # based on pvalue
        for cluster in clusterLis:
            curCoreLis = []
            for pattern in cluster:
                RNAPvalue = RNAPattern2pvalue[pattern]
                miRPvalue = miRPattern2pvalue[pattern]
                curCoreLis.append((pattern, RNAPvalue, miRPvalue))
            curCoreLis = sorted(curCoreLis, key=lambda x:x[1])[:maxNodeCnt]
            clusterCoreLis.append(curCoreLis)
        return clusterCoreLis
    elif fetchGroupMode == 2: # reconstruct cluster list and call findCore. In this mode, all one-node cluster is merged as one cluster. Others are not changed.
        newClusterLis = []
        oneNodeCluter = set()
        for cluster in clusterLis:
            if len(cluster) > 1:
                newClusterLis.append(cluster)
            else:
                node = cluster.pop()
                oneNodeCluter.add(node)
        if oneNodeCluter:
            newClusterLis.append(oneNodeCluter)
        return findCore(newClusterLis, userKmer2covIdxInt, seqCnt, pattern2IC, RNAPattern2pvalue, miRPattern2pvalue, maxNodeCnt, fetchGroupMode=0, covPrecThres=covPrecThres)
    elif fetchGroupMode == 3: # In this mode, all one-node cluster is merged as one cluster. all multi-node cluster is merged as another cluster.
        newClusterLis = []
        oneNodeCluter = set()
        multinodeCluster = set()
        for cluster in clusterLis:
            if len(cluster) > 1:
                # newClusterLis.append(cluster)
                multinodeCluster = multinodeCluster | cluster
            else:
                node = cluster.pop()
                oneNodeCluter.add(node)
        if multinodeCluster:
            newClusterLis.append(multinodeCluster)
        if oneNodeCluter:
            newClusterLis.append(oneNodeCluter)

        return findCore(newClusterLis, userKmer2covIdxInt, seqCnt, pattern2IC, RNAPattern2pvalue, miRPattern2pvalue, maxNodeCnt, fetchGroupMode=0, covPrecThres=covPrecThres)
#end_func


def searchForCoreInCluster(clusterLis, userKmer2seqIdSet, seqCnt, pattern2cov, pattern2IC, RNAPattern2pvalue, miRPattern2pvalue, fetchGroupMode=0, covPrecThres=0.1):
    clusterCoreLis = findCore(clusterLis, userKmer2seqIdSet, seqCnt, pattern2IC, RNAPattern2pvalue, miRPattern2pvalue, fetchGroupMode=fetchGroupMode, covPrecThres=covPrecThres)
    # print '**' * 20
    # print clusterCoreLis
    # print '**' * 20
    # print clusterCoreLis[0]
    # print '**' * 20

    coreMotifSet = set()
    motifCandSet = set()
    coreOutputLines = []
    for clusterCores in clusterCoreLis:
        for idx, (core, totalSeq, restSeq) in enumerate(clusterCores):
            if idx == 0:
                coreMotifSet.add(core)
            else:
                motifCandSet.add(core)
            outputLis = [core, str(pattern2cov[core]), str(totalSeq), str(restSeq)]
            outputLine = '\t'.join(outputLis)
            coreOutputLines.append(outputLine)
        coreOutputLines.append('==')

    uncoveredSeqIdSet = set(range(seqCnt)) - reduce(lambda x, y: x | y, map(lambda x: BioinfoComm.covIdxInt2covIdxSet(userKmer2seqIdSet[x]), coreMotifSet))
    # the top pattern in each cluster already cover all sequences
    while motifCandSet and uncoveredSeqIdSet:
        bestMotifCand = motifCandSet.pop()
        motifCandSet.add(bestMotifCand)
        bestCov = uncoveredSeqIdSet - BioinfoComm.covIdxInt2covIdxSet(userKmer2seqIdSet[bestMotifCand])
        for curMotifCand in motifCandSet:
            curCov = uncoveredSeqIdSet - BioinfoComm.covIdxInt2covIdxSet(userKmer2seqIdSet[curMotifCand])
            if len(curCov) < len(bestCov):  # cover more sequences in the rest
                bestMotifCand = curMotifCand
                bestCov = curCov
            elif len(curCov) == len(bestCov) and pattern2cov[curMotifCand] > pattern2cov[bestMotifCand]:
                bestMotifCand = curMotifCand  # cover equal sequences in the rest but have better overall coverage
        newUncoveredSeqIdSet = uncoveredSeqIdSet - BioinfoComm.covIdxInt2covIdxSet(userKmer2seqIdSet[bestMotifCand])
        if len(newUncoveredSeqIdSet) == len(uncoveredSeqIdSet): break
        uncoveredSeqIdSet = newUncoveredSeqIdSet
        coreMotifSet.add(bestMotifCand)
        motifCandSet.remove(bestMotifCand)

    return coreMotifSet, coreOutputLines
#end_func