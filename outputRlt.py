#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ProjectName     : ''
Author          : 'Tian Gao'
CreationDate    : '2017/5/14'
Description     :

"""

import os
import re
import math
import BioinfoComm
import Comm
import SignificaceEvaluation
from Comm import PrintWithTime
import bitarray
from joblib import Parallel, delayed
import logging

INFO = logging.getLogger('root').info

def FetchPatternSetCov(patternSet, userSeqLis, kmer2seqIdInt, seqCnt):
    pattern2cov = {}
    pattern2seq = {}
    for curPattern in patternSet:
        cov, covSeqIdLis, kmer2seqIdInt = BioinfoComm.FetchPatternCov(curPattern, userSeqLis, kmer2seqIdInt)
        pattern2cov[curPattern] = cov
        pattern2seq[curPattern] = covSeqIdLis
    return pattern2cov, pattern2seq
#end_func

# @profile()
def FetchPatternPvalue(minIC, seqCnt, pattern2IC, pattern2cov, sampledMatrix, refKmer2covIdxMatrix, samplingFreq, enableFilter=True, disabledInputPattern=set(), enableParallel=False, coreNm=8):
    pattern2pvalue = {}
    disabledOutputPattern = set()
    if not enableParallel:
        for curPattern, userCov in pattern2cov.iteritems():
            if curPattern in disabledInputPattern: continue
            IC = pattern2IC[curPattern]
            # IC = pattern2IC.get(curPattern, 0)

            # a rough filter of patterns
            # ==TODO
            if enableFilter and (IC < minIC or '[ACTG]' in curPattern or '[ACGU]' in curPattern): continue

            distriLis = FetchPatternCovDistri(curPattern, sampledMatrix, refKmer2covIdxMatrix, seqCnt, samplingFreq)
            # rawPvalue = BioinfoComm.FetchPvalueFromBG(distriLis, userCov)
            rawPvalue = BioinfoComm.FetchPvalueFromBG_FindNonZero(distriLis, userCov) # new pvalue method
            adjustedPvalue = min(rawPvalue * len(pattern2cov), 1.0) # adjust P-value by multiplying the length

            # ==TODO
            # print curPattern
            # print adjustedPvalue
            # ==

            pattern2pvalue[curPattern] = adjustedPvalue
            if enableFilter and adjustedPvalue > 0.05: disabledOutputPattern.add(curPattern) # mark the patterns with Pvalue(RNA) > 0.05 so that we don't need to calculate Pvalue(miR) any more
            PrintWithTime('%s, IC:%s, pvalue:%s, coverage:%s' % (curPattern, IC, adjustedPvalue, userCov), isLogging=True)
    else:
        pvalueInfoLis = Parallel(n_jobs=coreNm)(delayed(ParallelFetchPvalue)(curPattern, userCov, disabledInputPattern, pattern2IC, minIC, sampledMatrix, refKmer2covIdxMatrix, seqCnt, samplingFreq, pattern2cov) for curPattern, userCov in pattern2cov.iteritems())
        for curPattern, adjustedPvalue in pvalueInfoLis:
            if adjustedPvalue == -1: continue
            pattern2pvalue[curPattern] = adjustedPvalue
            if adjustedPvalue > 0.05: disabledOutputPattern.add(curPattern)

    return pattern2pvalue, disabledOutputPattern
#end_func

def ParallelFetchPvalue(curPattern, userCov, disabledInputPattern, pattern2IC, minIC, sampledMatrix,
                        refKmer2covIdxMatrix, seqCnt, samplingFreq, pattern2cov):
    if curPattern in disabledInputPattern: return curPattern, -1
    IC = pattern2IC[curPattern]

    # a rough filter of patterns
    if IC < minIC or '[ACTG]' in curPattern or '[ACUG]' in curPattern: return curPattern, -1

    distriLis = FetchPatternCovDistri(curPattern, sampledMatrix, refKmer2covIdxMatrix, seqCnt, samplingFreq)
    # rawPvalue = BioinfoComm.FetchPvalueFromBG(distriLis, userCov)
    rawPvalue = BioinfoComm.FetchPvalueFromBG_FindNonZero(distriLis, userCov)  # new pvalue method
    adjustedPvalue = min(rawPvalue * len(pattern2cov), 1.0)  # adjust P-value by multiplying the length
    PrintWithTime('%s, IC:%s, pvalue:%s, coverage:%s' % (curPattern, IC, adjustedPvalue, userCov), isLogging=True)
    return curPattern, adjustedPvalue
#end_func

def FetchPatternCovDistri(curPattern, sampledMatrix, refKmer2covIdxMatrix, seqCnt, samplingFreq):
    """
    fetch the distribution list of a certain pattern in the sampled sequences matrix. each kmer has a matrix.
    ***
    to speed up, parse pattern into kmer and store seq id in refKmer2covIdxMatrix
    ***
    """
    covLis = []
    kmerLis = BioinfoComm.FetchAllKmerFromPattern(curPattern)
    for i in xrange(samplingFreq):
        patternCovIdxInt = 0
        for kmer in kmerLis:
            if kmer in refKmer2covIdxMatrix:
                curPatternCovIdxInt = refKmer2covIdxMatrix[kmer][i]
            else:
                covIdxMatrix, covLis = SignificaceEvaluation.FetchCovSeqDetail(sampledMatrix, kmer, seqCnt)
                refKmer2covIdxMatrix[kmer] = covIdxMatrix
                curPatternCovIdxInt = refKmer2covIdxMatrix[kmer][i]
            patternCovIdxInt = patternCovIdxInt | long(curPatternCovIdxInt)
        curBitarray = Comm.int2bitarray(patternCovIdxInt)
        trueCnt = int(bitarray.bitarray.count(curBitarray))
        covLis.append(trueCnt)
    distriLis = FetchCovDistri(covLis, seqCnt)
    return distriLis
#end_func

def FetchCovFromSeqID(allCovIdMatrix):
    """
    merge all the sequence id list in the same sampling together and fetch coverage
    """
    mergeCovLis = []
    matrixCnt = len(allCovIdMatrix)
    samplingFreq = len(allCovIdMatrix[0])
    for i in xrange(samplingFreq):
        allSeqSet = set()
        for j in xrange(matrixCnt):
            curSeqSet = allCovIdMatrix[j][i]
            allSeqSet = allSeqSet | curSeqSet
        cov = len(allSeqSet)
        mergeCovLis.append(cov)
    return mergeCovLis
#end_func

def FetchCovDistri(covLis, samplingSize):
    """
    turn the covLis list into a special form(distribution list) to calculate pvalue and others
    :return
        a dict map from kmer to distribution list
        in the distribution list, index means coverage, value means how many time this coverage occur when sampling
        eg.
        {'ACT': [0, 0, 0, 3, ..., 0]}
    """
    freqOfCertainCoverage = [0] * (samplingSize + 1)
    for cov in covLis:
        freqOfCertainCoverage[cov] += 1
    return freqOfCertainCoverage
#end_func

def OutputPPM(motifLength, patternInfoLis, PMMFn, uniprobeFn, userSeqLis, alphabet, backgroundRatio):
    pattern2PPM = {}
    uniprobeFnLis = []

    # search match of the pattern in user data
    for item in patternInfoLis:
        if len(item) != 5:continue
        pattern, IC, cov, RNAPvalue, miRPvalue = item

        # match pattern to user sequences to get the exact segment to form PWM
        # weight is calculated based on how many times the segments are matched in a sequences
        matchMatrix = [re.findall(pattern, userSeq) for userSeq in userSeqLis]
        seqWeightLis = [1.0 / len(matchLis) if len(matchLis) > 0 else 0 for matchLis in matchMatrix]
        # PPM = [{curChar: 0 for curChar in alphabet}] * motifLength
        PPM = [{curChar: 0 for curChar in alphabet} for _ in xrange(motifLength)]

        for seqIdx, matchLis in enumerate(matchMatrix):
            if len(matchLis) == 0: continue
            weight = seqWeightLis[seqIdx]
            for matchedSegment in matchLis:
                for charIdx, curChar in enumerate(matchedSegment):
                    PPM[charIdx][curChar] += weight

        # output PMM
        fnLis = os.path.splitext(PMMFn)
        curPPMFn = '%s-%s%s' % (fnLis[0], pattern, fnLis[1])
        with open(curPPMFn, 'w') as PPMFileobj:
            headerLine = 'PO\t%s\n' % '\t'.join(alphabet)
            PPMFileobj.write(headerLine)

            for posIdx, posDic in enumerate(PPM):
                PPMValues = map(lambda x:str(posDic[x]), alphabet)
                outputLine = '%s\t%s\n' % (posIdx + 1, '\t'.join(PPMValues))
                PPMFileobj.write(outputLine)

        pattern2PPM[pattern] = PPM

        # output motif in MEME format to run uniprobe2meme and tomtom
        fnLis = os.path.splitext(uniprobeFn)
        curUniprobeFn = '%s-%s%s' % (fnLis[0], pattern, fnLis[1])
        lenPPM = len(PPM)
        uniprobeFnLis.append(curUniprobeFn)
        with open(curUniprobeFn, 'w') as curUniprobeFileobj:
            curUniprobeFileobj.write('%s\n' % pattern)

            for curChar in alphabet:
                charLis = [str(PPM[pos][curChar]) for pos in xrange(lenPPM)]
                outputLine = '%s:\t%s' % (curChar.replace('U', 'T'), '\t'.join(charLis)) # replace U with T to run uniprobe2meme and tomtom
                # outputLine = '%s:\t%s' % (curChar, '\t'.join(charLis)) # replace U with T to run uniprobe2meme and tomtom
                curUniprobeFileobj.write('%s\n' % outputLine)

    return pattern2PPM, uniprobeFnLis
#end_func

def OutputLogo(coreMotifLis, motifLength, ghostscriptProcessor, PPMFn, epsFn, pngFn):
    for coreMotif in coreMotifLis:
        fnLis = os.path.splitext(PPMFn)
        curPMMFn = '%s-%s%s' % (fnLis[0], coreMotif, fnLis[1])

        # output eps
        fnLis = os.path.splitext(epsFn)
        curEpsFn = '%s-%s%s' % (fnLis[0], coreMotif, fnLis[1])
        cmd = 'weblogo --format EPS < "%s" > "%s"' % (curPMMFn, curEpsFn)
        INFO(cmd)
        os.system(cmd)

        # output LOGO
        fnLis = os.path.splitext(pngFn)
        curPngFn = '%s-%s%s' % (fnLis[0], coreMotif, fnLis[1])
        cmd = './%s -dBATCH -dNOPAUSE -dEPSCrop -r300 -sDEVICE=png16m -sOutputFile="%s" "%s"' % (ghostscriptProcessor, curPngFn, curEpsFn)
        INFO(cmd)
        os.system(cmd)
#end_func

def GetPWM(strLis, alphabet, bgRateMap):
    # TODO
    cntMatrix = []
    for curChar in alphabet:
        curCntLis = []
        for curStr in strLis:
            cnt = curStr.count(curChar)
            rate = cnt * 1.0 / len(strLis[0])
            if rate == 0: rate = 0.01
            logRate = math.log(rate * 1.0 / bgRateMap[curChar], 2)
            logRate = '%.2f' % logRate # only keep the first two digits after point
            curCntLis.append(logRate)
        cntMatrix.append(curCntLis)
    return cntMatrix
#end_func
