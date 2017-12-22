#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/5/14
Description     :

"""
import BioinfoComm, Comm
import random
import DimerGraph
import bitarray
import numpy as np
from joblib import Parallel, delayed

def SampleRNARef(referenceFilenm, samplingSize, samplingFreq, minLen=18, maxLen=25, alphabet='ACGT'):
    """
    :return
        return sampled sequences from RNA reference
    """
    refSegLis = []
    refSeqLis, seqCnt, titleLis = BioinfoComm.loadMultilineSeq(referenceFilenm, maxLen)
    for refSeq in refSeqLis:
        samplingTestCnt = 0
        segLength = random.randint(minLen, maxLen)  # original length of the segment is 18 ~ 25
        while len(refSeq) - segLength - 1 <= 0 and samplingTestCnt < 3:
            segLength = random.randint(minLen, maxLen)
            samplingTestCnt += 1
        if samplingTestCnt >= 3: continue

        startPos = random.randint(0, len(refSeq) - segLength - 1)
        geneSeg = refSeq[startPos: startPos + segLength]
        if alphabet == 'ACGU': geneSeg = geneSeg.replace('T', 'U')
        while 'N' in geneSeg and samplingTestCnt < 5:  # if there is 'N' in the segment, redo sampling
            startPos = random.randint(0, len(refSeq) - segLength - 1)
            geneSeg = refSeq[startPos: startPos + segLength]
            samplingTestCnt += 1
            if alphabet == 'ACGU': geneSeg = geneSeg.replace('T', 'U')
        if samplingTestCnt >= 5: continue

        # shuffle here
        strLis = list(geneSeg)
        random.shuffle(strLis)
        geneSeg = ''.join(strLis)

        refSegLis.append(geneSeg)

    refSeqLis = refSegLis
    sampledSeqMatrix = [random.sample(refSeqLis, samplingSize) for _ in range(samplingFreq)]
    return sampledSeqMatrix
# end_func

def SampleMiRRef(referenceFilenm, samplingSize, samplingFreq, alphabet='ACGT'):
    """
    :return
       return sampled sequences from miRNA reference
    """
    refSeqLis, seqCnt, _, _ = BioinfoComm.loadSinglelineSeq(referenceFilenm)
    if alphabet == 'ACGU': refSeqLis = map(lambda x:x.replace('T', 'U'), refSeqLis)
    sampledSeqMatrix = [random.sample(refSeqLis, samplingSize) for _ in range(samplingFreq)]
    return sampledSeqMatrix
#end_func

# @profile
def EvalSignificance(K, sampledSeqMatrix, lastLayerKmerLis, samplingSize, samplingFreq, pathDic, inTreeKmerSet, alphabet, targetKmerLis=[], enableParallel=False, coreNm=8, covPrecThres=0.1):
    userKmer2Cov, kmer2VisDetail = DimerGraph.PathDic2PathCov(pathDic, K)

    allKmerSet = BioinfoComm.GenerateKmerLis(alphabet, K)
    enabledKmerSet = [kmer for kmer in allKmerSet if userKmer2Cov.get(kmer, 0) > covPrecThres * samplingSize or kmer in targetKmerLis]

    # === TODO
    # enabledKmerSet = {'AAAAAC', 'AAAAAG', 'AAAATC', 'AAAATG', 'AAACAC', 'AAACAG', 'AAACTC', 'AAACTG', 'AAAGAC',
    # 'AAAGAG', 'AAAGTC', 'AAAGTG', 'AAGAAC', 'AAGAAG', 'AAGATC', 'AAGATG', 'AAGCAC', 'AAGCAG', 'AAGCTC', 'AAGCTG',
    # 'AAGGAC', 'AAGGAG', 'AAGGTC', 'AAGGTG', 'AATAAC', 'AATAAG', 'AATATC', 'AATATG', 'AATCAC', 'AATCAG', 'AATCTC',
    # 'AATCTG', 'AATGAC', 'AATGAG', 'AATGTC', 'AATGTG', 'AGAAAC', 'AGAAAG', 'AGAATC', 'AGAATG', 'AGACAC', 'AGACAG',
    # 'AGACTC', 'AGACTG', 'AGAGAC', 'AGAGAG', 'AGAGTC', 'AGAGTG', 'AGGAAC', 'AGGAAG', 'AGGATC', 'AGGATG', 'AGGCAC',
    # 'AGGCAG', 'AGGCTC', 'AGGCTG', 'AGGGAC', 'AGGGAG', 'AGGGTC', 'AGGGTG', 'AGTAAC', 'AGTAAG', 'AGTATC', 'AGTATG',
    # 'AGTCAC', 'AGTCAG', 'AGTCTC', 'AGTCTG', 'AGTGAC', 'AGTGAG', 'AGTGTC', 'AGTGTG', 'ATAAAC', 'ATAAAG', 'ATAATC',
    # 'ATAATG', 'ATACAC', 'ATACAG', 'ATACTC', 'ATACTG', 'ATAGAC', 'ATAGAG', 'ATAGTC', 'ATAGTG', 'ATGAAC', 'ATGAAG',
    # 'ATGATC', 'ATGATG', 'ATGCAC', 'ATGCAG', 'ATGCTC', 'ATGCTG', 'ATGGAC', 'ATGGAG', 'ATGGTC', 'ATGGTG', 'ATTAAC',
    # 'ATTAAG', 'ATTATC', 'ATTATG', 'ATTCAC', 'ATTCAG', 'ATTCTC', 'ATTCTG', 'ATTGAC', 'ATTGAG', 'ATTGTC', 'ATTGTG',
    # 'GAAAAC', 'GAAAAG', 'GAAATC', 'GAAATG', 'GAACAC', 'GAACAG', 'GAACTC', 'GAACTG', 'GAAGAC', 'GAAGAG', 'GAAGTC',
    # 'GAAGTG', 'GAGAAC', 'GAGAAG', 'GAGATC', 'GAGATG', 'GAGCAC', 'GAGCAG', 'GAGCTC', 'GAGCTG', 'GAGGAC', 'GAGGAG',
    # 'GAGGTC', 'GAGGTG', 'GATAAC', 'GATAAG', 'GATATC', 'GATATG', 'GATCAC', 'GATCAG', 'GATCTC', 'GATCTG', 'GATGAC',
    # 'GATGAG', 'GATGTC', 'GATGTG', 'GGAAAC', 'GGAAAG', 'GGAATC', 'GGAATG', 'GGACAC', 'GGACAG', 'GGACTC', 'GGACTG',
    # 'GGAGAC', 'GGAGAG', 'GGAGTC', 'GGAGTG', 'GGGAAC', 'GGGAAG', 'GGGATC', 'GGGATG', 'GGGCAC', 'GGGCAG', 'GGGCTC',
    # 'GGGCTG', 'GGGGAC', 'GGGGAG', 'GGGGTC', 'GGGGTG', 'GGTAAC', 'GGTAAG', 'GGTATC', 'GGTATG', 'GGTCAC', 'GGTCAG',
    # 'GGTCTC', 'GGTCTG', 'GGTGAC', 'GGTGAG', 'GGTGTC', 'GGTGTG', 'GTAAAC', 'GTAAAG', 'GTAATC', 'GTAATG', 'GTACAC',
    # 'GTACAG', 'GTACTC', 'GTACTG', 'GTAGAC', 'GTAGAG', 'GTAGTC', 'GTAGTG', 'GTGAAC', 'GTGAAG', 'GTGATC', 'GTGATG',
    # 'GTGCAC', 'GTGCAG', 'GTGCTC', 'GTGCTG', 'GTGGAC', 'GTGGAG', 'GTGGTC', 'GTGGTG', 'GTTAAC', 'GTTAAG', 'GTTATC',
    # 'GTTATG', 'GTTCAC', 'GTTCAG', 'GTTCTC', 'GTTCTG', 'GTTGAC', 'GTTGAG', 'GTTGTC', 'GTTGTG', 'TAAAAC', 'TAAAAG',
    # 'TAAATC', 'TAAATG', 'TAACAC', 'TAACAG', 'TAACTC', 'TAACTG', 'TAAGAC', 'TAAGAG', 'TAAGTC', 'TAAGTG', 'TAGAAC',
    # 'TAGAAG', 'TAGATC', 'TAGATG', 'TAGCAC', 'TAGCAG', 'TAGCTC', 'TAGCTG', 'TAGGAC', 'TAGGAG', 'TAGGTC', 'TAGGTG',
    # 'TATAAC', 'TATAAG', 'TATATC', 'TATATG', 'TATCAC', 'TATCAG', 'TATCTC', 'TATCTG', 'TATGAC', 'TATGAG', 'TATGTC',
    # 'TATGTG', 'TGAAAC', 'TGAAAG', 'TGAATC', 'TGAATG', 'TGACAC', 'TGACAG', 'TGACTC', 'TGACTG', 'TGAGAC', 'TGAGAG',
    # 'TGAGTC', 'TGAGTG', 'TGGAAC', 'TGGAAG', 'TGGATC', 'TGGATG', 'TGGCAC', 'TGGCAG', 'TGGCTC', 'TGGCTG', 'TGGGAC',
    # 'TGGGAG', 'TGGGTC', 'TGGGTG', 'TGTAAC', 'TGTAAG', 'TGTATC', 'TGTATG', 'TGTCAC', 'TGTCAG', 'TGTCTC', 'TGTCTG',
    # 'TGTGAC', 'TGTGAG', 'TGTGTC', 'TGTGTG', 'TTAAAC', 'TTAAAG', 'TTAATC', 'TTAATG', 'TTACAC', 'TTACAG', 'TTACTC',
    # 'TTACTG', 'TTAGAC', 'TTAGAG', 'TTAGTC', 'TTAGTG', 'TTGAAC', 'TTGAAG', 'TTGATC', 'TTGATG', 'TTGCAC', 'TTGCAG',
    # 'TTGCTC', 'TTGCTG', 'TTGGAC', 'TTGGAG', 'TTGGTC', 'TTGGTG', 'TTTAAC', 'TTTAAG', 'TTTATC', 'TTTATG', 'TTTCAC',
    # 'TTTCAG', 'TTTCTC', 'TTTCTG', 'TTTGAC', 'TTTGAG', 'TTTGTC', 'TTTGTG'}
    # === TODO

    kmer2refDistri, kmer2refCovThreshold, refKmer2covIdxMatrix, signKmerSet = fetchSignKmerSet(K, sampledSeqMatrix, enabledKmerSet, samplingSize, userKmer2Cov, enableParallel=enableParallel, coreNm=coreNm)

    #TODO:
    # signKmerSet = enabledKmerSet

    if not signKmerSet: return signKmerSet, userKmer2Cov, kmer2VisDetail, kmer2refCovThreshold, {}, {}, {}, set()
    kmer2pvalue = FetchKmerPvalue(signKmerSet, kmer2refDistri, userKmer2Cov, samplingFreq)

    if K == 2:
        inTreeKmerSet = set(lastLayerKmerLis)
    else:
        for enabledKmer in lastLayerKmerLis:
            for curChar in alphabet:
                inTreeKmerSet.add('%s%s' % (curChar, enabledKmer))
                inTreeKmerSet.add('%s%s' % (enabledKmer, curChar))

    return signKmerSet, userKmer2Cov, kmer2VisDetail, kmer2refCovThreshold, kmer2pvalue, kmer2refDistri, refKmer2covIdxMatrix, inTreeKmerSet
#end_func

def fetchSignKmerSet(K, sampledSeqMatrix, enabledKmerSet, samplingSize, userKmer2Cov, distriShreshold=0.99, enableParallel=False, coreNm=8):
    """
    fetch significant kmer and their coverage distribution, threshold, coverage detail
    """
    # fetch significant kmer
    kmer2distri = {}
    kmer2covThreshold = {}
    refKmer2covIdxMatrix = {}
    signKmerSet = set()

    if not enableParallel:
        for kmer in enabledKmerSet:
            covIdxMatrix, covLis = FetchCovSeqDetail(sampledSeqMatrix, kmer, samplingSize)
            distriLis = FetchCovDistri(covLis, samplingSize)
            covShreshold = BioinfoComm.FetchVarConfIdx(distriLis, distriShreshold)
            userCov = userKmer2Cov.get(kmer, 0)
            kmer2distri[kmer] = distriLis
            refKmer2covIdxMatrix[kmer] = covIdxMatrix
            kmer2covThreshold[kmer] = covShreshold
            if userCov > covShreshold:
                signKmerSet.add(kmer)
    else:
        signRltLis = Parallel(n_jobs=coreNm)(delayed(ParallelSignKmer)(kmer, sampledSeqMatrix, samplingSize, distriShreshold, userKmer2Cov) for kmer in enabledKmerSet)
        for kmer, signRlt, distriLis, covIdxMatrix, covShreshold in signRltLis:
            kmer2distri[kmer] = distriLis
            kmer2covThreshold[kmer] = covShreshold
            refKmer2covIdxMatrix[kmer] = covIdxMatrix
            if signRlt:
                signKmerSet.add(kmer)

    return kmer2distri, kmer2covThreshold, refKmer2covIdxMatrix, signKmerSet
# end_func

def ParallelSignKmer(kmer, sampledSeqMatrix, samplingSize, distriShreshold, userKmer2Cov):
    covIdxMatrix, covLis = FetchCovSeqDetail(sampledSeqMatrix, kmer, samplingSize)
    distriLis = FetchCovDistri(covLis, samplingSize)
    covShreshold = BioinfoComm.FetchVarConfIdx(distriLis, distriShreshold)
    userCov = userKmer2Cov[kmer]
    signRlt = True if userCov > covShreshold else False
    return kmer, signRlt, distriLis, covIdxMatrix, covShreshold
#end_func

def FetchCovSeqDetail(seqMatrix, kmer, samplingSize):
    """
    fetch sequence id in seqMatrix covered by kmer     
    """
    covIdxMatrix = []
    covLis = []
    for seqLis in seqMatrix:
        covIdxLis = {seqId for seqId, seq in enumerate(seqLis) if kmer in seq}
        covIdxBitarray = BioinfoComm.IdxLis2bin(covIdxLis, samplingSize)
        covIdxInt = Comm.bitarray2int(covIdxBitarray)
        covIdxMatrix.append(covIdxInt)
        covLis.append(int(bitarray.bitarray.count(covIdxBitarray)))
    # to save memory
    covIdxMatrix = np.array(covIdxMatrix)
    return covIdxMatrix, covLis
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

def FetchKmerPvalue(possKmerSet, bgDistri, userKmer2Cov, samplingFreq):
    kmer2pvalue = {}
    for kmer in possKmerSet:
        distriLis = bgDistri[kmer]
        userCov = userKmer2Cov.get(kmer, 0)
        pvalue = BioinfoComm.FetchPvalueFromBG_FindNonZero(distriLis, userCov)
        kmer2pvalue[kmer] = pvalue
    return kmer2pvalue
#end_func
