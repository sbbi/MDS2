#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/6/16
Description     :

"""
import BioinfoComm
import scipy.stats.stats as stats
import logging

INFO = logging.getLogger('root').info

def TestPatternInNegSeq(posPatternCovLis, posSeqCnt, negSeqFn, allKmerSet):
    patternSet = set()
    patternCnt = len(posPatternCovLis)
    initPatternSet = map(lambda x:x[0], posPatternCovLis)
    negSeqLis, negSeqCnt, _, _ = BioinfoComm.loadSinglelineSeq(negSeqFn)
    negKmer2seqIdSet = BioinfoComm.FetchCovInSeqLisMutliKmer(negSeqLis, allKmerSet)
    negKmer2seqIdInt = BioinfoComm.formatCovId(negKmer2seqIdSet, negSeqCnt)

    for pattern, posCov in posPatternCovLis:
        posUncov = posSeqCnt - posCov
        negCov, _, negKmer2seqIdInt = BioinfoComm.FetchPatternCov(pattern, negSeqLis, negKmer2seqIdInt)
        negUncov = negSeqCnt - negCov
        dataTable = [[posCov, posUncov], [negCov, negUncov]]
        _, rawPValue = stats.fisher_exact(dataTable)
        adjustedPvalue = min(rawPValue * patternCnt, 1)
        if adjustedPvalue < 0.05: patternSet.add(pattern)

    INFO('pattern before negative filter')
    INFO(initPatternSet)
    INFO('pattern after negative filter')
    INFO(patternSet)
    return patternSet
#end_func
