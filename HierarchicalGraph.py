#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ProjectName     : ''
Author          : 'Tian Gao'
CreationDate    : '2017/5/14'
Description     :

"""
import os
from Comm import print0, print9
from graphviz import Digraph
import logging

INFO = logging.getLogger('root').info

def BuildTree(K, wholeTree, sigKmerSet, lastLayerKmerLis, userCoverageDict, userSeqLis, covPrecThres=0.1):
    """
    generate necessary information so that developing tree can be built in the following steps  
    """
    filtSigKmerSet = set()
    kmer2parentDic = {}
    curLayerP2Kmers = {}

    filtSigKmerSet.add('_' * K)
    if K > 2:
        parent = '_' * (K - 1)
        curLayerP2Kmers.setdefault(parent, [set(), set()])
        curLayerP2Kmers[parent][0].add('_' * K)

    for kmer in sigKmerSet:
        hasParent = False
        userCov = userCoverageDict[kmer]
        if userCov <= len(userSeqLis) * covPrecThres: continue # do filtering when coverage is too low

        for lastLayerKmer in lastLayerKmerLis:
            pos = kmer.find(lastLayerKmer)  # only two possible pos: 0 or 1
            if pos != -1:
                filtSigKmerSet.add(kmer)
                parent = lastLayerKmer if K > 2 else ''  # the first layer has no parents
                kmer2parentDic[kmer] = parent
                if parent:
                    curLayerP2Kmers.setdefault(parent, [set(), set()])  # two possible pos: 0 or 1
                    curLayerP2Kmers[parent][pos].add(kmer)
                    hasParent = True

        if not hasParent:
            parent = '_' * (K - 1)
            curLayerP2Kmers.setdefault(parent, [set(), set()])  # two possible pos: 0 or 1
            curLayerP2Kmers[parent][0].add(kmer)

    curLayer = {'root': [filtSigKmerSet, []]} if K == 2 else curLayerP2Kmers
    wholeTree.append(curLayer)

    return filtSigKmerSet, wholeTree
#end_func

def DisplayLayerInfo(title, origKmers, filteredKmers, wholeTree, kmer2Pvalue):
    """
    display information for each K
    """
    INFO(title)
    INFO("original significant kmer: %s" % origKmers)
    INFO("filtered significant kmer: %s" % filteredKmers)
    for kmer in filteredKmers:
        if kmer not in kmer2Pvalue: continue
        INFO('kmer: %s, pvalue: %s' % (kmer, kmer2Pvalue[kmer]))
    # print "tree structure in current layer %s" % wholeTree[-1]
#end_func


def DrawTree(outputFn, wholeTree, overallKmer2Cov, overallRefKmer2CovThreshold, enableMiRSampleing=True, enableView=False):
    dot = Digraph(format='png')
    node2idx = {}
    K2kmers = {}
    whiteEdgeLis = []
    normalEdgeLis = []

    # add root node and 2-mer nodes
    rootIdx = '1'
    node2idx['root'] = rootIdx
    dot.node(rootIdx, 'root', color='red')
    K2kmers[1] = (set('root'), set('root'))

    # add nodes in the deeper layer and edges
    RNATree = wholeTree[0]
    miRTree = wholeTree[1] if enableMiRSampleing else []

    RNADepth = len(RNATree)
    miRDepth = len(miRTree)
    maxDepth = max(RNADepth, miRDepth)

    for layerId in xrange(maxDepth):
        K = layerId + 2
        RNARefKmer2CovThreshold = overallRefKmer2CovThreshold[0]
        miRRefKmer2CovThreshold = overallRefKmer2CovThreshold[1] if enableMiRSampleing else {}
        RNALayerDic = RNATree[layerId] if len(RNATree) > layerId else {}
        miRLayerDic = miRTree[layerId] if len(miRTree) > layerId else {}

        RNAParent = set(RNALayerDic.keys())
        miRParent = set(miRLayerDic.keys())
        allParent = RNAParent | miRParent

        # merge nodes in left and right into one set
        RNAtmp = map(lambda x: set(x[0]) | set(x[1]), RNALayerDic.values())
        layerAllRNAKmers = reduce(lambda x, y: x | y, RNAtmp) if RNAtmp else set()
        miRtmp = map(lambda x:set(x[0]) | set(x[1]), miRLayerDic.values())
        layerAllmiRKmers = reduce(lambda x,y:x | y, miRtmp) if miRtmp else set()
        layerCommKmers = layerAllRNAKmers.intersection(layerAllmiRKmers)
        K2kmers[K] = (layerAllRNAKmers - {'_' * K}, layerAllmiRKmers - {'_' * K})

        for parentNode in allParent: # go through all the parent node
            parentNodeIdx = node2idx[parentNode]
            RNALeftKmers, RNARightKmers = RNALayerDic.get(parentNode, (set(), set()))
            miRLeftKmers, miRRightKmers = miRLayerDic.get(parentNode, (set(), set()))
            allRNAKmers = set(RNALeftKmers) | set(RNARightKmers)
            allMiRKmers = set(miRLeftKmers) | set(miRRightKmers)
            allKmers = allRNAKmers | allMiRKmers

            for kmer in allKmers: # go through all the children node
                node2idx.setdefault(kmer, str(len(node2idx) + 1))
                childNodeIdx = node2idx[kmer]
                if kmer != '_' * K and parentNode != '_' * (K - 1): # normal node
                    RNAThreshold = RNARefKmer2CovThreshold.get(kmer, 'N/A')
                    miRThreshold = miRRefKmer2CovThreshold.get(kmer, 'N/A')
                    nodeColor = 'red' if kmer in layerCommKmers else 'black' if kmer in allRNAKmers else 'green'
                    nodeLabel = '%s_%s,%s,%s' % (kmer, overallKmer2Cov[K][kmer], RNAThreshold, miRThreshold)
                    nodeStyle = 'solid'
                    edgeColor = nodeColor
                    edgestyle = 'solid'
                    edgeShape = 'normal'
                elif kmer == '_' * K: # hidden node
                    nodeColor = 'white'
                    nodeLabel = ""
                    nodeStyle = 'solid'
                    edgeColor = 'white'
                    edgestyle = 'dotted'
                    edgeShape = 'none'
                else: # separate node
                    RNAThreshold = RNARefKmer2CovThreshold.get(kmer, 'N/A')
                    miRThreshold = miRRefKmer2CovThreshold.get(kmer, 'N/A')
                    nodeColor = 'red' if kmer in layerCommKmers else 'black' if kmer in allRNAKmers else 'green'
                    nodeLabel = '%s_%s,%s,%s' % (kmer, overallKmer2Cov[K][kmer], RNAThreshold, miRThreshold)
                    nodeStyle = 'diagonals'
                    edgeColor = 'white'
                    edgestyle = 'dotted'
                    edgeShape = 'none'
                dot.node(childNodeIdx, nodeLabel, color=nodeColor, style=nodeStyle)
                dot.edge(parentNodeIdx, childNodeIdx, color=edgeColor, style=edgestyle, shape=edgeShape)
    try:
        dot.render(outputFn, view=enableView)
    except:
        pass
    return K2kmers
#end_func

def DrawCovTrend(outputDir, K2kmers, overallKmer2VisDetail):
    """
    draw the coverage trend for each layer 
    """
    allK2cov = []
    K = 0
    for K, (RNAKmers, miRKmers) in K2kmers.iteritems():
        if K == 1: continue
        allKmers = RNAKmers | miRKmers
        tmp = set()
        for kmer in allKmers:
            visSeq = set(map(lambda x:x[0], overallKmer2VisDetail[K][kmer]))
            tmp = tmp | visSeq
        allK2cov.append((K, len(tmp)))
    allK2cov.append((K + 1, 0))

    # draw layer coverage trend
    pngFn = os.path.join(outputDir, 'layer coverage.png')
    # DrawLayerCov(map(lambda x: x[1], allK2cov), 'coverage in each layer', pngFn)

    # draw layer coverage change trend
    layerCovChangeLis = GetCovChangeTrend([map(lambda x: x[1], allK2cov)])[0]
    pngFn = os.path.join(outputDir, 'layer coverage change.png')
    # DrawLayerCov(layerCovChangeLis, 'coverage change in each layer', pngFn)
    print9(isLogging=True)

    motifLength = layerCovChangeLis.index(max(layerCovChangeLis)) + 2

    return motifLength
#end_func

# def DrawLayerCov(layerCovLis, title, outputFn):
#     plt.clf()
#     treeDepth = len(layerCovLis)
#     xData = map(lambda x:x + 2, range(treeDepth))
#     yData = layerCovLis
#     plt.plot(xData, yData)
#     plt.title(title)
#     plt.savefig(outputFn, format='png')
# #end_func

def GetCovChangeTrend(covTrendLis):
    covChangeTrendLis = []
    for covTrend in covTrendLis:
        covChangeTrend = [covTrend[idx] - covTrend[idx + 1] for idx in xrange(len(covTrend) - 1)]
        covChangeTrendLis.append(covChangeTrend)
    return covChangeTrendLis
#end_func
