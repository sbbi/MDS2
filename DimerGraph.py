#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/5/14
Description     :

"""
import BioinfoComm
import CommStructure
import itertools
import kmerGraphStructure

def DivideSeq(seqLis, k=2):
    """
    divide sequences into kmers
    default k is 2
    :return:
    kmers matrix in each sequences
    """
    kmerMatrix = []
    for seq in seqLis:
        kmerLis = BioinfoComm.Seq2Kmer(seq, k)
        kmerMatrix.append(kmerLis)
    return kmerMatrix
# end_func

def BuildGraph(kmerMatrix, k=2, alphabet='ACTG'):
    # build a empty kmerGraph(all node and no edge)
    kmerGraph = kmerGraphStructure.KmerGrpah(True)
    for item in itertools.product(alphabet, repeat=k):
        kmer = ''.join(item)
        curNode = CommStructure.Node(kmer, [])
        kmerGraph.addNode(curNode)
    edgeTitle = None
    edgeWeight = 1
    edge = None

    # add edges(real data in sequences)
    for seqIdx, kmerLis in enumerate(kmerMatrix):
        for pos, kmer in enumerate(kmerLis[:-1]):
            startNode = CommStructure.Node(kmer)
            endNode = CommStructure.Node(kmerLis[pos + 1])
            edge = CommStructure.Edge(startNode, endNode, edgeTitle, edgeWeight, (seqIdx, pos))
            kmerGraph.addEdge(edge)
        # update the last kmer in a sequence
        kmerGraph.updateLastNode(edge)
        # sequence count + 1
        kmerGraph.addSeq()
    #end_for

    # update statistical info of the graph
    kmerGraph.updateStatDict()

    return kmerGraph
# end_func

def SearchPaths(dimerGraph, maxPathLength=20, enableFiltering=False, covPercThres=0.1):
    """
    Search all the possible paths of 2-mers
    :param covPrecThres:
    :param dimerGraph:
    :return:
        all the paths in the graph
    """
    dimerGraph.findAllPath(maxPathLength, enableFiltering, covPercThres)
    return dimerGraph.path
#end_func

def PathDic2PathCov(pathDic, pathLength):
    """
    format path info into two dict: kmer2coverage and kmer2VisitDetail
    kmer2VisitDetail includes the sequences that the kmer cover
    """
    kmer2coverage = {}
    kmer2VisitDetail = {}
    pathLis = pathDic[pathLength - 1]
    for path in pathLis:
        kmer = BioinfoComm.Kmer2Str(path.nodeLis)
        coverage = path.coverage
        visitDetail = path.visitLis
        kmer2coverage[kmer] = coverage
        kmer2VisitDetail[kmer] = visitDetail
    return kmer2coverage, kmer2VisitDetail
#end_func