#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
ProjectName     : ''
Author          : 'Tian Gao'
CreationDate    : '2017/5/14'
Description     :

"""

import Comm
import BioinfoComm
import graphviz
from joblib import Parallel, delayed
import logging

INFO = logging.getLogger('root').info

class Graph(object):
    def __init__(self,*args,**kwargs):
        self.node_neighbors = {}
        self.visited = {}
        self.allPattern = set()

    def add_nodes(self,nodelist):
        for node in nodelist:
            self.add_node(node)
    #end_func

    def add_node(self,node):
        if not node in self.nodes():
            self.node_neighbors[node] = set()
    #end_func

    def add_edge(self,edge):
        u,v = edge
        if(v not in self.node_neighbors[u]) and ( u not in self.node_neighbors[v]):
            self.node_neighbors[u].add(v)

            if(u!=v):
                self.node_neighbors[v].add(u)
    #end_func

    def nodes(self):
        return self.node_neighbors.keys()
    #end_func

    def FetchAllPath(self, stack):
        self.allPattern.add(BioinfoComm.SegmentMerge(sorted(stack)))
        lastNode = stack[-1]
        for node in self.node_neighbors[lastNode]:
            if node not in stack:
                stack.append(node)
                self.FetchAllPath(stack)
                stack.pop()
    #end_func
#end_class

def BuildSimiGraph(graphFn, segmentLis, enableView=False):
    g = Graph()
    g.add_nodes(segmentLis)
    dot = graphviz.Graph(format='png')
    nodeColor = 'white'
    for fragIdx1, fragment1 in enumerate(segmentLis):
        dot.node(fragment1, style="filled", fillcolor=nodeColor)
        for fragment2 in segmentLis:
            if fragment1 > fragment2 and BioinfoComm.IsSimi(fragment1, fragment2):
                g.add_edge((fragment1, fragment2))
                dot.edge(fragment1, fragment2)
    try:
        dot.render(graphFn, view=enableView)
    except:
        pass

    return g
#end_func

def ParallelDiffStartnode(startNode, bp, adjMap, startNode2rlt, minIC, maxIC, overallKmer2Cov):
    includedNodes = {startNode}
    availNodes = adjMap[startNode]
    disableNode = set()
    userCov, curPattern = startNode2rlt[startNode] if startNode in startNode2rlt and minIC <= startNode2rlt[startNode][
        0] <= maxIC else bp.searchPattern(includedNodes, availNodes, disableNode)
    # print 'current starting node: %s, current pattern: %s' % (startNode, curPattern)
    if not curPattern: return startNode, -1, '', 0, {}

    # calculate new information content
    curIC, bp.pattern2IC = FetchPatternWeighedIC(curPattern, bp.pattern2IC, overallKmer2Cov)
    return startNode, curIC, curPattern, userCov, bp.pattern2IC
#end_func

def FetchPatternWithICIter(allNode, adjMap, minIC, maxIC, seqLis, allPattern, kmer2seqIdSet, overallKmer2Cov, pattern2IC, startNode2rlt={}, ICminInterval=0.05, enableParallel=False, coreNm=8, pattern2kmerSet={}):
    Comm.PrintWithTime('current IC: min: %s, max: %s' % (minIC, maxIC), isLogging=True)

    # build the backpack model
    bp = Backpack(adjMap, minIC, maxIC, seqLis, kmer2seqIdSet, overallKmer2Cov, pattern2IC, pattern2kmerSet)

    # give different start node in case that there are multiple motif
    patternSet = set()
    ICLis = []

    if not enableParallel:
        for startNode in allNode:
            includedNodes = {startNode}
            availNodes = adjMap[startNode]
            disableNode = set()
            userCov, curPattern = startNode2rlt[startNode] if startNode in startNode2rlt and minIC <= startNode2rlt[startNode][0] <= maxIC else bp.searchPattern(includedNodes, availNodes, disableNode)
            startNode2rlt[startNode] = (userCov, curPattern)
            # INFO('current starting node: %s, current pattern: %s' % (startNode, curPattern))
            if not curPattern: continue

            # filter 'ACGU' and 'X[ACGU]X'
            parsedPatternLis = BioinfoComm.parsePatternStr(curPattern)
            charCntLis = map(lambda x: len(x), parsedPatternLis)
            if set(charCntLis) == {1} or 4 in charCntLis: continue

            # ==== TODO
            # kmerLis = BioinfoComm.FetchAllKmerFromPattern(curPattern)
            # patternLen = len(kmerLis[0])
            # if patternLen == 4:
            # if patternLen == 6:
            #     curPattern = '[GT]G[ACG][CG]'
            #     curPattern = '[CGT]CC[AGT]'
                # curPattern = '[AGT][AGT][AGT][ACG][AT][CG]'
            # ====

            # calculate new information content
            curIC, bp.pattern2IC = FetchPatternWeighedIC(curPattern, bp.pattern2IC, overallKmer2Cov)
            ICLis.append(curIC)
            patternSet.add(curPattern)
    else:
        startNodeInfoLis = Parallel(n_jobs=coreNm)(delayed(ParallelDiffStartnode)(startNode, bp, adjMap, startNode2rlt, minIC, maxIC, overallKmer2Cov) for startNode in allNode)
        for startNode, curIC, curPattern, userCov, curPattern2IC in startNodeInfoLis:
            if curIC == -1: continue
            pattern2IC = dict(pattern2IC, **curPattern2IC)
            startNode2rlt[startNode] = (userCov, curPattern)
            ICLis.append(curIC)
            patternSet.add(curPattern)

    kmer2seqIdSet = bp.kmer2seqIdSet
    pattern2kmerSet = bp.pattern2kmerSet
    # pattern2IC = bp.pattern2IC
    if not patternSet:return allPattern, kmer2seqIdSet, pattern2IC, pattern2kmerSet # in case there is no result
    allPattern = allPattern | patternSet
    nextIC = min(ICLis)

    # the left part
    minIC1 = minIC
    maxIC1 = min(nextIC, maxIC - 0.05) # at least move the bound for 0.05
    if maxIC1 < maxIC and maxIC1 - minIC1 > ICminInterval:
        allPattern, kmer2seqIdSet, pattern2IC, pattern2kmerSet = FetchPatternWithICIter(allNode, adjMap, minIC1, maxIC1, seqLis, allPattern, kmer2seqIdSet, overallKmer2Cov, pattern2IC, startNode2rlt, enableParallel=enableParallel, coreNm=coreNm, pattern2kmerSet=pattern2kmerSet)

    # the right part
    minIC2 = max(nextIC, minIC + 0.05) # at least move the bound for 0.05
    maxIC2 = maxIC
    if minIC2 > minIC and maxIC2 - minIC2 > ICminInterval:
        allPattern, kmer2seqIdSet, pattern2IC, pattern2kmerSet = FetchPatternWithICIter(allNode, adjMap, minIC2, maxIC2, seqLis, allPattern, kmer2seqIdSet, overallKmer2Cov, pattern2IC, startNode2rlt, enableParallel=enableParallel, coreNm=coreNm, pattern2kmerSet=pattern2kmerSet)

    INFO('iteration for IC [%s, %s]' % (minIC, maxIC))

    return allPattern, kmer2seqIdSet, pattern2IC, pattern2kmerSet
#end_func

class Backpack:
    def __init__(self, adjMap, minIC, maxIC, seqLis, kmer2seqIdSet, overallKmer2Cov, pattern2IC, pattern2kmerSet={}):
        self.adjMap = adjMap
        self.minIC = minIC
        self.maxIC = maxIC
        self.seqLis = seqLis
        self.nodesStr2rlt = {} # store calculated mediate result in BP problem
        self.kmer2seqIdSet = kmer2seqIdSet
        self.overallKmer2Cov = overallKmer2Cov
        self.pattern2IC = pattern2IC
        self.pattern2kmerSet = pattern2kmerSet
    # end_init

    def searchPattern(self, includedNodes, availNodes, disableNode):
        """
        backpack problem
        information content is considered as size of an item in pb.
        coverage is considered as value of an item in pb.
        """

        includedNodesStr = ''.join(includedNodes)
        disableNodesStr = ''.join(disableNode)
        if (includedNodesStr, disableNodesStr) in self.nodesStr2rlt:
            return self.nodesStr2rlt[(includedNodesStr, disableNodesStr)]

        if not availNodes:# basic case 1
            if not includedNodes:
                self.nodesStr2rlt[(includedNodesStr, disableNodesStr)] = (0, '')
                return 0, ''
            curPattern = BioinfoComm.SegmentMerge(list(includedNodes))
            curCov, curSeq, self.kmer2seqIdSet = BioinfoComm.FetchPatternCov(curPattern, self.seqLis, self.kmer2seqIdSet)

            self.nodesStr2rlt[(includedNodesStr, disableNodesStr)] = (curCov, curPattern)
            self.pattern2kmerSet[curPattern] = includedNodes
            return curCov, curPattern

        if len(availNodes) == 1: # basic case 2
            nextNodesSet = includedNodes | availNodes
            nextPattern = BioinfoComm.SegmentMerge(list(nextNodesSet))
            nextIC, self.pattern2IC = FetchPatternWeighedIC(nextPattern, self.pattern2IC, self.overallKmer2Cov)

            # if the last node is affordable, just take it.
            nextNodesSet = nextNodesSet if self.maxIC > nextIC > self.minIC else includedNodes
            if not nextNodesSet:
                self.nodesStr2rlt[(includedNodesStr, disableNodesStr)] = (0, '')
                return 0, ''
            nextPattern = BioinfoComm.SegmentMerge(list(nextNodesSet))
            nextCov, nextSeq, self.kmer2seqIdSet = BioinfoComm.FetchPatternCov(nextPattern, self.seqLis, self.kmer2seqIdSet)
            nextNodesStr = ''.join(nextNodesSet)
            self.nodesStr2rlt[(nextNodesStr, disableNodesStr)] = (nextCov, nextPattern)
            self.pattern2kmerSet[nextPattern] = nextNodesSet
            return nextCov, nextPattern

        # TODO: may need change. original: for node in availNodes
        node = availNodes.pop()

        # the case that current node is NOT selected
        nextWithoutIncludedNodes = includedNodes
        nextWithoutDisableNode = disableNode | {node}
        nextWithoutDisableNodeStr = ''.join(nextWithoutDisableNode)
        nextWithoutAvailNodes = availNodes - nextWithoutIncludedNodes - nextWithoutDisableNode
        withoutNodeCov, withoutNodePattern = self.searchPattern(nextWithoutIncludedNodes, nextWithoutAvailNodes, nextWithoutDisableNode)

        # the case that current node is selected
        nextWithNodesSet = includedNodes | {node}
        nextWithPattern = BioinfoComm.SegmentMerge(list(nextWithNodesSet))
        nextWithIC, self.pattern2IC = FetchPatternWeighedIC(nextWithPattern, self.pattern2IC, self.overallKmer2Cov)

        if not self.maxIC > nextWithIC > self.minIC:
            # if the new node is not affordable, return the result without the current node
            nextWithoutIncludedNodesStr = ''.join(nextWithoutIncludedNodes)
            self.nodesStr2rlt[(nextWithoutIncludedNodesStr, nextWithoutDisableNodeStr)] = (withoutNodeCov, withoutNodePattern)
            return withoutNodeCov, withoutNodePattern
        else:
            # if the new node is affordable, compare the result with current node and the result without current node, then return the best one
            nextWithIncludedNodes = includedNodes | {node}
            nextWithDisableNode = disableNode
            nextWithDisableNodeStr = disableNodesStr
            # nextWithAvailNodes = (availNodes | self.adjMap[node]) - (nextWithIncludedNodes | nextWithDisableNode) # this one is slower
            nextWithAvailNodes = reduce(lambda seqSet1, seqSet2: seqSet1 | seqSet2, map(lambda node: self.adjMap[node], nextWithIncludedNodes)) - nextWithIncludedNodes - nextWithDisableNode
            withNodeCov, withNodePattern = self.searchPattern(nextWithIncludedNodes, nextWithAvailNodes, nextWithDisableNode)
            withNodeCov, withNodeSeq, self.kmer2seqIdSet = BioinfoComm.FetchPatternCov(withNodePattern, self.seqLis, self.kmer2seqIdSet)
            withoutNodeCov, withoutNodeSeq, self.kmer2seqIdSet = BioinfoComm.FetchPatternCov(withoutNodePattern, self.seqLis, self.kmer2seqIdSet)

            if withNodeCov > withoutNodeCov: # compare the two cases, return the better one
                nextWithIncludedNodesStr = ''.join(nextWithIncludedNodes)
                self.nodesStr2rlt[(nextWithIncludedNodesStr, nextWithDisableNodeStr)] = (withNodeCov, withNodePattern)
                return withNodeCov, withNodePattern
            else:
                nextWithoutIncludedNodesStr = ''.join(nextWithoutIncludedNodes)
                self.nodesStr2rlt[(nextWithoutIncludedNodesStr, nextWithoutDisableNodeStr)] = (withoutNodeCov, withoutNodePattern)
                return withoutNodeCov, withoutNodePattern
    # end_func
#end_class

def FetchPatternWeighedIC(pattern, pattern2IC, overallKmer2Cov):
    """
    get information content of a pattern weighted by coverage of each kmer
    'A[CT]G', ACG:2, ATG:3 ->
    information content of ['AAAAA', 'CCTTT', 'GGGGG']
    """
    if pattern in pattern2IC:
        IC = pattern2IC[pattern]
    else:
        ICLis = []
        kmerId2cov = {}
        kmerLis = BioinfoComm.FetchAllKmerFromPattern(pattern)
        stringLen = len(kmerLis[0])
        stringCnt = len(kmerLis)
        for kmerId, kmer in enumerate(kmerLis):
            cov = overallKmer2Cov[stringLen].get(kmer, 1)
            kmerId2cov[kmerId] = cov
        for posId in xrange(stringLen):
            stringInEachPos = ''
            for kmerId in xrange(stringCnt):
                weight = kmerId2cov[kmerId]
                letter = kmerLis[kmerId][posId]
                stringInEachPos += letter * weight
            ICInCurPos = BioinfoComm.GetSingleInfoContentByString(stringInEachPos)
            ICLis.append(ICInCurPos)
        IC = sum(ICLis)
        pattern2IC[pattern] = IC
    return IC, pattern2IC
#end_func

def test():
    # segmentLis = set(['UUUCC', 'UUUCA', 'CACUA', 'UUUGU', 'GCACU', 'UUUCU', 'UAUUU', 'CAAAA', 'AUAAA', 'UUAAU', 'UAAUA', 'UUGUU', 'GUUUU', 'UUAAA', 'UUAAG', 'UACAC', 'UACAA', 'UAAUU', 'UUGUA', 'AUUCU', 'AAAUA', 'AGAUU', 'AAAUU', 'UUUUA', 'CAAAU', 'CAGCA', 'UUCUG', 'AACAC', 'AACAA', 'UUUUG', 'UUCUU', 'GUGAA', 'AACAU', 'ACUAC', 'GUGCA', 'UGCAG', 'AUACU', 'CUUCU', 'ACAUU', 'UGCAC', 'CACUU', 'UUUUU', 'ACUUU', 'CUCAA', 'CUGAU', 'UUUUC', 'UUACA', 'AAUUA', 'AAUUC', 'CUUUC', 'CUUUA', 'AUGUA', 'UCAUU', 'AAAAA', 'CUUUU', 'AAUUU', 'UAUGC', 'CUACU', 'AUGUU', 'AUUUU', 'AAAGU', 'CAUUA', 'AUUUG', 'CUCUC', 'ACUGU', 'CAUUG', 'UUGCA', 'AAGUG', 'CUUAA', 'UCACU', 'ACAAU', 'AUCAU', 'CACAU', 'AGUGC', 'UAGUG', 'CUGUG', 'UUCAC', 'AUUAU', 'CUACA', 'UUCAU', 'AGUGU', 'CACAA', 'ACACU', 'AUUGU', 'UGUAG', 'UGUAA', 'UUAUU', 'AAAAU', 'UGUAU', 'AAUGG', 'AAAAC', 'UUUAG', 'UUUAA', 'AAACA', 'AAACU', 'UCAAA', 'UGUUA', 'AGGAA', 'UCUUA', 'UCUUU', 'UGUUU', 'CUAAA', 'UAAAU', 'ACAAA', 'AUUUC', 'CAUUU', 'AUUAA', 'UGCAU', 'AAUAU', 'UAAAG', 'UAAAA', 'AUAUG', 'GUAGU', 'AAUAA', 'UAUCA', 'UACUG', 'GUUAA'])
    # segmentLis = set(['ACAG', 'AAAU', 'ACAA', 'UAUA', 'UAUC', 'UCAU', 'AAAA', 'AAAC', 'GUUU', 'AAAG', 'UCAC', 'GUGU', 'UCAA', 'UAUU', 'UAAU', 'UGUU', 'AAUU', 'GUGC', 'ACUU', 'UGUA', 'UGUG', 'UAAA','ACUA', 'ACUG', 'CUUG', 'UUUG', 'CUUC', 'CAUU', 'AUCA', 'AUGG', 'GCAC', 'UUUU', 'CUUU', 'UAUG', 'GCAU', 'CAAU', 'GUUA', 'AUGU', 'UCUG', 'UCUC', 'UCUA', 'AAUA', 'UGAU', 'UUAU', 'ACAU', 'CUAU', 'UCUU','UUAA', 'AUUU', 'CUAA', 'GAUU', 'UUAG', 'AUUA', 'AUUC', 'UGCA', 'AUUG', 'UGGA', 'AACA', 'UGCU', 'AAGU', 'CUCU', 'UUUA', 'UGGU', 'AACU', 'UUUC', 'UAGU', 'UCAG', 'UACA', 'ACAC', 'UACU', 'AAUG', 'CUUA','CACU', 'AGUG', 'AUAA', 'CAGU', 'AUAC', 'CUAC', 'CACA', 'CAAA', 'AUAU', 'AGUU', 'UUCU', 'CUGU', 'UUGC', 'AUCU', 'UUCA', 'UUGU', 'CUGC'])
    # print len(segmentLis)
    # BuildSimiGraph(r'C:\Users\Administrator\Desktop\2\1', segmentLis, True)
    pattern = '[CTG]'
    overallKmer2Cov = {}
    overallKmer2Cov[1] = {'C':10, 'T': 10, 'G': 10}
    print 1
    print FetchPatternWeighedIC(pattern, {}, overallKmer2Cov)
# end_test

def main():
    test()
# end_main

if __name__ == "__main__":
    main()
# end_if