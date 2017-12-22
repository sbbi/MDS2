#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     = 'MotifFinding'
Author          = 'Tian Gao'
CreationDate    = '2016/11/11'
Description     =
This file includes all the structure and related functions to build a k-mer graph
"""
import sys
import Comm
import CommStructure
import BioinfoComm
import logging

INFO = logging.getLogger('root').info

SEEDREGION = range(2, 8)
MINPATTERNCOVERAGERATE = 0.1

class Path():
    """
    path is a list of kmer that stored in the nodeLis.
    visitLis gives info about sequences that this path covers.
    """
    def __init__(self, nodeLis, visitLis):
        self.startNode = nodeLis[0]
        self.endNode = nodeLis[-1]
        self.nodeLis = nodeLis
        self.pathLength = len(nodeLis) - 1
        self.visitLis = visitLis
        self.coverage = len(set(seqId for seqId, pos in visitLis))

    # end_init

    @staticmethod
    def extendPath(path1, path2):
        """
        intersect path1 with path2
        ('AT', 'TG'), ('TG', 'GC') -> ('AT', 'TG', 'GC')

        note:
        endNode of path1 must be the same with startNode of path2

        :return:
        new path
        """
        newNodeLis = path1.nodeLis[:-1] + path2.nodeLis
        newVisitLis = []
        # TODO: this loop can be improved
        for seqId1, pos1 in path1.visitLis:
            for seqId2, pos2 in path2.visitLis:
                if seqId1 == seqId2 and pos2 == pos1 + path1.pathLength:
                    newVisitLis.append((seqId1, pos1))
        return Path(newNodeLis, newVisitLis)
        # end_func


# end_class

class Pattern():
    def __init__(self, patternStr, visitLis):
        self.patternStr = patternStr
        self.visitLis = visitLis
        self.patternLength = BioinfoComm.getPatternLength(patternStr)
        self.coverage = len(set(seqId for seqId, pos in visitLis))
        self.mergeCnt = BioinfoComm.getMergeCnt(patternStr)

    # end_init

    @staticmethod
    def mergePattern(pattern1, pattern2, maxPatternDist=1):
        """
        merge pattern1 with pattern2
        ATG, ACG -> A[TC]G

        note:
        the two paths must have equal length

        :return:
        new pattern
        """
        newStrLis = []
        patternDist = 0
        strLis1 = BioinfoComm.parsePatternStr(pattern1.patternStr)
        strLis2 = BioinfoComm.parsePatternStr(pattern2.patternStr)
        for strIdx, curStr1 in enumerate(strLis1):
            curStr2 = strLis2[strIdx]
            # calculate distance and generate a new string
            if curStr1 == curStr2:  # (A, A) -> A, (TC, TC) -> TC
                newStr = curStr1
            elif curStr1 in curStr2: # [CT] and C -> don't merge if occur twice
                patternDist += 0.6
                newStr = curStr2
            elif curStr2 in curStr1:
                patternDist += 0.6
                newStr = curStr1
            else: # (A, C) -> AC, (TC, TA) -> TCA, (CA, TG) -> X
                # TODO: the following command is slow, and both of them have bugs : [AC], [TG] / [ACG], [AG]
                # curStrSet1 = set(curStr1)
                # curStrSet2 = set(curStr2)
                # patternDist += min(len(curStrSet1.difference(curStrSet2)), len(curStrSet2.difference(curStrSet1)))
                patternDist += 1
                curLis = list(set(curStr1).union(set(curStr2)))
                newStr = ''.join(sorted(curLis))
            if patternDist > maxPatternDist or len(newStr) == 4: return None
            newStr = '[%s]' % newStr if len(newStr) >= 2 else newStr  # AT -> [AT]
            newStrLis.append(newStr)
        if patternDist == 0:return None
        newPatternStr = ''.join(newStrLis)
        visitLis = list(set(pattern1.visitLis + pattern2.visitLis))
        newPattern = Pattern(newPatternStr, visitLis)
        return newPattern
    # end_func

    def AnalyzePatternStr(self):
        """
        analyze pattern string and get three features
        fixedCnt: how many fixed letters in the string
        maxContinuedFixedLetterCnt: how long is the longest continued fixed letters
        varCnt: how many possible different path the string can have
        maxContinuedVarLetterCnt: how many continued positions have variables
        varPosRate: what is the percentage of positions that have variables
        """
        fixedStr, _, tmp = self.patternStr.partition('[')
        singlePattern, _, rest = tmp.partition(']')
        maxContinuedFixedLetterCnt = len(fixedStr)
        maxContinuedVarLetterCnt = 1
        curContinuedVarLetterCnt = 1
        fixedCnt = len(fixedStr)
        varCnt = max(1, len(singlePattern))
        varPosCnt = self.patternStr.count('[')
        varPosRate = varPosCnt * 1.0 / self.patternLength if self.patternLength > 0 else 0.1
        while rest:
            fixedStr, _, tmp = rest.partition('[')
            singlePattern, _, rest = tmp.partition(']')
            curContCnt = len(fixedStr)
            if curContCnt > maxContinuedFixedLetterCnt: maxContinuedFixedLetterCnt = curContCnt
            fixedCnt += len(fixedStr)
            varCnt *= max(1, len(singlePattern))
            if not fixedStr:
                curContinuedVarLetterCnt += 1
            else:
                maxContinuedVarLetterCnt = max(maxContinuedVarLetterCnt, curContinuedVarLetterCnt)
        maxContinuedVarLetterCnt = max(maxContinuedVarLetterCnt, curContinuedVarLetterCnt)
        return self.coverage, fixedCnt, maxContinuedFixedLetterCnt, varCnt, maxContinuedVarLetterCnt, varPosRate, self.patternLength
    #end_func

# end_class

class KmerGrpah(CommStructure.Graph):
    def __init__(self, isDirected=False):
        self.isDirected = isDirected
        self.nodeCnt = 0
        self.edgeCnt = 0
        self.seqCnt = 0
        self.seqLen = {}
        self.nodeDict = {}
        self.edgeDict = {}
        self.adjacentEdgeDict = {}
        self.adjacentNodeDict = {}

        self.path = {}
        self.pattern = {}

        # statistical dict
        self.coverageDict = {}
        self.visitFreqDict = {}
        # statistical info about edges between two nodes
        # (node1, node2) -> [pairVisitCnt, pairSeedVisitCnt, pairNonseedVisitCnt,\
        #                    pairVisitCoverage, pairSeedVisitCoverage, pairNonseedVisitCoverage]
        # TODO 1: position bias of motif
        # self.pairNodeStat = {}

    # end_init

    def addSeq(self):
        self.seqCnt += 1

    # end_func

    def getNextNodeLis(self, currentNodeTitle):
        nodeLis = []
        for adjacentNode, edgeLis in self.adjacentNodeDict[currentNodeTitle].iteritems():
            S2EEdgeLis = [isS2E for (seqId, pos, isS2E) in edgeLis if isS2E]
            if S2EEdgeLis: nodeLis.append(adjacentNode)
        return nodeLis

    # end

    def updateStatDict(self):
        node2CoveredSeq = {}
        for node1Title, adjacentNodeDict in self.adjacentNodeDict.iteritems():
            node2CoveredSeq.setdefault(node1Title, [])
            for node2Title, edgeLis in adjacentNodeDict.iteritems():
                # pair node statistics
                visitLis = [seqId for (seqId, pos, isNode1ToNode2) in edgeLis if isNode1ToNode2]
                seedVisitLis = [seqId for (seqId, pos, isNode1ToNode2) in edgeLis if
                    pos in SEEDREGION and isNode1ToNode2]
                nonseedVisitLis = [seqId for (seqId, pos, isNode1ToNode2) in edgeLis if
                    pos not in SEEDREGION and isNode1ToNode2]
                pairVisitCnt = len(visitLis)
                pairSeedVisitCnt = len(seedVisitLis)
                pairNonseedVisitCnt = len(nonseedVisitLis)
                pairVisitCoverage = len(set(visitLis))
                pairSeedVisitCoverage = len(set(seedVisitLis))
                pairNonseedVisitCoverage = len(set(nonseedVisitLis))

                # TODO 1: position bias of motif
                # self.pairNodeStat[(node1Title, node2Title)] = (
                #     pairVisitCnt, pairSeedVisitCnt, pairNonseedVisitCnt, pairVisitCoverage, pairSeedVisitCoverage,
                #     pairNonseedVisitCoverage)

                # node statistics
                self.visitFreqDict[node1Title] = self.visitFreqDict.get(node1Title, 0) + len(visitLis)
                node2CoveredSeq.setdefault(node2Title, [])
                node2CoveredSeq[node1Title] += visitLis
                node2CoveredSeq[node2Title] += visitLis

                # construct init path(with two nodes)
                pathNodeLis = [node1Title, node2Title]
                visitLis = [(seqId, pos) for seqId, pos, isNode1ToNode2 in edgeLis if isNode1ToNode2]
                if not visitLis: continue
                curPath = Path(pathNodeLis, visitLis)
                self.path.setdefault(2, [])
                self.path[2].append(curPath)

        for nodeTitle, visitLis in node2CoveredSeq.iteritems():
            self.coverageDict[nodeTitle] = len(set(visitLis))

    # end_def

    def addEdge(self, edge):
        seqIdx, pos = edge.edgeData

        # if edge title is not given, the default value is the total number of edge
        # edgeData : [(seqIdx, pos)]
        edge.edgeTitle = str(self.edgeCnt) if edge.edgeTitle is None else edge.edgeTitle
        self.edgeDict.setdefault(edge.edgeTitle, edge)
        self.edgeCnt += 1

        # nodeData : [(seqIdx, pos)]
        startNode = self.nodeDict[edge.startNodeTitle]
        startNode.nodeData.append(edge.edgeData)

        '''
        update adjacentNodeDict:
        key: Given a node A
        value: title of nodes that have edges connected with A and where the edges locate in sequences
        '''
        self.adjacentNodeDict.setdefault(edge.startNodeTitle, {})
        self.adjacentNodeDict[edge.startNodeTitle].setdefault(edge.endNodeTitle, [])
        self.adjacentNodeDict.setdefault(edge.endNodeTitle, {})
        self.adjacentNodeDict[edge.endNodeTitle].setdefault(edge.startNodeTitle, [])
        curStartNodeInfo = (seqIdx, pos, True)
        curEndNodeInfo = (seqIdx, pos, False)
        self.adjacentNodeDict[edge.startNodeTitle][edge.endNodeTitle].append(curStartNodeInfo)
        self.adjacentNodeDict[edge.endNodeTitle][edge.startNodeTitle].append(curEndNodeInfo)

        '''
        update adjacentEdgeDict:
        key: Given a node A
        value: id of edges connecting with A
        '''
        self.adjacentEdgeDict.setdefault(edge.startNodeTitle, [])
        self.adjacentEdgeDict.setdefault(edge.endNodeTitle, [])
        curOutEdgeInfo = (edge.edgeTitle, False) if self.isDirected else (edge.edgeTitle, None)
        curInEdgeInfo = (edge.edgeTitle, True) if self.isDirected else (edge.edgeTitle, None)
        self.adjacentEdgeDict[edge.startNodeTitle].append(curOutEdgeInfo)
        self.adjacentEdgeDict[edge.endNodeTitle].append(curInEdgeInfo)

        # update sequence length
        seqId = edge.edgeData[0]
        self.seqLen[seqId] = self.seqLen.get(seqId, 1) + 1

    # end_func

    def fetchNMaxCoverageNode(self, n):
        tmpLis = sorted(self.coverageDict.iteritems(), key=lambda x: x[1], reverse=True)
        return tmpLis[:n]

    # end_func

    def updateLastNode(self, edge):
        lastNodeTitle = edge.endNodeTitle
        lastNode = self.nodeDict[lastNodeTitle]
        lastNodeSeqId, lastNodePos = edge.edgeData
        lastNode.nodeData.append((lastNodeSeqId, lastNodePos + 1))
        self.visitFreqDict[lastNodeTitle] = self.visitFreqDict.get(lastNodeTitle, 0) + 1

    # end_func

    def outputEdge(self, outputTarget=sys.stdout, isFile=False):
        if isFile:
            outputTarget = open(outputTarget, 'w')
        for edgeTitle, edge in self.edgeDict.iteritems():
            outputTarget.write("%s\t%s\t%s\n" % (edge.startNodeTitle, edge.edgeData[0], edge.endNodeTitle))
        if isFile:
            outputTarget.close()

    # end_func

    def findAllPath(self, maxMotifLength, enableFiltering=False, covPercThres=0.1):
        """
        calculate all the path(length < maxMotifLength - 1) by extending iteratively from initial path(length == 2)
        add initial path(with two nodes) in each iteration
        path length : 2 -> 3 -> 4 ...
        """
        # pathLength = self.path.keys()[0] # path.keys() list all the length, path.keys()[0] == 2
        pathLength = 2 # review modified

        # add 2-mer as a path with length = 1
        self.path[1] = []
        for dimer, nodeData in self.nodeDict.iteritems():
            dimerPath = Path([dimer], nodeData.nodeData)
            self.path[1].append(dimerPath)

        while pathLength < maxMotifLength - 1:
            path1Lis = self.path.get(pathLength, None)
            path2Lis = self.path.get(2, None)
            pathLength += 1
            if not path1Lis: break
            extendedPathLis = []
            for path1 in path1Lis:
                nodeLis1 = path1.nodeLis
                for path2 in path2Lis:
                    nodeLis2 = path2.nodeLis
                    if nodeLis1[-path2.pathLength:] != nodeLis2[:path2.pathLength]: continue
                    extendedPath = Path.extendPath(path1, path2)
                    if not extendedPath.visitLis: continue
                    if enableFiltering and extendedPath.coverage < covPercThres * self.seqCnt: continue  # review modified
                    extendedPathLis.append(extendedPath)
            if not extendedPathLis: continue
            self.path[pathLength] = extendedPathLis
    # end_func

    def fetchMotifSeg(self, maxN=1):
        """
        fetch top maxN extended paths that cover most sequences
        """
        rltMatrix = []
        for k, pathLis in self.path.iteritems():
            bestPathLis = sorted(pathLis, key=lambda x: x.coverage, reverse=True)[:maxN]
            rltMatrix.append((k, bestPathLis))
        return rltMatrix
    # end_func

    def findAllMotifPattern(self, maxMergeCnt, hasFilter=False):
        """
        ABANDONED!!

        calculate all the Patterns(mergeCnt <= maxMergeCnt) by merging iteratively
            from initial pattern(mergeCnt = 0, i.e. path)

        The patterns are stored in self.pattern which is a dict
            key : (k, mergeCnt)
            value : [(newPattern, pattern1, pattern2)]
        :param hasFilter:
            whether the filter works
        :param maxMergeCnt:
            how many times do we merge
        """
        for k, pathLis in self.path.iteritems():
            patternLis2 = set([(BioinfoComm.Kmer2Str(path.nodeLis), tuple(path.visitLis)) for path in pathLis])
            for mergeCnt in range(maxMergeCnt):
                print "searching motif pattern: length = %s, mergeCnt = %s " % (k + 1, mergeCnt + 1)
                Comm.PrintTime()

                # fetch two pattern list(contain all the patterns with a certain length and mergeCnt) to merge
                # patternLis2 is a list of paths and never change
                if (k, mergeCnt) in self.pattern: # (after the first loop) patternLis1 to be merged is pattern
                    mergedPatternLis = map(lambda x:x[0], self.pattern[(k, mergeCnt)])
                    # 'tuple' conversion is necessary otherwise there will be an error because list is unhashable
                    patternLis1 = set([(x.patternStr, tuple(x.visitLis)) for x in mergedPatternLis])
                    patternLis1 = [(x[0], x[1]) for x in patternLis1]
                    patternLis2 = patternLis2.union(patternLis1)
                else: # (in the first loop) patternLis1 to be merged is paths
                    patternLis1 = patternLis2

                # merge two pattern lists
                for (patternStr1, visitLis1) in patternLis1:
                    pattern1 = Pattern(patternStr1, list(visitLis1))
                    if hasFilter and self.IsPatternFiltered(visitLis1): continue
                    for (patternStr2, visitLis2) in patternLis2:
                        # TODO: to be improved
                        if patternStr1 >= patternStr2: continue
                        if self.IsPatternFiltered(visitLis1): continue
                        visitLis2 = list(visitLis2)
                        pattern2 = Pattern(patternStr2, visitLis2)
                        newPattern = Pattern.mergePattern(pattern1, pattern2)
                        if not newPattern: continue
                        self.pattern.setdefault((k, mergeCnt + 1), [])
                        self.pattern[(k, mergeCnt + 1)].append((newPattern, pattern1, pattern2))
    # end_func

    def IsPatternFiltered(self, visitLis):
        """
        ABANDONED!!
        filter a pattern if its coverage is too low
        """
        coveredSeqCnt = len(set(map(lambda x:x[0], visitLis)))
        coverageRate = coveredSeqCnt * 1.0 / self.seqCnt
        rlt = True if coverageRate < MINPATTERNCOVERAGERATE else False
        return rlt
    # end_func

    def fetchMotifPatternInfo(self, maxN=1):
        """
        fetch motif pattern with top maxN coverage
        fetch a dict which maps pattern info to covered sequences
        """
        rltMatrix = []
        pattern2coverSeqLis = {}
        for (k, mergeCnt), patternLis in self.pattern.iteritems():
            curMergedPatternDict = {}
            for (newPattern, pattern1, pattern2) in patternLis:
                pattern2coverSeqLis[newPattern.patternStr] = set([seqId for seqId, _ in newPattern.visitLis])
                pattern1Info = (pattern1.patternStr, pattern1.coverage)
                pattern2Info = (pattern2.patternStr, pattern2.coverage)
                curMergedPatternDict.setdefault((newPattern.patternStr, newPattern.coverage), [])
                curMergedPatternDict[(newPattern.patternStr, newPattern.coverage)].append((pattern1Info, pattern2Info))
            topCoveragePatternLis = sorted(curMergedPatternDict.items(), key=lambda x:x[0][1], reverse=True)[:maxN]
            rltMatrix.append((k, mergeCnt, topCoveragePatternLis))
        return rltMatrix, pattern2coverSeqLis
    #end_func
# end_class

def test():
    # g = KmerGrpah(True)
    # n1 = CommStructure.Node('n1', [])
    # n2 = CommStructure.Node('n2', [])
    # e1 = CommStructure.Edge(n1, n2, None, 1, (1, 2))
    #
    # g.addNode(n1)
    # g.addNode(n2)
    # g.addEdge(e1)
    # g.updateLastNode(e1)
    # print g


    pattern1 = Pattern('[ACG]ACC', [(1,2)])
    pattern2 = Pattern('[AG][AT]CC', [(1,2)])
    newPattern = Pattern.mergePattern(pattern1, pattern2)
    print 1
# end_test

def main():
    test()


# end_main

if __name__ == "__main__":
    main()
# end_if
