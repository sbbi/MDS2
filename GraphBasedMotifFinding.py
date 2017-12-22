#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/5/14
Description     :
    cmd: python GraphBasedMotifFinding.py XX 6 ACGT 100000 False
    cmd: python GraphBasedMotifFinding.py XX 6 ACGT 100000 True

"""
import os
import pickle
import sys
import logging.config
logging.config.fileConfig('loggingConf.ini')
logger = logging.getLogger('root')
INFO = logger.info

if len(sys.argv) > 1:
    logFn = os.path.join(sys.argv[1], 'log.txt')
    logger.root.handlers[1].stream.close()
    logger.root.handlers[1].baseFilename = logFn
    logger.root.handlers[1].stream = open(logFn, 'w')

import BioinfoComm, Comm
from Comm import print0, print9, PrintWithTime
import DimerGraph
import SignificaceEvaluation
import HierarchicalGraph
import Backpack
import patternClustering
import outputRlt
import negativeTest

SINGLE_LENGTH_MODE = False
TARGET_LENGTH = 5

ENABLE_CPROFILE = False
ENABLE_PARALLEL = False
CORE_NUM = 16

# ENABLE_MIR_SAMPLING = False
ENABLE_MIR_SAMPLING = True

MIR_REF_FN = './data/allHumanMiRNA(HSA).fa'
# MIR_REF_FN = './data/miR_bta.fa'
RNA_REF_FN = './data/human_CDS.fasta'


def showPatternSimi(pattern2pwmSimi, pattern2cov, RNAPattern2pvalue, miRPattern2pvalue, patternPairInfo):
    patternWithEdgeSet = set()
    allPattern = set()
    for (pattern1, pattern2), simi in pattern2pwmSimi.iteritems():
        allPattern.add(pattern1)
        allPattern.add(pattern2)
        if simi > 0:
            patternWithEdgeSet.add(pattern1)
            patternWithEdgeSet.add(pattern2)

    for pattern in allPattern:
        cov = pattern2cov[pattern]
        RNAPvalue = RNAPattern2pvalue[pattern]
        miRPvalue = miRPattern2pvalue.get(pattern, 'N/A')
        outputLine = '%s\t%s\t%s\t%s\t%s' % (pattern, pattern in patternWithEdgeSet, cov, RNAPvalue, miRPvalue)
        INFO(outputLine)
#end_func

def rawFormat(fn):
    fn = fn.replace('[', '\[').replace(']', '\]').replace(' ', '\ ')
    return fn
#end_func


class GraphBasedMotifFinding:
    def __init__(self, baseDir, maxMotifLength, alphabet, samplingFreq, enableMiRSampling=True, enableNegSeq=False, fetchGroupMode=0, covPrecThres=0.1):
        self.maxMotifLength = int(maxMotifLength)
        self.alphabet = alphabet
        self.samplingFreq = int(samplingFreq)
        self.maxPathLength = 20 # max length when search for paths in 2mer graph
        self.backgroundRatio = {'A': 0.221, 'C': 0.235, 'G': 0.288, 'T': 0.256, 'U': 0.256}
        self.ghostscriptProcessor = 'gs-921-linux-x86_64'
        self.uniprobe2memeProcessor = './uniprobe2meme'
        self.tomtomProcessor = './tomtom'
        self.enableMiRSampling = enableMiRSampling
        self.enableNegSeq = enableNegSeq
        self.fetchGroupMode = fetchGroupMode
        self.covPrecThres = covPrecThres

        # dir info
        self.baseDir = baseDir.replace(' ', '\\ ')
        self.PPMDir = os.path.join(baseDir, 'PPM')
        self.uniprobeDir = os.path.join(baseDir, 'PPM')
        self.epsDir = os.path.join(baseDir, 'eps')
        self.pngDir = os.path.join(baseDir, 'png')
        folderLis = [self.PPMDir, self.epsDir, self.pngDir, self.uniprobeDir]

        # input file
        self.userDataFn = os.path.join(baseDir, 'userInput.fa')
        self.negSeqFn = os.path.join(baseDir, 'negSeq.fa')

        # output file
        self.inputInfoFn = os.path.join(baseDir, 'inputInfo.txt')
        self.pattern2pwmSimiFn = os.path.join(baseDir, 'pattern2pwmSimi')
        self.pattern2covSimiFn = os.path.join(baseDir, 'pattern2covSimi')
        self.pattern2covFn = os.path.join(baseDir, 'pattern2cov')
        self.pattern2kmerSetFn = os.path.join(baseDir, 'pattern2kmerSet.txt')
        self.patternInfoFn = os.path.join(baseDir, 'patternInfo.txt')
        self.clusterInfoFn = os.path.join(baseDir, 'clusterInfo.txt')
        self.finalRltFn = os.path.join(baseDir, 'finalRlt.txt')
        self.jobDoneFn = os.path.join(baseDir, 'job Done!!')

        # global variables
        self.overallSigKmer = [{}, {}] if self.enableMiRSampling else [{}]
        self.wholeTree = [[], []] if self.enableMiRSampling else [[]]
        self.allInTreeKmerSet = [set(), set()] if self.enableMiRSampling else [set()]
        self.overallKmer2Cov = {}
        self.overallRefKmer2CovThreshold = [{}, {}] if self.enableMiRSampling else [{}]
        # self.overallUserKmer2VisDetail = {}
        self.overallRefKmer2covIdxMatrix = [{}, {}] if self.enableMiRSampling else [{}]
        self.layerCovLis = []
        self.K2kmers = {}
        self.pattern2kmerSet = {}

        # open file
        self.patternInfoFile = open(self.patternInfoFn, 'w')
        self.clusterInfoFile = open(self.clusterInfoFn, 'w')
        self.finalRltFile = open(self.finalRltFn, 'w')

        # initialization functions
        if os.path.exists(self.jobDoneFn): os.remove(self.jobDoneFn)
        Comm.CreateFolder(folderLis)
        self.ShowBasicInfo()
    #end_init

    def __del__(self):
        self.patternInfoFile.close()
        self.clusterInfoFile.close()
        self.finalRltFile.close()
    # end_del

    def ShowBasicInfo(self):
        PrintWithTime("start...", isLogging=True)
        INFO("current data file: %s" % self.userDataFn)
        INFO("miRNA ref data file: %s" % MIR_REF_FN)
        INFO("RNA ref data file: %s" % RNA_REF_FN)
        INFO("sampling count: %s" % self.samplingFreq)
        INFO("enable MiR Sampling: %s" % self.enableMiRSampling)
        INFO("enable negative test: %s" % self.enableNegSeq)
        INFO("enable single length mode: %s" % SINGLE_LENGTH_MODE)
        INFO("mode ID for fetching core in groups: %s" % self.fetchGroupMode)
        INFO("coverage percentage threshold: %s" % self.covPrecThres)
        if SINGLE_LENGTH_MODE: INFO("single length: %s" % TARGET_LENGTH)
        print9(isLogging=True)
    #end_func

    # @profile
    def pipeline(self):
        # TODO change
        if SINGLE_LENGTH_MODE:
            motifLengthLis = [TARGET_LENGTH]
        else:
            motifLengthLis = range(self.maxMotifLength + 1)[2:]  # 2 <= K < maxMotifLength

        # load user data, build 2-mer graph and fetch path
        userSeqLis, seqCnt, minSeqLen, maxSeqLen = BioinfoComm.loadSinglelineSeq(self.userDataFn)
        kmerMatrix = DimerGraph.DivideSeq(userSeqLis)
        dimerGraph = DimerGraph.BuildGraph(kmerMatrix, alphabet=self.alphabet)
        K2Paths = DimerGraph.SearchPaths(dimerGraph, self.maxPathLength, enableFiltering=False, covPercThres=self.covPrecThres)
        self.overallKmer2Cov[1] = {'root': seqCnt}
        self.layerCovLis.append(('root', seqCnt))  # root layer covers all sequences

        # load ref data
        sampledSeqMatrix = []
        RNASampledSeqMatrix = SignificaceEvaluation.SampleRNARef(RNA_REF_FN, seqCnt, self.samplingFreq, minSeqLen, maxSeqLen, self.alphabet)
        sampledSeqMatrix.append(RNASampledSeqMatrix)
        if self.enableMiRSampling:
            miRSampledSeqMatrix = SignificaceEvaluation.SampleMiRRef(MIR_REF_FN, seqCnt, self.samplingFreq, self.alphabet)
            sampledSeqMatrix.append(miRSampledSeqMatrix)

        PrintWithTime("finish loading data", isLogging=True)
        INFO("number of sequences: %s" % seqCnt)
        INFO("min length: %s" % minSeqLen)
        INFO("max length: %s" % maxSeqLen)

        print9(isLogging=True)

        kmer2visDetail = {}
        for motifLength in motifLengthLis:
            PrintWithTime("begin to calculate on layer motifLength = %s" % motifLength, isLogging=True)

            for typeId, allSigKmer in enumerate(self.overallSigKmer):
                title = 'RNA' if typeId == 0 else 'miR'
                lastLayerKmers = allSigKmer if motifLength > 2 else BioinfoComm.GenerateKmerLis(self.alphabet, 2)
                inTreeKmerSet = self.allInTreeKmerSet[typeId]

                # not all the k-mers are sampled in order to increase the speed
                sigKmerSet, userKmer2cov, kmer2visDetail, kmer2refCovThreshold, kmer2Pvalue, _, refKmer2covIdxMatrix, inTreeKmerSet = SignificaceEvaluation.EvalSignificance(motifLength, sampledSeqMatrix[typeId], lastLayerKmers, seqCnt, self.samplingFreq, K2Paths, inTreeKmerSet, self.alphabet, enableParallel=ENABLE_PARALLEL, coreNm=CORE_NUM, covPrecThres=self.covPrecThres)
                filtSigKmerSet, self.wholeTree[typeId] = HierarchicalGraph.BuildTree(motifLength, self.wholeTree[typeId], sigKmerSet, lastLayerKmers, userKmer2cov, userSeqLis, covPrecThres=self.covPrecThres)
                HierarchicalGraph.DisplayLayerInfo(title, sigKmerSet, filtSigKmerSet, self.wholeTree[typeId], kmer2Pvalue)

                # update overall dict
                self.overallKmer2Cov.setdefault(motifLength, userKmer2cov)
                self.overallRefKmer2covIdxMatrix[typeId] = refKmer2covIdxMatrix  # [RNA, miR]
                self.overallRefKmer2CovThreshold[typeId] = dict(kmer2refCovThreshold, **self.overallRefKmer2CovThreshold[typeId]) # merge dict

                # TODO: change
                if SINGLE_LENGTH_MODE:
                    self.overallSigKmer[typeId] = filtSigKmerSet
                else:
                    self.overallSigKmer[typeId] = sigKmerSet

                self.allInTreeKmerSet[typeId] = inTreeKmerSet

            if not self.overallSigKmer[0] and (self.enableMiRSampling and not self.overallSigKmer[1]): continue

            treeStructureFn = os.path.join(self.baseDir, 'currentTree_K=%s' % motifLength)

            # TODO: change
            if SINGLE_LENGTH_MODE:
                segmentLis = sigKmerSet
            else:
                self.K2kmers = HierarchicalGraph.DrawTree(treeStructureFn, self.wholeTree, self.overallKmer2Cov, self.overallRefKmer2CovThreshold, enableMiRSampleing=self.enableMiRSampling)
                print0(isLogging=True)

                if motifLength < 2: continue

                # display all the segments
                segmentLis = self.K2kmers[motifLength][0] | self.K2kmers[motifLength][1]

            PrintWithTime("all the segments", isLogging=True)
            INFO(segmentLis)
            print9(isLogging=True)

            # build segment similarity graph
            INFO('fetching segment similarity graph')
            graphFn = os.path.join(self.baseDir, 'segmentSimiGraph-K=%s' % motifLength)
            simiGraph = Backpack.BuildSimiGraph(graphFn, segmentLis)

            # format covered sequences into binary number
            allNode = set(simiGraph.node_neighbors.keys())
            kmer2seqIdSet = {kmer: set(map(lambda x: x[0], visLis)) for kmer, visLis in kmer2visDetail.iteritems()}
            userKmer2seqIdInt = BioinfoComm.formatCovId(kmer2seqIdSet, seqCnt)

            # build backpack model to find the path based on segments
            INFO('building backpack problem')
            initMinIC = (0.61 * 0.25 + 0.3 * 0.5 + 0.125 * 0.25) * motifLength
            initMaxIC = 0.61 * motifLength # 0 variables in each position
            allPatterns, userKmer2seqIdInt, pattern2IC, pattern2kmerSet = Backpack.FetchPatternWithICIter(allNode, simiGraph.node_neighbors, initMinIC, initMaxIC, userSeqLis, set(), userKmer2seqIdInt, self.overallKmer2Cov, pattern2IC={}, enableParallel=ENABLE_PARALLEL, coreNm=CORE_NUM, pattern2kmerSet={})

            self.pattern2kmerSet = dict(self.pattern2kmerSet, **pattern2kmerSet)

            PrintWithTime("all the patterns", isLogging=True)
            INFO(allPatterns)
            print0(isLogging=True)

            # fetch pattern's IC, coverage and p-value to do filtering and sorting
            pattern2cov, _ = outputRlt.FetchPatternSetCov(allPatterns, userSeqLis, userKmer2seqIdInt, seqCnt)
            PrintWithTime('calculating Pvalue of RNA', isLogging=True)
            RNAPattern2pvalue, disabledRNAPattern = outputRlt.FetchPatternPvalue(initMinIC, seqCnt, pattern2IC, pattern2cov, sampledSeqMatrix[0], self.overallRefKmer2covIdxMatrix[0], self.samplingFreq, enableParallel=ENABLE_PARALLEL, coreNm=CORE_NUM)
            self.overallRefKmer2covIdxMatrix[0] = {}
            print9(isLogging=True)

            miRPattern2pvalue = {}
            if self.enableMiRSampling:
                PrintWithTime('calculating Pvalue of miR', isLogging=True)
                miRPattern2pvalue, _ = outputRlt.FetchPatternPvalue(initMinIC, seqCnt, pattern2IC, pattern2cov, sampledSeqMatrix[1], self.overallRefKmer2covIdxMatrix[1], self.samplingFreq, disabledInputPattern=disabledRNAPattern, enableParallel=ENABLE_PARALLEL, coreNm=CORE_NUM)
                print9(isLogging=True)
            self.overallRefKmer2covIdxMatrix = [{}, {}] if self.enableMiRSampling else [{}]

            # format the pattern result and merge all the dict
            patternInfoLis = []
            for sourcePattern, userCov in pattern2cov.iteritems():
                if sourcePattern not in RNAPattern2pvalue or sourcePattern not in RNAPattern2pvalue: continue
                IC = pattern2IC[sourcePattern]
                RNAPvalue = 1.0 if RNAPattern2pvalue[sourcePattern] > 1 else RNAPattern2pvalue[sourcePattern]
                if RNAPvalue > 0.05: continue
                miRPvalue = 'N/A' if not self.enableMiRSampling else 1.0 if miRPattern2pvalue[sourcePattern] > 1 else miRPattern2pvalue[sourcePattern]
                patternInfoLis.append((sourcePattern, IC, userCov, RNAPvalue, miRPvalue))

            if not patternInfoLis: continue
            patternInfoLis.sort(key=lambda x: (x[2], x[1]), reverse=True)

            # filter negative pattern
            if self.enableNegSeq:
                patternCovLis = map(lambda x:(x[0], x[2]), patternInfoLis)
                allKmerSet = BioinfoComm.GenerateKmerLis(self.alphabet, motifLength)
                signKmerSet = negativeTest.TestPatternInNegSeq(patternCovLis, seqCnt, self.negSeqFn, allKmerSet)
                patternInfoLis = filter(lambda x:x[0] in signKmerSet, patternInfoLis)

            # display motif pattern result
            PrintWithTime("IC, cov and pvalue of all the patterns", isLogging=True)
            Comm.showList(patternInfoLis, isLogging=True)
            self.patternInfoFile.write('=====motif length: %s=====\n' % motifLength)
            titleLine = '\t'.join(['pattern', 'IC', 'coverage(in %s)' % len(userSeqLis), 'RNA_Pvalue', 'miR_Pvalue'])
            self.patternInfoFile.write('%s\n' % titleLine)
            for item in patternInfoLis:
                curLine = '\t'.join(map(lambda x: str(x), item))
                self.patternInfoFile.write('%s\n' % curLine)
            print0(isLogging=True)

            # fetch PPM
            PrintWithTime("begin to fetch PPM", isLogging=True)
            curPPMDir = os.path.join(self.PPMDir, str(motifLength))
            curUniprobeDir = os.path.join(self.uniprobeDir, str(motifLength))
            curEpsDir = os.path.join(self.epsDir, str(motifLength))
            curPngDir = os.path.join(self.pngDir, str(motifLength))
            Comm.CreateFolder([curPPMDir, curEpsDir, curPngDir])

            # output PPM and uniprobe for tomtom
            PPMFn = os.path.join(curPPMDir, 'PPM.txt')
            uniprobeFn = os.path.join(curUniprobeDir, 'uniprobe.txt')
            pattern2PPM, uniprobeFnLis = outputRlt.OutputPPM(motifLength, patternInfoLis, PPMFn, uniprobeFn, userSeqLis, self.alphabet, self.backgroundRatio)

            if len(patternInfoLis) <= 1:continue

            # output pattern similarity to do clustering
            # run uniprobe2meme
            memeFnLis = []
            for uniprobeFn in uniprobeFnLis:
                basename, ext = os.path.splitext(uniprobeFn)
                memeFn = '%s.meme' % basename
                memeFnLis.append(memeFn)
                bgFn = 'data/background_rna' if 'U' in self.alphabet else 'data/background_dna'
                rnaPara = '-rna' if 'U' in self.alphabet else ''
                cmd = '%s %s -bg %s %s > %s' % (self.uniprobe2memeProcessor, rnaPara, bgFn, rawFormat(uniprobeFn), rawFormat(memeFn))
                INFO(cmd)
                os.system(cmd)

            # merge all the motif info(*.meme) into one file
            allMotifFn = '%s/allMotif.meme' % os.path.dirname(uniprobeFn)
            hasHeader = False
            patternMap = {}
            mergedPatternLis = []
            with open(allMotifFn, 'w') as allMotifFileobj:
                for memeFn in memeFnLis:
                    sourcePattern = os.path.basename(memeFn).strip('uniprobe-').partition('.')[0]
                    with open(memeFn) as memeFileobj:
                        isMotifSection = False
                        for line in memeFileobj:
                            if line[:5] == 'MOTIF':
                                isMotifSection = True
                                encodePattern = line.strip().partition(' ')[2]

                                mergedPatternLis.append(encodePattern)

                                patternMap[encodePattern] = sourcePattern
                            if not hasHeader or isMotifSection:
                                allMotifFileobj.write(line)
                        hasHeader = True

            """
            # run tomtom 1 vs. 1, this section is replaced by running tomtom 1 vs. all, in which pvalue is more accurate 
            patternPairInfo = {}
            for memeFn1 in memeFnLis:
                for memeFn2 in memeFnLis:
                    if memeFn1 <= memeFn2: continue
                    pattern1 = os.path.basename(memeFn1).strip('uniprobe-').partition('.')[0]
                    pattern2 = os.path.basename(memeFn2).strip('uniprobe-').partition('.')[0]
                    cmd = '%s -text -no-ssc -oc . -verbosity 1 -min-overlap 2 -mi 1 -dist pearson -evalue -thresh 10.0 %s %s' % (self.tomtomProcessor, rawFormat(memeFn1), rawFormat(memeFn2))
                    INFO(cmd)
                    rltLine = os.popen(cmd).readlines()[-1]
                    items = rltLine.strip().split('\t')
                    offset = int(items[2])
                    pvalue = float(items[3])
                    evalue = float(items[4])
                    qvalue = float(items[5])
                    patternPairInfo[(pattern1, pattern2)] = (offset, pvalue, evalue, qvalue)
                    patternPairInfo[(pattern2, pattern1)] = (-1 * offset, pvalue, evalue, qvalue)
            """

            # run tomtom 1 vs. all
            allMotifTomtomRltFn = '%s/allMotifTomtomRlt.meme' % os.path.dirname(uniprobeFn)
            with open(allMotifTomtomRltFn, 'w') as allMotifTomtomRltFileobj:
                patternPairInfo = {}
                for memeFn in memeFnLis:
                    sourcePattern = os.path.basename(memeFn).strip('uniprobe-').partition('.')[0]
                    cmd = '%s -text -no-ssc -oc . -verbosity 1 -min-overlap %s -mi 1 -dist pearson -evalue -thresh 10.0 %s %s' % (self.tomtomProcessor, motifLength / 2, rawFormat(memeFn), allMotifFn)
                    INFO(cmd)
                    rltLines = os.popen(cmd).readlines()
                    isMotifSection = False
                    for rltLine in rltLines:
                        INFO(rltLine)
                        allMotifTomtomRltFileobj.write('%s' % rltLine)
                        if rltLine[0] == '#': # only read motif section which start with '#'
                            isMotifSection = True
                            continue
                        elif not isMotifSection:
                            continue
                        items = rltLine.strip().split('\t')
                        targetPattern = items[1]
                        offset = int(items[2])
                        pvalue = float(items[3])
                        evalue = float(items[4])
                        qvalue = float(items[5])
                        patternPairInfo[(sourcePattern, targetPattern)] = (offset, pvalue, evalue, qvalue)
                        patternPairInfo[(targetPattern, sourcePattern)] = (-1 * offset, pvalue, evalue, qvalue)

            # node's weight: pattern and their coverage
            PrintWithTime("begin to fetch node weight(pattern coverage)", isLogging=True)
            pattern2covFn = '%s-K=%s.txt' % (self.pattern2covFn, motifLength)
            pattern2cov = patternClustering.outputPatternCov(patternInfoLis, pattern2covFn)

            # edge's weight: pattern and their similarity based on tomtom
            PrintWithTime("begin to fetch edge weight(pattern similarity)", isLogging=True)
            pattern2simiFn = '%s-K=%s.txt' % (self.pattern2pwmSimiFn, motifLength)
            pattern2pwmSimi = patternClustering.CalcPatternPwmSimi(patternInfoLis, pattern2PPM, patternPairInfo, pattern2simiFn, self.alphabet, doNormalization=False)

            # print data to show pattern simi edge, cov and pvalue
            INFO('** pattern info: simi, coverage and pvalue **')
            showPatternSimi(pattern2pwmSimi, pattern2cov, RNAPattern2pvalue, miRPattern2pvalue, patternPairInfo)
            INFO('***' * 20)

            if not pattern2pwmSimi: continue

            """
            # build pattern similarity graph and do clustering, this section is replace by R script below
            PrintWithTime("begin to build pattern similarity graph", isLogging=True)
            clusterLis = patternClustering.patternClustering(pattern2cov, pattern2pwmSimi)
            """
            # run R script to build similarity graph and do clustering
            PrintWithTime("begin to build pattern similarity graph", isLogging=True)
            clusterRltFn = '%s/clusterRlt-K=%s.txt' % (os.path.dirname(pattern2simiFn), motifLength)
            cmd = "Rscript clustering.R %s %s" % (pattern2simiFn, clusterRltFn)
            INFO(cmd)
            os.system(cmd)
            clusterLis = patternClustering.loadClusterRlt(clusterRltFn)
            if not clusterLis: continue

            # show cluster info
            outputLine = '==== clusters for motif with length: %s ====' % motifLength
            self.clusterInfoFile.write('%s\n' % outputLine)
            INFO(outputLine)
            for cluster in clusterLis:
                outputLine = '\t'.join(cluster)
                self.clusterInfoFile.write('%s\n' % outputLine)
                INFO(outputLine)

            # for each cluster, find some cores
            PrintWithTime("begin to fetch cluster core", isLogging=True)
            coreMotifSet, coreOutputLines = patternClustering.searchForCoreInCluster(clusterLis, userKmer2seqIdInt, seqCnt, pattern2cov, pattern2IC, RNAPattern2pvalue, miRPattern2pvalue, self.fetchGroupMode, covPrecThres=self.covPrecThres)

            # show cores in each cluster, which are the final motif to report
            outputLine = '== core pattern in each cluster =='
            self.clusterInfoFile.write('%s\n' % outputLine)
            INFO(outputLine)
            for coreOutputLine in coreOutputLines:
                self.clusterInfoFile.write('%s\n' % coreOutputLine)
                INFO(coreOutputLine)

            # output final result
            INFO('== final result ==')
            outputLine = 'pattern length: %s' % motifLength
            self.finalRltFile.write('%s\n' % outputLine)
            INFO(outputLine)
            coreMotifLis = sorted(coreMotifSet, key=lambda x:(pattern2cov[x], pattern2IC[x]), reverse=True)
            for motif in coreMotifLis:
                outputLis = [motif, str(pattern2cov[motif]), str(pattern2IC[motif]), str(RNAPattern2pvalue[motif]), str(miRPattern2pvalue.get(motif, 'N/A'))]
                outputLine = '\t'.join(outputLis)
                self.finalRltFile.write('%s\n' % outputLine)
                INFO(outputLine)

            # output logo
            PrintWithTime("begin to output logo", isLogging=True)
            epsFn = os.path.join(curEpsDir, 'logo.eps')
            pngFn = os.path.join(curPngDir, 'logo.png')
            outputRlt.OutputLogo(coreMotifLis, motifLength, self.ghostscriptProcessor, PPMFn, epsFn, pngFn)

            # showing finishing a loop
            PrintWithTime("finish the loop for K=%s" % motifLength, isLogging=True)
        # end_for for current K

        # write pattern2kmerSet dict
        with open(self.pattern2kmerSetFn, 'w') as pattern2kmerSetFileobj:
            pickle.dump(self.pattern2kmerSet, pattern2kmerSetFileobj)

        PrintWithTime("finish the pipeline", isLogging=True)
        with open(self.jobDoneFn, 'w') as f: pass
    # end_func
#end_class

def main():
    if len(sys.argv) == 1: # test in Windows
        # baseDir = './output/testT'
        baseDir = './output/testU'
        # baseDir = r'D:\project\MotifFinding\kmerGraph\GraphBasedMotifFinding\output\par_clip'
        # baseDir = r'D:\project\MotifFinding\kmerGraph\GraphBasedMotifFinding\output\sw620_39'
        # baseDir = r'D:\project\MotifFinding\kmerGraph\GraphBasedMotifFinding\output\cell2016_seq103'
        # baseDir = r'D:\project\MotifFinding\kmerGraph\GraphBasedMotifFinding\output\PAR-CLIP\head_50'
        # baseDir = r'D:\project\MotifFinding\kmerGraph\GraphBasedMotifFinding\output\mRNA_negTest'
        maxMotifLength = 6
        alphabet = 'ACGU'
        samplingFreq = 100
        fetchGroupMode = 2
        enableNegSeq = False
        enableMiRSampling = ENABLE_MIR_SAMPLING
        covPrecThres = 0.1
        g = GraphBasedMotifFinding(baseDir, maxMotifLength, alphabet, samplingFreq, enableMiRSampling, enableNegSeq, fetchGroupMode, covPrecThres)
        g.pipeline()
    else:
        _, baseDir, maxMotifLength, alphabet, samplingFreq, enableMiRSampling, enableNegSeq, fetchGroupMode, covPrecThres = sys.argv
        g = GraphBasedMotifFinding(baseDir, maxMotifLength, alphabet, samplingFreq, enableMiRSampling=eval(enableMiRSampling), enableNegSeq=eval(enableNegSeq), fetchGroupMode=int(fetchGroupMode), covPrecThres=float(covPrecThres))
        g.pipeline()
#end_main

if __name__ == "__main__":
    if ENABLE_CPROFILE:
        import cProfile
        import pstats

        cp = cProfile.Profile()
        cp.runcall(main)
        ps = pstats.Stats(cp)
        sorted_stats = ps.sort_stats('time')
        sorted_stats.print_stats(10)
    else:
        main()
#end_if
