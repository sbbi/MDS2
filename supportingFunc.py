#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     :
Author          : Tian Gao
CreationDate    : 2017/6/12
Description     :

"""
import sys
import BioinfoComm, Comm
import os
import re
import patternClustering
import DimerGraph
import SignificaceEvaluation
import outputRlt
import patternSimi

ENABLE_CPROFILE = False

def FetchMilkRNAData():
    idxFn = r"D:\Dropbox\Motif analysis\2017-6-16 cow's milk\segment idx.txt"
    oriSeqFn = r"C:\Users\Administrator\Desktop\bta_rna.fa"
    outputFn = r"C:\Users\Administrator\Desktop\milkSeq.fa"

    # load idx file
    allIdxLis = []
    with open(idxFn) as idxFileobj:
        idxFileobj.readline()
        for line in idxFileobj:
            curIdxLis = []
            items = line.strip().split('\t')[1:]
            for indexLisStr in items:
                indexLis = indexLisStr.split('|')
                for indexStr in indexLis:
                    nm, _, posStr = indexStr.partition(':')
                    startPos, _, endPos = posStr.partition('-')
                    curIdxLis.append((nm, int(startPos), int(endPos)))
            minLenIdx = min(curIdxLis, key=lambda x:x[2]-x[1])
            allIdxLis.append(minLenIdx)
    print max(allIdxLis, key=lambda x:x[2]-x[1])
    print min(allIdxLis, key=lambda x:x[2]-x[1])

    # load whole sequences file and format it
    seqLis, seqCnt, titleLis = BioinfoComm.loadMultilineSeq(oriSeqFn)
    # Comm.showList(titleLis)
    seqNm2seq = {}
    for id, title in enumerate(titleLis):
        seq = seqLis[id]
        seqNm = title.split('|')[3].partition('.')[0]
        seqNm2seq[seqNm] = seq

    # get segment
    with open(outputFn, 'w') as outputFileobj:
        for nm, startPos, endPos in allIdxLis:
            if nm not in seqNm2seq:
                print nm
                continue
            curSeq = seqNm2seq[nm][startPos: endPos + 1]
            outputLine = '>%s-%s\n%s\n' % (nm, len(curSeq), curSeq)
            outputFileobj.write(outputLine)
#end_func

def GenerateNegTestData():
    inputBaseDir = r'D:\Dropbox\Motif analysis\2017-6-7 RNA seq\miR_mRNA_binding site'
    outputBaseDir = r'C:\Users\Administrator\Desktop\1212312\1'
    clusterNm = ['seqInfo_I', 'seqInfo_II', 'seqInfo_III', 'seqInfo_IV', 'seqInfo_V']
    userDataNmLis = map(lambda x:os.path.join(inputBaseDir, x, 'userInput.fa'), clusterNm)
    userDataLis = map(lambda x:open(x).readlines(), userDataNmLis)

    for posId, posUserData in enumerate(userDataLis):
        negIdLis = [str(x + 1) for x in xrange(len(clusterNm)) if x != posId]
        tmpDir = 'pos_%s-neg_%s' % (posId + 1, ','.join(negIdLis))
        outputDir = os.path.join(outputBaseDir, tmpDir)
        if not os.path.exists(outputDir): os.makedirs(outputDir)
        posDataFn = os.path.join(outputDir, 'userInput.fa')
        negDataFn = os.path.join(outputDir, 'negSeq.fa')
        with open(posDataFn, 'w') as posDataFileobj, open(negDataFn, 'w') as negDataFileobj:
            for line in posUserData:
                posDataFileobj.write(line)
            for negId, negUserData in enumerate(userDataLis):
                if negId == posId: continue
                for line in negUserData:
                    negDataFileobj.write(line)
#end_func

def SeqDiffExo():
    allExoFn = r'D:\Dropbox\project\MotifFinding\sources\old\allmiRNA.fa'
    btaOutputFn = r'C:\Users\Administrator\Desktop\1212312\1\miR_bta.fa'
    hsaOutputFn = r'C:\Users\Administrator\Desktop\1212312\1\miR_hsa.fa'
    mmuOutputFn = r'C:\Users\Administrator\Desktop\1212312\1\miR_mmu.fa'

    with open(btaOutputFn, 'w') as btaOutputFileobj, open(hsaOutputFn, 'w') as hsaOutputFileobj, open(mmuOutputFn, 'w') as mmuOutputFileobj:
        seqLis, _, titleLis = BioinfoComm.loadMultilineSeq(allExoFn)
        for seqId, title in enumerate(titleLis):
            if title[:3] == 'bta':
                curFileobj = btaOutputFileobj
            elif title[:3] == 'hsa':
                curFileobj = hsaOutputFileobj
            elif title[:3] == 'mmu':
                curFileobj = mmuOutputFileobj
            else:
                continue
            seq = seqLis[seqId]
            outputLine = '>%s\n%s\n' % (title, seq)
            curFileobj.write(outputLine)
#end_func




def formatJiangData():
    inputFn = r'C:\Users\Administrator\Desktop\jiangData\all.txt'
    posSeqFn = r'C:\Users\Administrator\Desktop\jiangData\userInput.fa'
    negSeqFn = r'C:\Users\Administrator\Desktop\jiangData\negSeq.fa'

    with open(inputFn) as inputFileobj, open(posSeqFn, 'w') as posSeqFileobj, open(negSeqFn, 'w') as negSeqFileobj:
        inputFileobj.readline()
        for line in inputFileobj:
            items = line.strip().split('\t')
            title = items[0]
            seq = items[1]
            weight = int(items[2])
            if weight >= 30:
                outputLine = '>%s\n%s\n' % (title, seq)
                posSeqFileobj.write(outputLine)
            elif weight < 30:
                outputLine = '>%s\n%s\n' % (title, seq)
                negSeqFileobj.write(outputLine)
#end_func

def FetchPPMDic(seqLis, patternLis, motifLength, alphabet='ACGU'):
    pattern2PPM = {}
    for pattern in patternLis:
        pattern = pattern.replace('T', 'U')
        matchMatrix = [re.findall(pattern, userSeq.replace('T', 'U')) for userSeq in seqLis]
        seqWeightLis = [1.0 / len(matchLis) if len(matchLis) > 0 else 0 for matchLis in matchMatrix]
        PPM = [{curChar: 0 for curChar in alphabet} for _ in xrange(motifLength)]

        for seqIdx, matchLis in enumerate(matchMatrix):
            if len(matchLis) == 0: continue
            weight = seqWeightLis[seqIdx]
            for matchedSegment in matchLis:
                for charIdx, curChar in enumerate(matchedSegment):
                    PPM[charIdx][curChar] += weight

        pattern2PPM[pattern] = PPM
    return pattern2PPM
#end_func

def FetchDiffClusterSimi():
    """
    fetch similarity of same RNA with different clusters
    """
    titleLis = ['cluster1', 'cluster2', 'cluster3', 'cluster4', 'cluster5']
    inputDir = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result(unfinished)\result\normal'
    outputDir = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result(unfinished)\pattern similarity'

    for title1 in titleLis:
        for title2 in titleLis:
            if title1 <= title2: continue
            dataset1Dir = os.path.join(inputDir, title1)
            dataset2Dir = os.path.join(inputDir, title2)
            curOutputDir = os.path.join(outputDir, '%s_%s' % (title1, title2))
            if not os.path.exists(curOutputDir): os.makedirs(curOutputDir)
            FetchTwoDatasetSimi(4, dataset1Dir, title1, dataset2Dir, title2, curOutputDir)
#end_func

def FetchDiffClusterSimi_WithTomtom():
    """
    fetch similarity of RNA with same length in different dataset
    """
    inputInfoLis = [ # (title, path)
        # ('sw620_39', r'motifSimi/sw620_seq39'),
        # ('sw620_112', r'motifSimi/sw620_seq112'),
        ('nature', r'motifSimi/nature'),
        # ('cell', r'motifSimi/cell'),
        ('cowMilk_RNA', r'motifSimi/cowMilk_RNA'),
        ('cowMilk_miR', r'motifSimi/cowMilk_miR'),
    ]

    lenLis = [3, 4, 5, 6]

    outputSimiRltFn = r'motifSimi/allSimiRlt.txt'

    for motifLen in lenLis:
        pattern2PPM = {}

        # load PPM and merge all the motif in the same length
        allMotifFn = r'motifSimi/allMotif_%s.txt' % motifLen
        allMotifTomtomRltFn = r'motifSimi/allMotifTomtomRlt_%s.txt' % motifLen
        motifSimiFn = r'motifSimi/motifSimi_%s.txt' % motifLen
        hasHeader = False
        with open(allMotifFn, 'w') as allMotifFileobj:
            mergedPatternLis = []
            for title, basePath in inputInfoLis:
                rltFn = os.path.join(basePath, "finalRlt.txt")
                seqFn = os.path.join(basePath, "userInput.fa")
                seqLis, seqCnt, _ = BioinfoComm.loadMultilineSeq(seqFn)
                pattern2cov = patternSimi.loadRlt(rltFn, motifLen)
                for pattern, cov in pattern2cov.iteritems():
                    encodedPattern = '%s_%s' % (pattern, title)
                    memeFn = os.path.join(basePath, 'PPM', str(motifLen), 'uniprobe-%s.meme' % pattern)
                    PPMFn = os.path.join(basePath, 'PPM', str(motifLen), 'PPM-%s.txt' % pattern)

                    # load meme file
                    with open(memeFn) as memeFileobj:
                        isMotifSection = False
                        for line in memeFileobj:
                            if line[:5] == 'MOTIF':
                                isMotifSection = True
                                mergedPatternLis.append(encodedPattern)
                                allMotifFileobj.write('MOTIF %s' % encodedPattern)
                                continue
                            if not hasHeader or isMotifSection:
                                allMotifFileobj.write(line)
                        hasHeader = True

                    # load PPM file
                    with open(PPMFn) as PPMFileobj:
                        totalPPM = []
                        charLis = PPMFileobj.readline().strip().split('\t')[1:]
                        print charLis
                        for line in PPMFileobj:
                            linePPM = {}
                            items = line.strip().split('\t')[1:]
                            for idx, cnt in enumerate(items):
                                linePPM[charLis[idx]] = float(cnt)
                            totalPPM.append(linePPM)
                        pattern2PPM[encodedPattern] = totalPPM

        # run tomtom 1 vs. all
        with open(allMotifTomtomRltFn, 'w') as allMotifTomtomRltFileobj:
            patternPairInfo = {}
            for title, basePath in inputInfoLis:
                rltFn = os.path.join(basePath, "finalRlt.txt")
                seqFn = os.path.join(basePath, "userInput.fa")
                seqLis, seqCnt, _ = BioinfoComm.loadMultilineSeq(seqFn)
                pattern2cov = patternSimi.loadRlt(rltFn, motifLen)
                for pattern, cov in pattern2cov.iteritems():
                    singleMotifFn = os.path.join(basePath, 'PPM', str(motifLen), 'uniprobe-%s.meme' % pattern)
                    cmd = './tomtom -text -no-ssc -oc . -verbosity 1 -min-overlap %s -mi 1 -dist pearson -evalue -thresh 10.0 %s %s' % (motifLen / 2, Comm.rawFormat(singleMotifFn), allMotifFn)
                    print cmd
                    rltLines = os.popen(cmd).readlines()
                    isMotifSection = False
                    for rltLine in rltLines:
                        if rltLine[0] == '#' and not isMotifSection:  # only read motif section which start with '#'
                            isMotifSection = True
                            allMotifTomtomRltFileobj.write('%s' % rltLine)
                            continue
                        elif not isMotifSection:
                            continue
                        items = rltLine.strip().split('\t')
                        sourcePattern = items[0]
                        sourcePattern = '%s_%s' % (sourcePattern, title)
                        targetPattern = items[1]
                        offset = int(items[2])
                        pvalue = float(items[3])
                        evalue = float(items[4])
                        qvalue = float(items[5])
                        allMotifTomtomRltFileobj.write('%s\t%s\n' % (sourcePattern, '\t'.join(items[1:])))
                        patternPairInfo[(sourcePattern, targetPattern)] = (offset, pvalue, evalue, qvalue)
                        patternPairInfo[(targetPattern, sourcePattern)] = (-1 * offset, pvalue, evalue, qvalue)

            # edge's weight: pattern and their similarity based on tomtom
            pattern2pwmSimi = patternClustering.CalcPatternPwmSimi([], pattern2PPM, patternPairInfo, motifSimiFn, 'ACGU', doNormalization=False)
            print pattern2pwmSimi


#end_func

def FetchTwoDatasetSimi(motifLength, dataset1Dir, dataset1Title, dataset2Dir, dataset2Title, outputDir):
    # output file
    pattern2covFn = os.path.join(outputDir, "pattern2cov_%s-%s_length%s.txt" % (dataset1Title, dataset2Title, motifLength))
    pattern2simiFn = os.path.join(outputDir, "pattern2simi_%s-%s_length%s.txt" % (dataset1Title, dataset2Title, motifLength))

    # fetch rlt and ppm from dataset1
    seqFn1 = os.path.join(dataset1Dir, "userInput.fa")
    rltFn1 = os.path.join(dataset1Dir, "finalRlt.txt")
    seqLis1, seqCnt1, _ = BioinfoComm.loadMultilineSeq(seqFn1)
    pattern2cov1 = patternSimi.loadRlt(rltFn1, motifLength)
    rltLis1 = pattern2cov1.keys()
    pattern2PPM1 = FetchPPMDic(seqLis1, rltLis1, motifLength)
    pattern2PPM1 = {'%s_%s' % (pattern, dataset1Title): ppm for pattern, ppm in pattern2PPM1.iteritems()}
    pattern2cov1 = {'%s_%s' % (pattern, dataset1Title): cov * 1.0 / seqCnt1 for pattern, cov in pattern2cov1.iteritems()}

    # fetch rlt and ppm from dataset2
    seqFn2 = os.path.join(dataset2Dir, "userInput.fa")
    rltFn2 = os.path.join(dataset2Dir, "finalRlt.txt")
    seqLis2, seqCnt2, _ = BioinfoComm.loadMultilineSeq(seqFn2)
    pattern2cov2 = patternSimi.loadRlt(rltFn2, motifLength)
    rltLis2 = pattern2cov2.keys()
    pattern2PPM2 = FetchPPMDic(seqLis2, rltLis2, motifLength)
    pattern2PPM2 = {'%s_%s' % (pattern, dataset2Title): ppm for pattern, ppm in pattern2PPM2.iteritems()}
    pattern2cov2 = {'%s_%s' % (pattern, dataset2Title): cov * 1.0 / seqCnt2 for pattern, cov in pattern2cov2.iteritems()}

    # merge cov dic and ppm dic
    pattern2cov = dict(pattern2cov1, **pattern2cov2)
    pattern2ppm = dict(pattern2PPM1, **pattern2PPM2)

    # fetch pattern pair similarity
    pattern2simi = {}
    for pattern1, PPM1 in pattern2ppm.iteritems():
        for pattern2, PPM2 in pattern2ppm.iteritems():
            if pattern1 <= pattern2: continue
            simi = patternClustering.PCC(PPM1, PPM2)
            pattern2simi[(pattern1, pattern2)] = simi

    # normalization
    pattern2simi = patternClustering.NormalizeSimi(pattern2simi)

    # output
    inClusterInfo = []
    transClusterInfo = []
    with open(pattern2simiFn, 'w') as pattern2simiFileobj, open(pattern2covFn, 'w') as pattern2covFileobj:
        for pattern, cov in pattern2cov.iteritems():
            outputLine = '%s\t%s' % (pattern, cov)
            pattern2covFileobj.write('%s\n' % outputLine)

        for (pattern1, pattern2), simi in pattern2simi.iteritems():
            outputLine = '%s\t%s\t%s' % (pattern1, pattern2, simi)
            pattern2simiFileobj.write('%s\n' % outputLine)
            p1, _, rnaCluster1 = pattern1.partition('_')
            p2, _, rnaCluster2 = pattern2.partition('_')
            if rnaCluster1 == rnaCluster2:
                inClusterInfo.append((pattern1, pattern2, simi))
            else:
                transClusterInfo.append((pattern1, pattern2, simi))

    # analyze in-cluster and trans-cluster patterns to fetch most different pair
    transClusterInfo.sort(key=lambda x:x[2], reverse=True)
    transAverageSimi = reduce(lambda x,y:x+y, map(lambda x:x[2], transClusterInfo)) * 1.0 / len(transClusterInfo)
    inClusterInfo.sort(key=lambda x:x[2], reverse=True)
    inAverageSimi = reduce(lambda x,y:x+y, map(lambda x:x[2], inClusterInfo)) * 1.0 / len(transClusterInfo)
    # if transAverageSimi >= inAverageSimi: return
    print 'transClusterInfo:'
    print 'average simi: %.2f' % transAverageSimi
    Comm.showList(transClusterInfo)
    print 'inClusterInfo:'
    print 'average simi: %.2f' % inAverageSimi
    Comm.showList(inClusterInfo)
    Comm.print9()






    # # fetch top simi pattern
    # miR2RNASimiLis = {}
    # for (miRPattern, RNAPattern), simi in pattern2simi.iteritems():
    #     miR2RNASimiLis.setdefault(miRPattern, set())
    #     simi = pattern2simi[(miRPattern, RNAPattern)]
    #     miR2RNASimiLis[miRPattern].add((RNAPattern, simi))
    #
    # for miRPattern, simi in seq39Pattern2PPM.iteritems():
    #     miR2RNASimiLis[miRPattern] = sorted(miR2RNASimiLis[miRPattern], key=lambda x:x[1], reverse=True)
    #
    # # write into file
    # patternSimiLis = [(pattern[0], pattern[1], simi) for pattern, simi in pattern2simi.iteritems()]
    # patternSimiLis.sort(key=lambda x:x[2], reverse=True)
    # outputFn = r'C:\Users\Administrator\Desktop\allPairSimi.txt'
    # with open(outputFn, 'w') as outputFileobj:
    #     for pattern1, pattern2, simi in patternSimiLis:
    #         outputLine = '%s\t%s\t%s' % (pattern1, pattern2, simi)
    #         outputFileobj.write('%s\n' % outputLine)
    #
    # outputFn = r'C:\Users\Administrator\Desktop\miRPatternTopSimi.txt'
    # with open(outputFn, 'w') as outputFileobj:
    #     for miRPattern, patternLis in miR2RNASimiLis.iteritems():
    #         outputLine = '%s\t%s\t%s' % (miRPattern, patternLis[0][0], patternLis[0][1])
    #         outputFileobj.write('%s\n' % outputLine)
#end_func

def FetchGivenPatternPvalue():
    """
    given a pattern, fetch its pvalue
    """
    if len(sys.argv) > 1:
        _, pattern, userDataFn, samplingFreq, alphabet = sys.argv
        samplingFreq = int(samplingFreq)
    else:
        # pattern = '[AGU][AGU][AGU][ACG][AU][CG]'
        pattern = 'GG[ACGU][ACGU][ACGU][ACGU]G[UG][ACG][AC]'

        # userDataFn = r'D:\Dropbox\Motif analysis\final result\2. result part\2.1 validation result\nature\userInput.fa'
        # userDataFn = r'C:\Dropbox\Motif analysis\final result\2. result part_original\2.1 paper validation result\cell2016\userInput.fa'
        userDataFn = r'D:\Dropbox\Motif analysis\final result\2. result part_1st revision\2.3_miR-mRNA_binding_site_result\cluster2\userInput.fa'
        # userDataFn = r'D:\Dropbox\Motif analysis\final result\2. result part_1st revision\2.3_miR-mRNA_binding_site_result\cluster3\userInput.fa'
        # userDataFn = r'D:\Dropbox\Motif analysis\final result\2. result part_1st revision\2.3_miR-mRNA_binding_site_result\cluster4\userInput.fa'

        # samplingFreq = 100000
        samplingFreq = 100
        alphabet = r'ACGU'

    targetPatternLis = [pattern]
    RNA_REF_FN = r'./data/human_CDS.fasta'
    motifLength = BioinfoComm.getPatternLength(pattern)
    maxPathLength = motifLength

    # load user data, build 2-mer graph and fetch path
    userSeqLis, seqCnt, minSeqLen, maxSeqLen = BioinfoComm.loadSinglelineSeq(userDataFn)
    kmerMatrix = DimerGraph.DivideSeq(userSeqLis)
    dimerGraph = DimerGraph.BuildGraph(kmerMatrix, alphabet=alphabet)
    K2Paths = DimerGraph.SearchPaths(dimerGraph, maxPathLength)

    # load ref data
    RNASampledSeqMatrix = SignificaceEvaluation.SampleRNARef(RNA_REF_FN, seqCnt, samplingFreq, minSeqLen,
                                                             maxSeqLen, alphabet)

    print "finish loading data"
    print "number of sequences: %s" % seqCnt
    print "min length: %s" % minSeqLen
    print "max length: %s" % maxSeqLen
    print "begin to calculate on layer motifLength = %s" % motifLength

    userKmer2Cov, kmer2VisDetail = DimerGraph.PathDic2PathCov(K2Paths, motifLength)
    enabledKmerSet = set()
    for curPattern in targetPatternLis:
        kmerLis = BioinfoComm.FetchAllKmerFromPattern(curPattern)
        enabledKmerSet = set(kmerLis) | enabledKmerSet
    kmer2refDistri, kmer2refCovThreshold, refKmer2covIdxMatrix, signKmerSet = SignificaceEvaluation.fetchSignKmerSet(motifLength, RNASampledSeqMatrix, enabledKmerSet, seqCnt, userKmer2Cov)
    userKmer2Cov, kmer2VisDetail = DimerGraph.PathDic2PathCov(K2Paths, motifLength)

    kmer2seqIdSet = {kmer: set(map(lambda x: x[0], visLis)) for kmer, visLis in kmer2VisDetail.iteritems()}
    userKmer2seqIdInt = BioinfoComm.formatCovId(kmer2seqIdSet, seqCnt)
    pattern2cov, _ = outputRlt.FetchPatternSetCov(targetPatternLis, userSeqLis, userKmer2seqIdInt, seqCnt)

    for curPattern in targetPatternLis:
        distriLis = outputRlt.FetchPatternCovDistri(curPattern, RNASampledSeqMatrix, refKmer2covIdxMatrix, seqCnt, samplingFreq)
        # rawPvalue = BioinfoComm.FetchPvalueFromBG(distriLis, pattern2cov[curPattern])
        rawPvalue = BioinfoComm.FetchPvalueFromBG_FindNonZero(distriLis, pattern2cov[curPattern])  # new pvalue method

        adjustedPvalue = min(rawPvalue * len(pattern2cov), 1.0)  # adjust P-value by multiplying the length
        print curPattern, adjustedPvalue


#end_func



def sepBindingSiteCluster():
    """
    generate data of different cluster as input for pipeline from raw data
    """
    inputFn = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result\more_mirnas\let7a\let7a_binding_sites_input.txt'
    # outputDir = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result\more_mirnas\mir16'
    outputDir = r'C:\Users\Administrator\Desktop\33'

    clusterIdxMap = {'I':1, 'II':2, 'III':3, 'IV': 4, 'V':5}
    tmpDic = {}
    with open(inputFn) as inputFileobj:
        for line in inputFileobj:
            items = line.strip().split('\t')
            seqNm = items[0]
            seq = items[1]
            cluster = clusterIdxMap[items[2]]
            tmpDic.setdefault(cluster, [])
            tmpDic[cluster].append((seqNm, seq))

    for cluster, seqLis in tmpDic.iteritems():
        if len(seqLis) < 30: print cluster
        curOutputDir = os.path.join(outputDir, 'cluster_%s' % cluster)
        if not os.path.exists(curOutputDir): os.makedirs(curOutputDir)
        curOutputFn = os.path.join(curOutputDir, 'userInput.fa')
        with open(curOutputFn, 'w') as curOutputFileobj:
            for (seqNm, seq) in seqLis:
                outputLine = '>%s\n%s' % (seqNm, seq)
                print outputLine
                curOutputFileobj.write('%s\n' % outputLine)
#end_func

def FetchDIffRNASimi():
    """
    fetch similarity of different RNA with different clusters
    """
    rnaLis = ['mir16', 'let7a', 'mir10a']
    for rna1 in rnaLis:
        for rna2 in rnaLis:
            if rna1 >= rna2: continue
            baseDataset1Dir = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result\more_mirnas\%s\clusterRlt' % (rna1)
            baseDataset2Dir = r'D:\Dropbox\Motif analysis\final result\2. result part\2.3 miR-mRNA binding site result\more_mirnas\%s\clusterRlt' % (rna2)
            baseOutputDir = r'C:\Users\Administrator\Desktop\similarity'
            motifLength = 4
            clusterLis = [1, 3, 5]
            for cluster1 in clusterLis:
                for cluster2 in clusterLis:
                    title1 = '%s_cluster%s' % (rna1, cluster1)
                    title2 = '%s_cluster%s' % (rna2, cluster2)
                    dataset1Dir = os.path.join(baseDataset1Dir, 'cluster_%s' % cluster1)
                    dataset2Dir = os.path.join(baseDataset2Dir, 'cluster_%s' % cluster2)
                    curOutputDir = os.path.join(baseOutputDir, '%s_%s' % (title1, title2))
                    if not os.path.exists(curOutputDir): os.makedirs(curOutputDir)
                    FetchTwoDatasetSimi(motifLength, dataset1Dir, title1, dataset2Dir, title2, curOutputDir)
#end_func


def AnalyzeDiffClusterCoexist():
    """
    fetch same coexist from raw result from different cluster
    input: raw coexist result
    output: same coexist from different cluster
    """
    inputDir = r'C:\Users\Administrator\Desktop\coexist\coexist rlt\mir16'
    outputFn = r'C:\Users\Administrator\Desktop\coexist\coexist rlt\mir16\diffClusterCoexist.txt'
    patternPair2order = {}

    for localDir in os.listdir(inputDir):
        coexistFn = os.path.join(inputDir, localDir, 'patternPair2mutualInfo.txt')
        with open(coexistFn) as coexistFileobj:
            lines = coexistFileobj.readlines()
            patternPairCnt = len(lines)
            for lineIdx, line in enumerate(lines):
                line = line.strip()
                items = line.split('\t')
                p1 = items[0]
                p2 = items[2]
                patternPair2order.setdefault((p1, p2), [])
                patternPair2order[(p1, p2)].append((localDir, patternPairCnt, lineIdx))

    with open(outputFn, 'w') as outputFileobj:
        for (p1, p2), orderInfoLis in patternPair2order.iteritems():
            if len(orderInfoLis) > 1:
                orderInfoLis = [(items[0], '%s(in %s)' % (items[2], items[1])) for items in orderInfoLis]
                outputLine = '%s\t%s\t%s' % (p1, p2, orderInfoLis)
                outputFileobj.write('%s\n' % outputLine)
#end_func

def FormatwebAllCellLine():
    fnLis = ['evpedia_hsa_Caput_epithelial_cell_input', 'evpedia_hsa_Caput_luminal_fluid_input', 'evpedia_hsa_Cauda_epithelial_cell_input', 'evpedia_hsa_Cauda_luminal_fluid_input', 'evpedia_hsa_Colorectal_cancer_cell_SW480__input', 'evpedia_hsa_Lung_cancer_cell_DMS563__input', 'evpedia_hsa_Lung_cancer_cell_NCI_H69__input', 'evpedia_hsa_Seminal_plasma_input', 'exocarta_hsa_B_cells', 'exocarta_hsa_Colon_carcinoma_cells', 'exocarta_hsa_Colorectal_cancer_cells', 'exocarta_hsa_Endothelial_cells', 'exocarta_hsa_Plasma', 'exocarta_hsa_Serum', 'exocarta_hsa_T_cells', 'vesiclepedia_hsa_Breast_cancer_cells', 'vesiclepedia_hsa_B_cells', 'vesiclepedia_hsa_Lung_cancer_cells', 'vesiclepedia_hsa_Serum', 'vesiclepedia_hsa_T_cells', 'vesiclepedia_hsa_Urine']
    outputDic = {}
    for fn in fnLis:
        category, _, cellNm = fn.partition('_')
        outputDic.setdefault(category, [])
        outputDic[category].append(cellNm)
    print outputDic
#end_func

def FormatKEGGData():
    inputFn = r'C:\Users\Administrator\Desktop\2\kegg_path.txt'
    outputFn = r'C:\Users\Administrator\Desktop\2\kegg_path_rlt.txt'
    with open(inputFn) as inputFileobj, open(outputFn, 'w') as outputFileobj:
        lastTitle = ''
        lastSubtitle = ''
        nameLis = []
        for line in inputFileobj:
            line = line.rstrip()
            if not line: continue
            isTitle = len(line) - len(line.strip()) == 1
            isSubtitle = len(line) - len(line.strip()) == 3
            line = line.strip()
            if isTitle:
                lastTitle = line
            else:
                if isSubtitle:
                    if nameLis:
                        for curNm, subSubtitle in nameLis:
                            outputLine = '%s, %s, %s, %s' % (curNm, subSubtitle, lastSubtitle, lastTitle)
                            outputFileobj.write('%s\n' % outputLine)
                    nameLis = []
                    lastSubtitle = line
                else:
                    item0, _, restStr = line.partition(' ')
                    subSubtitle, _, pathId = restStr.strip(']').partition('[PATH:')
                    if not pathId: continue
                    subSubtitle = subSubtitle.strip()
                    nameLis.append((pathId, subSubtitle))
        for curNm, subSubtitle in nameLis:
            outputLine = '%s, %s, %s, %s' % (curNm, subSubtitle, lastSubtitle, lastTitle)
            outputFileobj.write('%s\n' % outputLine)
#end_func

def sendEmail():
    import smtplib
    import email.mime.text
    # my test mail
    mail_username = 'gtfish1987@gmail.com'
    mail_password = 'Gaomg2162319'
    from_addr = mail_username
    to_addrs = ('tgaochn@gmail.com')

    # HOST & PORT
    HOST = 'smtp.gmail.com'
    PORT = 25

    # Create SMTP Object
    smtp = smtplib.SMTP()
    print 'connecting ...'

    # show the debug log
    smtp.set_debuglevel(1)

    # connet
    try:
        print smtp.connect(HOST, PORT)
    except:
        print 'CONNECT ERROR ****'
        # gmail uses ssl
    smtp.starttls()
    # login with username & password
    try:
        print 'loginning ...'
        smtp.login(mail_username, mail_password)
    except:
        print 'LOGIN ERROR ****'
        # fill content with MIMEText's object
    msg = email.mime.text.MIMEText('Hi ,I am leehark')
    msg['From'] = from_addr
    msg['To'] = ';'.join(to_addrs)
    msg['Subject'] = 'hello , today is a special day'
    print msg.as_string()
    smtp.sendmail(from_addr, to_addrs, msg.as_string())
    smtp.quit()
#end_func

def FetchSeqCnt():
    fn = r'C:\Users\Administrator\Desktop\new 3.txt'
    allExpID = {15, 45, 82, 83, 84, 85, 86, 87, 88, 108, 109, 110, 111, 112, 115, 116, 117, 162, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 186, 203, 205, 218, 222, 223, 226, 228, 229, 230, 261, 265, }
    tmpSet = set()
    with open(fn) as f:
        f.readline()
        for line in f:
            items = line.strip().split('\t')
            miRID = items[0].lower()
            spec = items[1].lower()
            # expID = items[2] if len(items) == 3 else 15
            expID = items[2]
            if int(expID) not in allExpID:continue
            seqNm = (miRID, spec)
            tmpSet.add(seqNm)

    tmpSet = sorted(tmpSet, key=lambda x:x[0])
    tmpSet = ['%s\t%s' % (miRID, spec) for miRID, spec in tmpSet]
    print len(tmpSet)
    # Comm.showList(tmpSet)
#end_func

def FetchMIDistri():
    dataLis1 = [0.727825704, 0.727825704, 0.727825704, 0.727825704, 0.726415211, 0.726415211, 0.726415211, 0.721316496, 0.721316496, 0.721316496, 0.712753238, 0.706915995, 0.706915995, 0.706915995, 0.699832962, 0.699832962, 0.699832962, 0.699832962, 0.699832962, 0.691281037, 0.691281037, 0.691281037, 0.691281037, 0.680958806, 0.668452087, 0.653178417, 0.653178417, 0.653178417, 0.634293739, 0.634293739, 0.634293739, 0.610527204, 0.610527204, 0.579869232, 0.579869232, 0.579869232, 0.579869232, 0.538932567, 0.538932567, 0.538932567, 0.538932567, 0.538932567, 0.538932567, 0.538932567, 0.538932567, 0.503269341, 0.498371029, 0.498371029, 0.493523194, 0.481498508, 0.481498508, 0.481498508, 0.481498508, 0.481498508, 0.481498508, 0.481498508, 0.481498508, 0.480091737, 0.477835, 0.467305156, 0.457893173, 0.457893173, 0.439380054, 0.43064049, 0.421068615, 0.398963678, 0.396658952, 0.394718806, 0.394718806, 0.394718806, 0.372059744, 0.372059744, 0.372059744, 0.370497339, 0.353376107, 0.353376107, 0.353376107, 0.353376107, 0.346917484, 0.346917484, 0.342883086, 0.338912858, 0.338912858, 0.338912858, 0.338912858, 0.337600453, 0.337600453, 0.312507936, 0.297385289, 0.297385289, 0.297385289, 0.297385289, 0.297385289, 0.297385289, 0.297385289, 0.297385289, 0.296654216, 0.27640448, 0.267133439, 0.252857084, 0.252857084, 0.252857084, 0.252857084, 0.24741261, 0.24741261, 0.24741261, 0.24741261, 0.24741261, 0.24741261, 0.24741261, 0.244216948, 0.224533619, 0.219558022, 0.219558022, 0.219558022, 0.219558022, 0.218139157, 0.19023545, 0.19023545, 0.19023545, 0.19023545, 0.19023545, 0.19023545, 0.188795248, 0.174410947, 0.174410947, 0.174410947, 0.174410947, 0.161496257, 0.157438215, 0.152445331, 0.148391422, 0.138307904, 0.138307904, 0.136594328, 0.129499291, 0.128225173, 0.111129566, 0.111129566, 0.111129566, 0.111129566, 0.100926143, 0.09724841, 0.09724841, 0.09724841, 0.09724841, 0.09724841, 0.09724841, 0.09724841, 0.09724841, 0.09724841, 0.09724841, 0.086112948, 0.081761102, 0.081761102, 0.081761102, 0.081761102, 0.081761102, 0.061402251, 0.047700227, 0.047700227, 0.047700227, 0.047700227, 0.047700227, 0.047700227, 0.047700227, 0.046989042, 0.037134756, 0.037134756, 0.037134756, 0.037134756, 0.037134756, 0.037134756, 0.037134756, 0.034962368, 0.029046142, 0.028892268, 0.028079906, 0.022720604, 0.022720604, 0.022720604, 0.022720604, 0.021572784, 0.021572784, 0.019053147, 0.017695924, 0.017695924, 0.017695924, 0.017695924, 0.014938641, 0.014938641, 0.013660567, 0.013660567, 0.013660567, 0.01253319, 0.010397829, 0.008177459, 0.008177459, 0.007753212, 0.007753212, 0.007753212, 0.007753212, 0.007753212, 0.007118724, 0.006825674, 0.006323881, 0.005676703, 0.005614498, 0.004965171, 0.004965171, 0.004965171, 0.004965171, 0.004965171, 0.004965171, 0.004843738, 0.002788461, 0.002546153, 0.001584022, 0.001331702, 0.001199781, 0.001179359, 0.001179359, 0.000487602, 0.000487602, 0.000487602, 0.000388231, 0.000388231, 2.31E-05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    dataLis2 = [0.903912472, 0.903912472, 0.89802353, 0.897127128, 0.889530625, 0.885390009, 0.875590162, 0.851257178, 0.836058686, 0.836058686, 0.832970904, 0.82139741, 0.81442712, 0.811429134, 0.757842895, 0.716516067, 0.716463834, 0.71140344, 0.699632512, 0.686259204, 0.679745078, 0.676934142, 0.672183393, 0.670204145, 0.667705483, 0.659805306, 0.659805306, 0.659439813, 0.657297811, 0.654535791, 0.642899139, 0.641643262, 0.639127605, 0.638816, 0.6322536, 0.62850795, 0.621690099, 0.621646166, 0.616575381, 0.614280024, 0.609105497, 0.609105497, 0.607952018, 0.596635593, 0.591498527, 0.584209734, 0.578968576, 0.57056626, 0.555680243, 0.555431956, 0.550354303, 0.545466202, 0.534291383, 0.526218787, 0.524638009, 0.513199702, 0.50810948, 0.505079286, 0.498909999, 0.498222544, 0.462164087, 0.451236932, 0.443566146, 0.443385266, 0.442379576, 0.442077925, 0.435558071, 0.432016706, 0.424046147, 0.419792258, 0.413219867, 0.412528827, 0.407177857, 0.407176172, 0.386416832, 0.364639495, 0.364263852, 0.362593502, 0.362593502, 0.326296115, 0.296807981, 0.291485431, 0.280290133, 0.272933174, 0.254235678, 0.222133829, 0.168246757, 0.166418119, 0.156807289, 0.15240755, 0.146491159, 0.134837834, 0.123252009, 0.121399576, 0.120563941, 0.108489021, 0.107867416, 0.107140394, 0.095937774, 0.095874585, 0.094597243, 0.092699738, 0.088760181, 0.088524288, 0.086941082, 0.084380249, 0.083542211, 0.083063136, 0.081946839, 0.081946839, 0.078217419, 0.077730558, 0.074981343, 0.07356444, 0.07211326, 0.071322449, 0.067230345, 0.066229963, 0.065183738, 0.064057077, 0.06207164, 0.055080699, 0.052426556, 0.052132309, 0.050508534, 0.049743805, 0.047414877, 0.047414877, 0.046272755, 0.045212988, 0.043602841, 0.042287458, 0.042240505, 0.040757505, 0.039865334, 0.038247111, 0.038117268, 0.038117268, 0.036541865, 0.036250644, 0.034336854, 0.033039178, 0.032607719, 0.03153218, 0.031157738, 0.030749479, 0.03033531, 0.029943207, 0.028670151, 0.025699885, 0.024839205, 0.023164301, 0.021954724, 0.021161995, 0.020451149, 0.020069056, 0.020002729, 0.019913949, 0.019913949, 0.018174693, 0.017454756, 0.017264006, 0.016997496, 0.015996958, 0.014339308, 0.014326699, 0.014092858, 0.014092858, 0.014092858, 0.013870801, 0.013852761, 0.013852761, 0.01379089, 0.013238532, 0.012444624, 0.012389736, 0.012271114, 0.012249321, 0.01223055, 0.011705435, 0.011129913, 0.010371759, 0.010371759, 0.010000859, 0.00976893, 0.00966896, 0.009203982, 0.008352607, 0.008347742, 0.007706966, 0.007706966, 0.006921946, 0.006921946, 0.006604032, 0.006009635, 0.005550717, 0.005380741, 0.005358747, 0.004966401, 0.004559148, 0.004199709, 0.004087569, 0.004087569, 0.00401884, 0.003961995, 0.003814612, 0.003814612, 0.003512019, 0.002896513, 0.002699816, 0.0024785, 0.002380974, 0.002227305, 0.002182415, 0.00203885, 0.001768637, 0.001683253, 0.001532931, 0.001532931, 0.001513701, 0.001451587, 0.001451587, 0.001351885, 0.001264446, 0.001186517, 0.001186517, 0.001169422, 0.001071998, 0.001062621, 0.000983298, 0.000893748, 0.000716318, 0.000659315, 0.000592222, 0.000447026, 0.000380654, 0.000380654, 0.000380654, 0.00030995, 0.000294256, 0.000294256, 0.000253279, 0.000192325, 0.0001675, 0.000163618, 0.00011129, 9.56E-05, 8.40E-05, 8.15E-05, 7.81E-05, 6.98E-05, 6.34E-05, 3.97E-05, 2.12E-05, 2.12E-05, 1.45E-05, 8.36E-06, 8.36E-06, 1.68E-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    dataLis3 = [0.944229063, 0.943480686, 0.943101012, 0.942199976, 0.942081469, 0.93998476, 0.937159546, 0.935787495, 0.934551914, 0.933177244, 0.93243211, 0.927496555, 0.927496555, 0.926416464, 0.925851272, 0.914711216, 0.90941717, 0.908221583, 0.904304929, 0.899807496, 0.894587924, 0.888452849, 0.834825413, 0.827771414, 0.812303469, 0.812000379, 0.807316679, 0.800648823, 0.792471173, 0.787656955, 0.787191193, 0.786779272, 0.784005212, 0.766282095, 0.76081034, 0.757456385, 0.744300286, 0.744286953, 0.741966645, 0.740106925, 0.73905676, 0.738582402, 0.732181368, 0.730099204, 0.726912363, 0.713846904, 0.713846904, 0.713673815, 0.708815959, 0.708640932, 0.695847549, 0.690064367, 0.689554947, 0.684015805, 0.681843977, 0.680762959, 0.67748376, 0.676795553, 0.673875361, 0.671324892, 0.667219618, 0.665805463, 0.661892209, 0.661041684, 0.652275946, 0.647399603, 0.646419936, 0.645657231, 0.643890533, 0.643862852, 0.641171197, 0.630095165, 0.627755987, 0.627186424, 0.625890837, 0.619835714, 0.616828104, 0.594964539, 0.592630556, 0.590358815, 0.586372797, 0.585822135, 0.581245187, 0.578545945, 0.578454269, 0.572877108, 0.567932351, 0.545860355, 0.5457971, 0.544750491, 0.53760374, 0.537466305, 0.535784243, 0.533136886, 0.528353399, 0.526170001, 0.525611849, 0.525611849, 0.523575535, 0.519674068, 0.515796012, 0.512949389, 0.504938053, 0.49253258, 0.488282379, 0.486760577, 0.483248413, 0.481208471, 0.47217457, 0.464065361, 0.45605657, 0.455105256, 0.455105256, 0.454443, 0.44120044, 0.440089911, 0.439705001, 0.428076005, 0.422967488, 0.422953538, 0.422415864, 0.418119985, 0.414042142, 0.411273945, 0.404022238, 0.404022238, 0.402768177, 0.393113342, 0.391505763, 0.389817118, 0.38862, 0.381521737, 0.380654564, 0.377902069, 0.374927465, 0.374922684, 0.371078565, 0.37014891, 0.367654267, 0.364272218, 0.355166912, 0.352195168, 0.352101609, 0.343832707, 0.332782706, 0.332499142, 0.331352374, 0.32853938, 0.321868203, 0.315622044, 0.315512054, 0.310477171, 0.307241555, 0.30040231, 0.299129335, 0.297845602, 0.296145434, 0.295917958, 0.294877471, 0.292979543, 0.285130164, 0.283819051, 0.283819051, 0.280642325, 0.278964816, 0.278649658, 0.278647672, 0.272027958, 0.267169253, 0.265542901, 0.263591047, 0.262726318, 0.262479139, 0.260097449, 0.25725806, 0.249123298, 0.239516507, 0.23904382, 0.237148515, 0.236930567, 0.23568736, 0.235200907, 0.231519154, 0.228670525, 0.225788179, 0.225788179, 0.224373509, 0.216999074, 0.216999074, 0.215933484, 0.211368559, 0.211096802, 0.207069549, 0.207016849, 0.202566827, 0.19891734, 0.195209563, 0.193679222, 0.186974966, 0.185654862, 0.179999525, 0.173865476, 0.169389165, 0.169329705, 0.167829485, 0.166069774, 0.161700892, 0.154652878, 0.148666961, 0.142192664, 0.137552146, 0.122413646, 0.121115108, 0.116431864, 0.108597925, 0.105247751, 0.104992143, 0.104034522, 0.103134017, 0.099889602, 0.09729375, 0.090170906, 0.08396233, 0.074559525, 0.06275104, 0.058377919, 0.054585676, 0.04647512, 0.041705621, 0.039586797, 0.038290202, 0.037180118, 0.036293187, 0.031201462, 0.031157908, 0.029317476, 0.028774739, 0.028006391, 0.0258709, 0.025855941, 0.024562859, 0.023882343, 0.023568882, 0.023286416, 0.022878122, 0.02100021, 0.020950052, 0.019957557, 0.019272729, 0.018857392, 0.017690019, 0.017489379, 0.016164118, 0.015214537, 0.015173626, 0.014740415, 0.014104682, 0.014095551, 0.013686978, 0.01274818, 0.012147709, 0.011693661, 0.011420507, 0.009760427, 0.009521526, 0.009281978, 0.009272037, 0.008988433, 0.008988433, 0.008896761, 0.008826629, 0.008823679, 0.008636544, 0.007907434, 0.007841799, 0.007803955, 0.007623624, 0.007343521, 0.007334731, 0.007313471, 0.006739323, 0.00663047, 0.006503394, 0.0061014, 0.006073842, 0.005909982, 0.005285316, 0.005139564, 0.005139564, 0.005067713, 0.004893438, 0.004727981, 0.004727981, 0.004710022, 0.004694408, 0.004349103, 0.004275305, 0.003773251, 0.003661958, 0.003453733, 0.003318077, 0.003192679, 0.003192679, 0.003108921, 0.002912836, 0.002790155, 0.00272235, 0.002599307, 0.002505888, 0.002459329, 0.002436984, 0.002415383, 0.002375719, 0.002315796, 0.002315796, 0.002175374, 0.002105625, 0.002104693, 0.002072968, 0.002039572, 0.002039572, 0.002020298, 0.001999349, 0.001889662, 0.001730589, 0.001720548, 0.001632466, 0.001632466, 0.001565235, 0.001565235, 0.001482382, 0.001456489, 0.00138537, 0.001368996, 0.00125985, 0.001145847, 0.001134407, 0.001051084, 0.001029904, 0.000895973, 0.000886683, 0.000885188, 0.000879461, 0.000879461, 0.000873826, 0.000850671, 0.000837393, 0.000837393, 0.000824656, 0.000736483, 0.000736483, 0.000736483, 0.000716549, 0.000568309, 0.000548659, 0.000547581, 0.00049418, 0.000472463, 0.000466289, 0.000432179, 0.000431925, 0.000384836, 0.000296914, 0.000293845, 0.000293845, 0.000280083, 0.000266564, 0.000240561, 0.000220992, 0.000211789, 0.000184249, 0.000175277, 0.000157739, 0.000152938, 0.000152468, 0.000140861, 0.000140861, 9.23E-05, 8.67E-05, 5.58E-05, 5.58E-05, 5.35E-05, 5.27E-05, 4.74E-05, 3.70E-05, 3.69E-05, 3.69E-05, 3.64E-05, 2.90E-05, 2.82E-05, 2.39E-05, 2.37E-05, 9.80E-06, 3.69E-06, 2.06E-06, 6.46E-07, 3.13E-07, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    dataLis4 = [0.956301083, 0.954073006, 0.953882517, 0.951913517, 0.946355364, 0.944779411, 0.940561083, 0.939392621, 0.935810298, 0.934790471, 0.922259517, 0.914836303, 0.875808667, 0.872167207, 0.854146132, 0.842252679, 0.835283696, 0.808470085, 0.735640126, 0.722431865, 0.711206783, 0.690299815, 0.672921369, 0.669599791, 0.664329695, 0.663879801, 0.654427857, 0.653882295, 0.642561392, 0.631754543, 0.62551608, 0.624195384, 0.614455281, 0.610768796, 0.60737304, 0.590364115, 0.586708059, 0.574793627, 0.569046184, 0.561314154, 0.554562352, 0.537717118, 0.537717118, 0.537248443, 0.531633244, 0.524717998, 0.521887581, 0.518666215, 0.493167149, 0.491505176, 0.488968594, 0.480616419, 0.471471228, 0.469734811, 0.467122941, 0.457761136, 0.456675347, 0.453427568, 0.452202683, 0.442970925, 0.439119157, 0.433795604, 0.431587593, 0.416473002, 0.413551634, 0.394409724, 0.390635372, 0.389111634, 0.383919451, 0.382471152, 0.365534037, 0.360353452, 0.345476267, 0.340934998, 0.338615756, 0.336712786, 0.331827554, 0.324736414, 0.307008927, 0.288169672, 0.28350163, 0.28063034, 0.27981144, 0.275860462, 0.265792835, 0.264925731, 0.247439431, 0.237580234, 0.226646629, 0.211389708, 0.211381767, 0.210502613, 0.208052743, 0.204595801, 0.203811115, 0.179864287, 0.178733542, 0.176389173, 0.173295494, 0.168063568, 0.166215733, 0.162962115, 0.16202154, 0.161448583, 0.157134864, 0.151179844, 0.150523044, 0.147662622, 0.136797677, 0.132417033, 0.12936071, 0.12568684, 0.119019357, 0.113493013, 0.100762599, 0.092909823, 0.092024829, 0.08871204, 0.082155238, 0.07390684, 0.071241888, 0.068765267, 0.059535747, 0.058939595, 0.056680813, 0.054629464, 0.05368967, 0.04774099, 0.041869205, 0.041052995, 0.038359566, 0.038073464, 0.032371542, 0.03113917, 0.029676184, 0.028827314, 0.028561702, 0.025286782, 0.023790993, 0.022422568, 0.022396887, 0.01891546, 0.018910412, 0.018646607, 0.016498928, 0.016494399, 0.016174512, 0.015960851, 0.015389192, 0.014718127, 0.014473391, 0.014392371, 0.013549765, 0.013301458, 0.013063603, 0.012875005, 0.01280039, 0.012707509, 0.012571704, 0.012312008, 0.011614067, 0.011500362, 0.011148447, 0.011076368, 0.010131862, 0.009260477, 0.008062325, 0.006181801, 0.005742292, 0.005671618, 0.005424634, 0.004793402, 0.004793402, 0.004745146, 0.003118317, 0.001775734, 0.001707604, 0.001562969, 0.00141857, 0.001247352, 0.001025988, 0.000867825, 0.000630019, 0.000606549, 0.000542139, 0.000317456, 0.000203257, 0.000109629, 0.0000655, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    dataLis5 = [0.217602157, 0.167020894, 0.13683663, 0.102722446, 0.056085354, 0.056085354, 0.040919364, 0.026718501, 0.02374534, 0.021384269, 0.015753937, 0.012088521, 0.010838805, 0.010019793, 0.007848506, 0.005028917, 0.002847369, 0.001568658, 0.000326566, 0.000112746, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ]
    dataLis = [dataLis1, dataLis2, dataLis3, dataLis4, dataLis5]
    rltLis = []
    maxLen = max(map(lambda x:len(x), dataLis))
    for curDatalis in dataLis:
        boundryLis = range(25)
        boundryLis = map(lambda x:x*1.0/25, boundryLis)
        countLis, boundaryLis = Comm.CountAmountInRanges(curDatalis, boundryLis)
        rltLis.append(countLis)

    for i in xrange(maxLen):
        curLine = []
        for curDatalis in rltLis:
            if len(curDatalis) <= i:
                curStr = ''
            else:
                curStr = '%2f' % curDatalis[i]
            curLine.append(curStr)
        curLine = '\t'.join(curLine)
        print curLine
#end_func

def formatFasta():
    inputFn = r'C:\Users\Administrator\Desktop\new verifying data\all data.txt'
    outputFn = r'C:\Users\Administrator\Desktop\new verifying data\userInput.fa'

    with open(inputFn) as inputFileobj, open(outputFn, 'w') as outputFileobj:
        # seqLis = []
        # titleLis = []
        geneSeq = ''
        lastTitle = ''
        for line in inputFileobj:
            line = line.strip()
            if not line or line[0] == 'N' or line[0] == 'X' or line[0] == 'h': continue
            if line[0] == '>':  # title line
                if geneSeq and 'S' not in geneSeq and 'N' not in geneSeq:  # finish a line
                    outputFileobj.write('>%s\n' % lastTitle)
                    outputFileobj.write('%s\n' % geneSeq)
                    # seqLis.append(geneSeq)
                    # titleLis.append(lastTitle)
                geneSeq = ''
                lastTitle = line[1:]
            else:
                geneSeq += line

        if geneSeq and 'S' not in geneSeq and 'N' not in geneSeq:
            outputFileobj.write('>%s\n' % lastTitle)
            outputFileobj.write('%s\n' % geneSeq)
            #
    # return seqLis, len(seqLis), titleLis
#end_func


def runMEME():
    inputFnInfoLis = [ # (title, filename)
        ('nature', 'memeRlt/input_fa/NatComm.fa'),
        ('cell', 'memeRlt/input_fa/cell2016.fa'),
        ('sw620_seq39', 'memeRlt/input_fa/sw620-39.fa'),
        ('sw620_seq112', 'memeRlt/input_fa/sw620-112.fa'),
        ('cow_milk_miR', 'memeRlt/input_fa/cow-miR.fa'),
        ('cow_milk_RNA', 'memeRlt/input_fa/cow_milk_RNA.fa'),
        ('Par-Clip', 'memeRlt/input_fa/Par-Clip.fa'),
        ('miR-binding-1', 'memeRlt/input_fa/miR-binding-1.fa'),
        ('miR-binding-2', 'memeRlt/input_fa/miR-binding-2.fa'),
        ('miR-binding-3', 'memeRlt/input_fa/miR-binding-3.fa'),
        ('miR-binding-4', 'memeRlt/input_fa/miR-binding-4.fa'),
        ('miR-binding-5', 'memeRlt/input_fa/miR-binding-5.fa'),
    ]

    for title, inputFn in inputFnInfoLis:
        outputDir = 'memeRlt/output/MEME/%s' % title
        memeCmd = 'meme %s -rna -oc . -nostatus -time 18000 -maxsize 60000 -mod zoops -nmotifs 10 -minw 3 -maxw 50 -oc %s' % (inputFn, outputDir)
        print memeCmd
        os.system(memeCmd)

        outputDir = 'memeRlt/output/DREME/%s' % title
        # dremeCmd = 'dreme -v 1 -oc . -rna -p cell2016.fa -t 18000 -e 0.05 -dfile description'
        dremeCmd = 'dreme -v 1 -oc . -rna -p %s -t 18000 -e 0.05 -oc %s' % (inputFn, outputDir)
        print dremeCmd
        os.system(dremeCmd)
#end_func


def func():
    # FetchMilkRNAData()
    # GenerateNegTestData()
    # ExoName2seq()
    # SeqDiffExo()
    # formatJiangData()
    # FetchTargetPatternPvalue1()
    # FormatKEGGData()
    # FormatwebAllCellLine()

    # fetch co-exist pattern info
    # FetchAllClusterCoexist()
    # FetchOneDatasetCoexist()
    # AnalyzeDiffClusterCoexist()
    # FetchOneDatasetCoexist()

    # generate data of different cluster from raw data
    # sepBindingSiteCluster()

    # fetch similarity, output two file to visualize cluster: pattern2cov, pattern2simi
    # FetchDIffRNASimi()
    # FetchDiffClusterSimi()

    # test for sending emails
    # sendEmail()

    # FetchSeqCnt()

    # FetchMIDistri()
    # formatFasta()


    # FetchGivenPatternPvalue()
    FetchDiffClusterSimi_WithTomtom()
    # runMEME()
#end_test

def main():
    func()
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
