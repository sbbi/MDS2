# coding:utf8

"""
Author          :   Tian Gao
CreationDate    :   2016-9-20 20:57:55
Description     :
"""

import math
import os
import itertools
import re
import bitarray
import Comm
import random

COLORBLACK = "#000000"
COLORWHITE = "#FFFFFF"
COLORRED = "#FF0000"
COLORBLUE = "#0000FF"

def formatCovId(kmer2seqIdSet, totalSeqCnt):
    """
    format covered sequences Id set into binary num
    """
    userKmer2seqIdInt = {}
    for kmer, seqIdxSet in kmer2seqIdSet.iteritems():
        covIdxBitArray = IdxLis2bin(seqIdxSet, totalSeqCnt)
        covIdxBitInt = Comm.bitarray2int(covIdxBitArray)
        userKmer2seqIdInt[kmer] = covIdxBitInt
    return userKmer2seqIdInt
#end_func

def FetchCovWithRE(rePattern, seqLis):
    matchedSeqId = [seqId for seqId, seq in enumerate(seqLis) if re.findall(rePattern, seq)]
    cov = len(matchedSeqId)
    return cov
#end_func

def FetchSeqVis(seqLis, kmerSet):
    """
    Fetch all the visit of kmer in each sequence
    """
    kmer2visit = {}
    for kmer in kmerSet:
        kmer2visit[kmer] = set()
        for seqId, seq in enumerate(seqLis):
            kmerPos = seq.rfind(kmer)
            while kmerPos != -1:
                kmer2visit[kmer].add((seqId, kmerPos))
                seq = seq[:kmerPos]
                kmerPos = seq.rfind(kmer)
    return kmer2visit
#end_func

def FetchPatternCov(curPattern, seqLis, kmer2seqIdInt):
    """
    fetch the sequences ID in seqLis that curPattern covers
    ***
    to speed up 
    1. check whether current has already be calculated 
    2. parse pattern into kmer and store seq id in kmer2seqIdLis
    ***
    """
    seqCnt = len(seqLis)
    if curPattern in kmer2seqIdInt:
        patternCovSeqIdxInt = kmer2seqIdInt[curPattern]
    else:
        allCovIdLis = []
        kmerLis = FetchAllKmerFromPattern(curPattern)
        for kmer in kmerLis:
            if kmer in kmer2seqIdInt:
                allCovIdLis.append(kmer2seqIdInt[kmer])
            else:
                kmerCovIdxSet = {seqId for seqId, seq in enumerate(seqLis) if kmer in seq}
                kmerCovIdxBitarray = IdxLis2bin(kmerCovIdxSet, seqCnt)
                kmerCovIdxInt = Comm.bitarray2int(kmerCovIdxBitarray)
                kmer2seqIdInt[kmer] = kmerCovIdxInt
                allCovIdLis.append(kmerCovIdxInt)
        patternCovSeqIdxInt = reduce(lambda x, y: x | y, allCovIdLis)
        kmer2seqIdInt[curPattern] = patternCovSeqIdxInt
    cov = covIdxInt2covCnt(patternCovSeqIdxInt)
    return cov, patternCovSeqIdxInt, kmer2seqIdInt
#end_func

def IdxLis2bin(idxLis, seqCnt):
    """
    convert from coverage seq idx list(or set) to bitarray
    """
    bitLis = [False] * seqCnt
    for i in idxLis:
        bitLis[i] = True
    curBitarray = bitarray.bitarray(bitLis)
    return curBitarray
#end_func

def bitarray2covIdxSet(curBitarray):
    """
    convert from bitarray to coverage seq idx set
    """
    valueLis = bitarray.bitarray.tolist(curBitarray)
    covIdxSet = {covIdx for covIdx, value in enumerate(valueLis) if value}
    return covIdxSet
#end_func

def covIdxInt2covIdxSet(covIdxInt):
    """
    covert from binary coverage seq idx(form: int) to seq idx set
    """
    covIdxBitArray = Comm.int2bitarray(covIdxInt)
    covIdxSet = bitarray2covIdxSet(covIdxBitArray)
    return covIdxSet
#end_func

def covIdxInt2covCnt(covIdxInt):
    """
    covert from binary coverage seq idx(form: int) to coverage
    """
    covIdxBitArray = Comm.int2bitarray(covIdxInt)
    covCnt = bitarray.bitarray.count(covIdxBitArray)
    return int(covCnt)
#end_func

def PlotAlignment(alignLis, colorFormatLis):
    """
    plot a picture to demonstrate alignment with red color
    """
    im = Image.new("RGB", (330, len(alignLis) * 17 + 5), COLORWHITE)
    dr = ImageDraw.Draw(im)
    font = ImageFont.truetype(os.path.join("fonts", "msyh.ttf"), 14)

    for rowIdx, alignmentString in enumerate(alignLis):
        colorFormat = colorFormatLis[rowIdx]
        formatIdx = 0
        colorRange, curColor = colorFormat[formatIdx]
        for colIdx, curChar in enumerate(alignmentString):
            if colIdx > colorRange and formatIdx < len(colorFormat) - 1:
                formatIdx += 1
            colorRange, curColor = colorFormat[formatIdx]
            dr.text((2 + 8 * colIdx, 2 + 17 * rowIdx), curChar, font=font, fill=curColor)
    im.show()
# end_func

def GetSingleInfoContentByString(string, numOfPosibleString=4):
    """
    get information content of a string
    'ACCG', 4 -> 0.1505
    """
    cntDic = {}
    stringLen = len(string)
    IC = 0
    for letter in string:
        cntDic.setdefault(letter, 0)
        cntDic[letter] += 1
    for letter, cnt in cntDic.iteritems():
        prob = cnt * 1.0 / stringLen
        IC += prob * math.log(prob / (1.0 / numOfPosibleString), 10)
    return IC
# end_func

def GetLisInfoContent(seqLis):
    """
    get information content of several sequences with same length
    ['ACCG', 'ACCT'] -> [0.6020599913279623, 0.6020599913279623, 0.6020599913279623, 0.30102999566398114]
    """
    ICLis = []
    stringLen = len(seqLis[0])
    for idxY in range(stringLen):
        stringInPos = ''
        stringNum = len(seqLis)
        for idxX in range(stringNum):
            letter = seqLis[idxX][idxY]
            stringInPos += letter
        ICInPos = GetSingleInfoContentByString(stringInPos)
        ICLis.append(ICInPos)
    return ICLis
# end_func

def GetSingleInfoContentByLis(distriLis):
    """
    get information content of a list of integrals
    [1, 2, 1] -> 0.1505
    """
    totalSum = sum(distriLis)
    IC = 0
    for curCnt in distriLis:
        prob = curCnt * 1.0 / totalSum
        IC += prob * math.log(prob / 0.25, 10)
    return IC
# end_func

def FetchICFromPattern(patternStr, numOfPosibleString=4):
    """
    fetch information content from pattern
    A[TG]A[CG] -> [0.6020599913279623, 0.30102999566398114, 0.6020599913279623, 0.30102999566398114]
    """
    patternLis = parsePatternStr(patternStr)
    ICLis = [GetSingleInfoContentByString(singlePattern, numOfPosibleString) for singlePattern in patternLis]
    return sum(ICLis)
#end_func

def Seq2Kmer(seq, k):
    """
    divide the sequence by k-mer
    ('ATG', 2) -> ['AT', 'TG']
    :return:
    a list containing all k-mers in this sequence
    """
    kmerLis = []
    seqLen = len(seq)
    for startIdx in xrange(seqLen - k + 1):
        endIdx = startIdx + k
        kmer = seq[startIdx: endIdx]
        kmerLis.append(kmer)
    return kmerLis
# end_func

def Kmer2Str(kmerLis):
    """
    combine k-mer list into a string
    ['AT', 'TG'] -> 'ATG'
    :return:
    k-mer string
    """
    headKmer = kmerLis[0]
    charLis = [kmer[-1] for kmer in kmerLis[1:]]
    return headKmer + ''.join(charLis)
#end_func

def GenerateKmerLis(Base, k=2):
    """
    based on Base and k, generate kmer list
    eg.
    'ATCG', 2 -> ['AA', 'AT', 'AC', 'AG', 'TA', 'TT', 'TC', 'TG', 'CA', 'CT', 'CC', 'CG', 'GA', 'GT', 'GC', 'GG']
    """
    kmerLis = [''.join(item) for item in itertools.product(Base, repeat=k)]
    return kmerLis
#end_func

def GetSeqMatch(seq1, seq2):
    """
    get the match number of two sequences
    :return:
    match number
    """
    matchCnt = 0
    for idx, curChar1 in enumerate(seq1):
        curChar2 = seq2[idx]
        if curChar1 == curChar2: matchCnt += 1
    return matchCnt
#end_func

def FetchAllKmerFromPattern(patternStr):
    """
    parse a pattern and return all the relative kmers
    eg:
        A[CG][TC] -> ACT, ACC, AGT, AGC
    """
    rltLis = []
    parsedPatternLis = parsePatternStr(patternStr)
    for rawKmer in itertools.product(*parsedPatternLis):
        kmer = ''.join(rawKmer)
        rltLis.append(kmer)
    return rltLis
#end_func

def parsePatternStr(patternStr):
    """
    AT[CG]A[TC]A -> [A, T, CG, A, TC, A]
    :return:
    parsed pattern list
    """
    strLis = []
    while patternStr.find('[') != -1:
        head, _, rest = patternStr.partition('[')
        segPattern, _, patternStr = rest.partition(']')
        strLis.extend(list(head))
        strLis.append(segPattern)
    strLis.extend(list(patternStr))
    return strLis
#end_func

def parsePatternStrAdv(patternStr):
    """
    AT?[CG]?A[TC]A -> [A, T?, CG?, A, TC, A]
    """
    if '?' not in patternStr: return parsePatternStr(patternStr)

    # only handle with the string that include '?'
    parsedPatternLis = []
    subPatternLis = patternStr.split('?')
    for idx, subPattern in enumerate(subPatternLis[:-1]):
        parsedSubPatternLis = parsePatternStr(subPattern)
        parsedSubPatternLis[-1] += '?'
        parsedPatternLis += parsedSubPatternLis
    parsedPatternLis += parsePatternStr(subPatternLis[idx + 1])
    return parsedPatternLis
#end_func

def AssemblePatternLis(patternLis):
    """
    [A, T?, CG?, A, TC, A] -> 'AT?[CG]?A[TC]A'
    """
    formattedPatternLis = []
    for origPattern in patternLis:
        QFlag = True if origPattern[-1] == '?' else False
        origPattern = origPattern.rstrip('?')
        formattedPattern = origPattern if len(origPattern) == 1 else '[%s]' % origPattern
        if QFlag: formattedPattern += '?'
        formattedPatternLis.append(formattedPattern)
    assembledPattern = ''.join(formattedPatternLis)
    return assembledPattern
#end_func

def getMergeCnt(patternStr):
    """
    return how many merging a pattern does
    eg.
    'A[CT]A' -> 1
    'AT[CG]A[TC]A' -> 2
    'A[CTG]A' -> 2
    """
    strLis = parsePatternStr(patternStr)
    lengthLis = map(lambda x:len(x), strLis)
    return sum(lengthLis) - len(lengthLis)
#end_func

def getPatternLength(patternStr):
    """
    return the length of a pattern
    'AT[CG]A[TC]A' -> 6
    """
    return len(parsePatternStr(patternStr))
#end_func

def getPatternDist(patternStr1, patternStr2):
    """
    return similarity distance of two patterns
    (ATC, AGC) -> 1
    (A[TG]C, A[TCG]C) -> 0
    (A[TG]C, A[CA]C) -> 1
    """
    dist = 0
    strLis1 = parsePatternStr(patternStr1)
    strLis2 = parsePatternStr(patternStr2)
    for strIdx, curStr1 in enumerate(strLis1):
        curStr2 = strLis2[strIdx]
        if not isSinglePatternSimi(curStr1, curStr2): dist += 1
    return dist
#end_func

def isSinglePatternSimi(singlePatternStr1, singlePatternStr2):
    """
    return whether two pattern with only one digit similar to each other
    (A, A) -> True          (A, T) -> False
    (T, TC) -> True         (TC, TA) -> False
    (ATG, ACG) -> False
    """
    if len(singlePatternStr1) == 1 or len(singlePatternStr2) == 1: # most common case -> running faster
        if singlePatternStr1 in singlePatternStr2 or singlePatternStr2 in singlePatternStr1:
            return True
    else:
        charset1 = set(singlePatternStr1)
        charset2 = set(singlePatternStr2)
        if charset1.issubset(charset2) or charset2.issubset(charset1):
            return True
    return False
#end_func

def FetchDistriPercentage(distriLis, userCov):
    """
    given distribution and user coverage, return the percentage of that cover less 
    [10, 20, 50, 10, 6, 4], 4 -> 0.96 since 10 + 20 + 50 + 10 + 6 / sum = 96%
    """
    lessCovered = distriLis[:userCov + 1]
    return sum(lessCovered) * 1.0 / sum(distriLis)
#end_func

def FetchVarConfIdx(dataLis, percentage):
    """
    given a data list and a percentage, return the index of the first element before which sum of elements are 95%
    [10, 20, 50, 10, 6, 4] -> 4 since 10 + 20 + 50 + 10 + 6 / sum = 96%
    """
    # TODO: conbine the two function
    dataSum = sum(dataLis)
    accumulationPerc = 0.0
    for idx, data in enumerate(dataLis):
        curPerc = data * 1.0 / dataSum
        accumulationPerc += curPerc
        if accumulationPerc >= percentage:
            break
    return idx
#end_func

def FetchConf95Idx(dataLis):
    """
    given a data list, return the index of the first element before which sum of elements are 95%
    [10, 20, 50, 10, 6, 4] -> 4 since 10 + 20 + 50 + 10 + 6 / sum = 96%
    """
    dataSum = sum(dataLis)
    accumulationPerc = 0.0
    for idx, data in enumerate(dataLis):
        curPerc = data * 1.0 / dataSum
        accumulationPerc += curPerc
        if accumulationPerc >= 0.95:
            break
    return idx
#end_func

def GetKmerLis(base, k=2):
    """
    give a base (eg 'ATCG'), return all the k-mer in a list
    ('ATCG', 2) -> ['AT', 'AC', 'AG', ..., 'CG'] (16 elements in total)
    """
    kmerLis = []
    for item in itertools.product(base, repeat=k):
        kmer = ''.join(item)
        kmerLis.append(kmer)
    return kmerLis
#end_func

def loadSinglelineSeq(fastaFn):
    """
    load sequences from fasta file
    each sequence takes one line only
    eg. miRNA sequence
        >HSA-MIR-520E
        AAAGTGCTTCCTTTTTGAGGG
        >HSA-MIR-520B
        AAAGTGCTTCCTTTTAGAGGG
    :return:
    sequences list
    """
    seqLis = []
    minLen = float('inf')
    maxLen = 0

    with open(fastaFn) as fastaF:
        for line in fastaF:
            line = line.strip()
            if not line or line[0] == '>': continue
            line = line.replace('N', '')
            line = line.replace('Y', '')
            line = line.replace('S', '')
            seqLis.append(line)
            seqLen = len(line)
            if seqLen > maxLen: maxLen = seqLen
            if seqLen < minLen: minLen = seqLen
    return seqLis, len(seqLis), minLen, maxLen
# end_func

def loadMultilineSeq(fastaFn, minLength=0):
    """
    load fasta file and combine sequences in different lines together
    each sequence may take multiple lines
    eg. DNA sequence
        >ENSG00000000003|TSPAN6|ENST00000496771
        Sequence unavailable
        
        >ENSG00000000003|TSPAN6|ENST00000612152
        ATGCTAAAACTGTATGCAATGTTTCTGACTCTCGTTTTTTTGGTCGAACTGGTCGCTGCC
        ATCGTAGGATTTGTTTTCAGACATGAGATTAAGAACAGCTTTAAGAATAATTATGAGAAG
        GCTTTGAAGCAGTATAACTCTACAGGAGATTATAGAAGCCATGCAGTAGACAAGATCCAA
        AATACGTTGCATTGTTGTGGTGTCACCGATTATAGAGATTGGACAGATACTAATTATTAC
        TCAGAAAAAGGATTTCCTAAGAGTTGCTGTAAACTTGAAGATTGTACTCCACAGAGAGAT
        GCAGACAAAGTAAACAATGAACTGATTGGAATCTTTCTCGCCTACTGCCTCTCTCGTGCC
        ATAACAAATAACCAGTATGAGATAGTGTAA
    NOTICE: short lines are discarded
    """
    seqLis = []
    titleLis = []
    with open(fastaFn) as rndSeqFile:
        geneSeq = ''
        lastTitle = ''
        for line in rndSeqFile:
            line = line.strip()
            if not line: continue
            if line[0] == '>': # title line
                if geneSeq and len(geneSeq) >= minLength and 'S' not in geneSeq: # finish a line
                    seqLis.append(geneSeq)
                    titleLis.append(lastTitle)
                geneSeq = ''
                lastTitle = line[1:]
            else:
                geneSeq += line

        if len(geneSeq) >= minLength and 'S' not in geneSeq:
            seqLis.append(geneSeq)
            titleLis.append(lastTitle)

    return seqLis, len(seqLis), titleLis
#end_func

def loadOnlineInput(fastaFn):
    """
    load fasta file and combine sequences in different lines together
    each sequence may take multiple lines
    eg. DNA sequence
        >ENSG00000000003|TSPAN6|ENST00000496771
        Sequence unavailable

        >ENSG00000000003|TSPAN6|ENST00000612152
        ATGCTAAAACTGTATGCAATGTTTCTGACTCTCGTTTTTTTGGTCGAACTGGTCGCTGCC
        ATCGTAGGATTTGTTTTCAGACATGAGATTAAGAACAGCTTTAAGAATAATTATGAGAAG
        GCTTTGAAGCAGTATAACTCTACAGGAGATTATAGAAGCCATGCAGTAGACAAGATCCAA
        AATACGTTGCATTGTTGTGGTGTCACCGATTATAGAGATTGGACAGATACTAATTATTAC
        TCAGAAAAAGGATTTCCTAAGAGTTGCTGTAAACTTGAAGATTGTACTCCACAGAGAGAT
        GCAGACAAAGTAAACAATGAACTGATTGGAATCTTTCTCGCCTACTGCCTCTCTCGTGCC
        ATAACAAATAACCAGTATGAGATAGTGTAA
    NOTICE: short lines are discarded
    """
    seqLis = []
    titleLis = []
    minLen = float('inf')
    maxLen = 0
    with open(fastaFn) as fastaFileobj:
        geneSeq = ''
        lastTitle = ''
        for line in fastaFileobj:
            line = line.strip()
            if not line: continue
            if line[0] == '>':  # title line
                if geneSeq:  # finish a line, not the first line
                    seqLen = len(geneSeq)
                    if seqLen > maxLen: maxLen = seqLen
                    if seqLen < minLen: minLen = seqLen
                    seqLis.append(geneSeq.upper())
                    titleLis.append(lastTitle)
                geneSeq = ''
                lastTitle = line[1:]
            else:
                geneSeq += line

        seqLis.append(geneSeq.upper())
        titleLis.append(lastTitle)

    return seqLis, len(seqLis), titleLis, minLen, maxLen
# end_func

def loadSinglePPM(ppmFn):
    ppm = []
    with open(ppmFn) as ppmFileobj:
        alphabetLis = ppmFileobj.readline().strip().split('\t')[1:]
        for line in ppmFileobj:
            tmpDic = {}
            items = line.strip().split('\t')[1:]
            for idx, value in enumerate(items):
                key = alphabetLis[idx]
                tmpDic[key] = value
            ppm.append(tmpDic)
    return ppm
#end_func

def loadMultiPPM(patternLis, ppmDir):
    patter2ppm = {}
    for pattern in patternLis:
        ppmFn = 'PPM-%s.txt' % pattern
        ppmFn = os.path.join(ppmDir, ppmFn)
        ppm = loadSinglePPM(ppmFn)
        patter2ppm[pattern] = ppm
    return patter2ppm
#end_func

def FetchSinglePatternCovInSeqLis(pattern, seqLis):
    """
    Given a list of sequences(seqLis), fetch coverage of the pattern
    return the total coverage and which sequences the pattern covers
    """
    totalCoverage = 0
    seqSet = set()
    for seqId, seq in enumerate(seqLis):
        coverageCnt = 1 if re.search(pattern, seq) else 0
        totalCoverage += coverageCnt
        seqSet.add(seqId)
    return totalCoverage, seqSet
#end_func

def FetchSingleKmerCovInSeqLis(kmer, seqLis):
    """
    Given a list of sequences(seqLis), fetch coverage of the pattern
    return the total coverage and which sequences the pattern covers
    """
    totalCoverage = 0
    for seqId, seq in enumerate(seqLis):
        coverageCnt = 1 if kmer in seq else 0
        totalCoverage += coverageCnt
    return totalCoverage
#end_func

def FetchPatternCovInSeqLis(patternLis, seqLis):
    """
    Given a list of sequences(seqLis), fetch coverage of each pattern in patternLis
    """
    pattern2cov = {}
    pattern2seq = {}
    for curPattern in patternLis:
        totalCoverage, seqSet = FetchSinglePatternCovInSeqLis(curPattern, seqLis)
        pattern2cov[curPattern] = totalCoverage
        pattern2seq[curPattern] = seqSet
    return pattern2cov, pattern2seq
#end_func

def FetchPatternDistri(pattern, seqMatrix):
    """
    fetch the distribution list of a certain pattern in the sampled sequences matrix
    """
    distriLis = [0] * (len(seqMatrix[0]) + 1)
    for seqLis in seqMatrix:
        totalCoverage = 0
        for seq in seqLis:
            coverageCnt = 1 if re.search(pattern, seq) else 0
            totalCoverage += coverageCnt
        distriLis[totalCoverage] += 1
    return distriLis
#end_func

def FetchPvalueFromBG(bgDistri, curValue):
    # take only the larger ones and drop the current slot
    samplingFreq = sum(bgDistri)
    largerSampleCnt = sum(bgDistri[curValue:])
    pvalue = min(largerSampleCnt * 1.0 / samplingFreq, 1.0)
    return pvalue
#end_func

def FetchPvalueFromBG_TakeHalf(bgDistri, curValue):
    # take half of current slot
    samplingFreq = sum(bgDistri)
    largerSampleCnt = sum(bgDistri[curValue:]) + bgDistri[curValue - 1] / 2 if curValue - 1 >= 0 else sum(bgDistri[curValue:])
    pvalue = min(largerSampleCnt * 1.0 / samplingFreq, 1.0)
    return pvalue
#end_func

def FetchPvalueFromBG_FindNonZero(bgDistri, curValue):
    # this method is used in the 1st reviewed version
    # find the non-zero slot forward, add a random weight around 0.5 for each step
    samplingFreq = sum(bgDistri)
    largerSampleCnt = sum(bgDistri[curValue:])
    epsilon = (random.random() - 0.5) * 0.1
    curWeight = 0.5 + epsilon
    while largerSampleCnt == 0:
        curValue -= 1
        largerSampleCnt = bgDistri[curValue] * curWeight
        curWeight *= 0.5 + epsilon
    pvalue = min(largerSampleCnt * 1.0 / samplingFreq, 1.0)
    return pvalue
#end_func

def FetchCovInSeqLis(seqLis, kmer):
    """
    fetch seq idx covered by kmer
    """
    seqSet = set()
    for seqIdx, seq in enumerate(seqLis):
        if kmer in seq:seqSet.add(seqIdx)
    return seqSet
#end_func

def FetchCovInSeqLisMutliKmer(seqLis, kmerLis):
    """
    fetch seq idx covered by kmer
    """
    kmer2seqSet = {}
    for kmer in kmerLis:
        seqSet = FetchCovInSeqLis(seqLis, kmer)
        kmer2seqSet[kmer] = seqSet
    return kmer2seqSet
#end_func

def FetchVisCovInSeqLis(seqLis, kmerLis):
    """
    Given a list of sequences(seqLis), fetch visit count and coverage of each kmer in kmerLis
    """
    tmpDict = {}
    for kmer in kmerLis:
        totalVisitCnt = 0
        totalCoverage = 0
        for seq in seqLis:
            visitCnt = seq.count(kmer)
            coverageCnt = 1 if visitCnt > 0 else 0
            totalVisitCnt += visitCnt
            totalCoverage += coverageCnt
        tmpDict[kmer] = (totalVisitCnt, totalCoverage)
    return tmpDict
#end_func

def FetchCovInSeqMatrix(seqMatrix, kmer):
    """
    Given a matrix of sequences(from multiple sampling), fetch coverage of a kmer in kmerLis
    seqMatrix:
        a list of sequences
    kmerLis:
        a list of target kmer
    return:
        a dict that maps kmer to (visitCnt, Coverage)
    """
    tmpDict = {}
    for kmer in kmerLis:
        tmpDict.setdefault(kmer, [])
        for seqLis in seqMatrix:
            totalVisitCnt = 0
            totalCoverage = 0
            for seq in seqLis:
                visitCnt = seq.count(kmer)
                coverageCnt = 1 if visitCnt > 0 else 0
                totalVisitCnt += visitCnt
                totalCoverage += coverageCnt
            tmpDict[kmer].append((totalVisitCnt, totalCoverage))
    return tmpDict
# end_func

def FetchVisCovInSeqMatrix(seqMatrix, kmerLis):
    """
    Given a matrix of sequences(from multiple sampling), fetch visit count and coverage of each kmer in kmerLis
    seqMatrix:
        a list of sequences
    kmerLis:
        a list of target kmer
    return:
        a dict that maps kmer to (visitCnt, Coverage)
    """
    tmpDict = {}
    for kmer in kmerLis:
        tmpDict.setdefault(kmer, [])
        for seqLis in seqMatrix:
            totalVisitCnt = 0
            totalCoverage = 0
            for seq in seqLis:
                visitCnt = seq.count(kmer)
                coverageCnt = 1 if visitCnt > 0 else 0
                totalVisitCnt += visitCnt
                totalCoverage += coverageCnt
            tmpDict[kmer].append((totalVisitCnt, totalCoverage))
    return tmpDict
# end_func

def SingleMerge(char1, char2):
    """
    merge two single char into a pattern
    'A' + ' ' -> 'A?'
    'AG' + ' ' -> 'AG?'
    'AG' + 'AT' -> 'AGT'
    """
    if char1 == ' ':
        if char2[-1] != '?':
            return char2 + '?'
        else:
            return char2

    if char2 == ' ':
        if char1[-1] != '?':
            return char1 + '?'
        else:
            return char1

    if char1 in char2:
        return char2

    if char2 in char1:
        return char1

    return ''.join(set(char1) | set(char2))
#end_func

def SegmentMerge(segmentLis):
    """
    merge a list of segment with same length
    [ACC, ATG] -> 'A[CT][CG]'
    """
    patternLis = []
    segmentLength = len(segmentLis[0])
    for pos in xrange(segmentLength):
        charLis = map(lambda x:x[pos], segmentLis)
        patternPos = ''.join(set(charLis))
        patternPos = patternPos if len(patternPos) == 1 else '[%s]' % patternPos
        patternLis.append(patternPos)
    pattern = ''.join(patternLis)
    return pattern
#end_func

def MergePatternAndSegment(pattern, segment):
    """
    merge a pattern and a segment into a new pattern
    'AC[AT][CG]', 'AGTT' -> A[CG][AT][CTG]
    """
    patternLis = []
    segmentLength = len(segment)
    parsedLis = parsePatternStr(pattern)
    for pos in xrange(segmentLength):
        patternPos = set(parsedLis[pos])
        charPos = segment[pos]
        patternPos.add(charPos)
        mergedPatternPos = ''.join(patternPos)
        mergedPatternPos = mergedPatternPos if len(mergedPatternPos) == 1 else '[%s]' % mergedPatternPos
        patternLis.append(mergedPatternPos)
    pattern = ''.join(patternLis)
    return pattern
#end_func

def MergeIntoPattern(str1, str2):
    """
    merge two string into a pattern
    each char in the same position of the string will be merged
    the strings must have the same length
    """
    if len(str1) != len(str2): return ''
    singlePatternLis = []
    for idx, char1 in enumerate(str1):
        char2 = str2[idx]
        singlePattern = SingleMerge(char1, char2)
        singlePatternLis.append(singlePattern)
    pattern = AssemblePatternLis(singlePatternLis)
    return pattern
#end_func

def IsSimi(fragment1, fragment2):
    """
    whether the two fragments are similar or not
    simi method 2: same 2kmers count: n-3
    """
    segLen = len(fragment1)

    dimerLis1 = Seq2Kmer(fragment1, 2)
    dimerLis2 = Seq2Kmer(fragment2, 2)

    same2merCnt = 0
    # print dimerLis1, dimerLis2
    for dimerIdx, _ in enumerate(dimerLis1):
        if dimerLis1[dimerIdx] == dimerLis2[dimerIdx]:
            same2merCnt += 1
    if same2merCnt >= segLen - 3:
        return True
    else:
        return False
#end_func

def IsSimi1(fragment1, fragment2):
    """
    whether the two fragments are similar or not
    simi method 1: half of characters are the same and if the length is odd, one extra position should be the same
    """
    segLen = len(fragment1)
    commSegLen = segLen / 2

    # check common segment
    hasCommSeg = False
    dimerLis1 = Seq2Kmer(fragment1, commSegLen)
    dimerLis2 = Seq2Kmer(fragment2, commSegLen)
    for dimerIdx, _ in enumerate(dimerLis1):
        if dimerLis1[dimerIdx] == dimerLis2[dimerIdx]:
            hasCommSeg = True
            break

    if segLen % 2 == 1:
        fragment1 = '%s%s' % (fragment1[:dimerIdx], fragment1[dimerIdx + commSegLen:])
        fragment2 = '%s%s' % (fragment2[:dimerIdx], fragment2[dimerIdx + commSegLen:])

    # check separate common char when segment has odd char
    hasCommSepChar = True if (segLen % 2 == 0 or HasCommSepChar(fragment1, fragment2)) else False

    isSimi = True if hasCommSeg and hasCommSepChar else False
    return isSimi
#end_func

def HasCommSepChar(fragment1, fragment2):
    dimerLis1 = Seq2Kmer(fragment1, 1)
    dimerLis2 = Seq2Kmer(fragment2, 1)
    for dimerIdx, _ in enumerate(dimerLis1):
        if dimerLis1[dimerIdx] == dimerLis2[dimerIdx]: return True
    return False
#end_func

def func1():
    # seq1 = 'ATTCGG'
    # seq2 = 'CCTTAG'
    # patternStr = 'AT[CG]AA[AT]TC'
    #print Seq2Kmer('ATCGG', 3)
    #print GetSeqMatch(seq1, seq2)
    # print getPatternLength(patternStr)
    # print parsePatternStr(patternStr)
    # patternStr1 = 'A[GT][ATCG]'
    # patternStr2 = 'A[CT]'
    # print getPatternDist(patternStr1, patternStr2)
    # print getMergeCnt(patternStr1)
    # print FetchConf95Idx([10, 20, 50, 10, 4, 6])
    # print parsePatternStrAdv('  A?T?[CG]?A[TC]A?')
    # print SingleMerge('TG', 'TA')
    # print AssemblePatternLis(['A', 'T?', 'CG?', 'A', 'TC', 'A'])
    # print MergeIntoPattern(seq1, seq2)
    # print SegmentMerge([seq1, seq2])
    # print GetSingleInfoContentByString('AA', 4)
    # FetchICFromPattern('A[TG]A[CG]')
    # print MergePatternAndSegment('AC[AT][CG]', 'AGTT')
    # print FetchPvalueFromBG([1,2,3,4,5], 4)
    # print FetchDistriPercentage([10, 20, 50, 10, 6, 4], 4)

    # refFn = r'D:\project\MotifFinding\sources\old\human_CDS.fasta'
    # refFn = r'C:\Users\Administrator\Desktop\ref.txt'
    # refSeqLis = LoadGeneSeq(refFn)
    # for seq in refSeqLis:
    #     print seq

    # print FetchAllKmerFromPattern('[ABC][123]')
    # print IdxLis2bi([1,2,4], 10)
    # a = covIdxInt2covIdxSet(1077521343739520)
    # print map(lambda x:x+1, a)
    # print len(a)

    # a = set([0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 26, 28, 29, 30, 31, 32, 34, 36, 37, 38, 39, 40, 41, 42, 43, 45, 46, 47, 48, 49, 50, 51, 52, 53, 55, 56, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 92, 93, 95, 96, 98, 99, 100, 101, 102, 103, 104, 106, 107, 108, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 121, 122, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 154, 155, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 179, 180, 181, 182, 183, 184, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 199, 200, 201, 202, 203, 204, 205, 206, 208, 210, 211, 212, 213, 214, 215, 216, 217, 219, 221, 222, 223, 224, 225, 226, 227, 228, 229, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 301, 302, 303, 304, 305, 306, 307, 308, 310, 311, 312, 313, 314, 315, 316, 318, 319, 320, 321, 322, 323, 324, 326, 327, 328, 329, 330, 331, 332, 334, 335, 336, 337, 338, 340, 341, 343, 344, 345, 346, 347, 348, 349, 350, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 390, 391, 392, 394, 395, 397, 398, 399, 400, 401, 402, 403, 404, 406, 407, 408, 409, 410, 411, 412, 413, 415, 416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428, 429, 430, 432, 433, 434, 435, 436, 437, 439, 440, 441, 442, 443, 444, 445, 447, 448, 449, 450, 451, 452, 455, 456, 457, 458, 459, 461, 462, 463, 465, 466, 467, 468, 469, 470, 471, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 487, 488, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499, 500, 501, 502, 503, 504, 505, 506, 507, 508, 509, 511, 513, 514, 515, 516, 517, 518, 520, 522, 523, 524, 525, 526, 527, 528, 529, 530, 532, 533, 534, 535, 536, 538, 539, 540, 541, 542, 543, 544])
    # IdxLis2bi(a, 114)
    # print IsSimi('ACGAA', 'ACTAA')
    # print IsSimi('ACGAA', 'ACGTT')
    # print IsSimi('ACAAA', 'ACGAT')
    # print parsePatternStr('AAAC')

    # rePattern = '[AG]G[AG]\w{2,5}GU'
    # seqLis = ['AGGCAAGCCAGGACAAAGUGUG', 'UUUAAUAUCCUUGAGCCUGGGCAAGUGCACAAGU', 'AUGGACUUGCCAUGACAGGUGCAGUUC', 'AAGACCGGAGGGAUGAGUGCUU', 'UUCCCGGUAGGGCGAGUGCAU', 'ACCCCGGCGGGGACACGGUGGCAAU', 'AGGGACUGGUGCAAAUGGGUCCAGAAUUUUCAAAUCGAAUGCUCUGUGUU', 'AAGGGACAGCUGCAGUAGCU', 'AAAAGCCACUGGGGACGAGACAGGUGCUAAAGU', 'AGGGGAUAUUUGU', 'AGGGACCUGUGCAGUGGGCUC', 'AGGGGCUGGUGCAAAGGAU', 'AGAGGAAGGGCAAGUGUGCUGU', 'AUGAGAAAUCAAGGGAUUAGUGCAACCAGU', 'UGACCAGGACAAACGUGCAAUAAUGCC', 'CCCUCAAGAAGGAGGGCAGGUGU', 'UGGACAAGUGCACUGAACUA', 'UGGAGGACUGGUGCAAUCAUCC', 'UGAUGUAGAGAAACGCUCCAGAGGACAAGUGCUGUUUGAU', 'CUCCAAUUUCCUGUAGGACGAGUGCACCGC', 'UUCAAGAGCGUGUGCAGGGCAAGUGC', 'AAAAGAGGGACAAGUGGCUGG', 'GCUCAUUUACCAGGGACGAGUUCUGCAAGAUGAU', 'CAGGGACAUUUGAAGUGCAAAUU', 'AUUUUAAGGGACAGGUGAAUUU', 'GUAACUUUAAGAGGGCAUUGUGCAAUAGUU', 'AAGAAGGGACAAGGUGCUU', 'UUCGAGGGACAUGUGUCAGC', 'UCCACAGGGACGGUGCGCUCAC', 'AUGCCAGUGGGCAGUGCAUGUGGAAAGU', 'UAUUGCCUGGGCAGUGCCC', 'GGACAUGAGGGGACAGUGCUCAAUAAC', 'CUGUGUGGGGACGGUGCUGGCCAGCAGAC', 'AAAUGGGAGACCGAGCGAGCGCGGCAAGUGCUGGAAC', 'GACACCACCCAGGGACAGUGCCUAUGU', 'AGGGACUGUGCCAGAAAAAC', 'AGGGACGUGCAACAUACAGCUU', 'AGGGACAAGUAGGUGCCUUCGGU', 'UACGGGGGCAAUUGCUGCAAUGC', 'UUAGACAAGCCAGUGUGGGUGCAGGAAUUC', 'GUGGACAAGAUGCAGCUGCUGGAGAUU', 'GCAGGGCAGUGCACCCUG', 'UGGGGGAGUGCAUCAUCGCU']
    # print FetchCovWithRE(rePattern, seqLis)

    # bgDistri = [10,8,5,0,0,0,0]
    # curValue = 3
    # for curValue in range(8):
    #     print FetchPvalueFromBG_FindNonZero(bgDistri, curValue)

    # print parsePatternStr('ACTG')
    # motifLength = 2
    # initMinIC = (0.61 * 0.25 + 0.3 * 0.5 + 0.125 * 0.25) * motifLength
    # initMaxIC = 0.61 * motifLength  # 0 variables in each position
    # print initMinIC, initMaxIC, FetchICFromPattern('[AC]G')

    # MEME
    # curStr = 'A[AG][AG][ACG][AUC]G[GCU]U[CGU]U[ACGU][CU]C[AUG]UU[AU][ACU]' # nature
    # curStr = '[AU][AG]GC[ACU]C[AU][CG]' # cell
    # curStr = '[AU][AG]C[CGU][AG][UG][CUG]U[UG][CUG]A[CG][AU]GG[AG]' # sw620

    # COSMO
    # curStr = '[ACGU][CU]C[AU]' # nature
    # curStr = '[ACGU][ACU]C[ACU]' # cell
    # curStr = '[ACU][ACU][ACU][ACU]' # sw620

    # Improbizer
    # curStr = '[AG][AU][AC][CG]' # nature
    # curStr = '[ACG][ACG]U[GU]' # cell
    # curStr = '[AUG]A[UG][UCG]' # sw620

    # DMINDA2
    # curStr = '[AT]C[ACT][AG]GG' # nature
    # curStr = 'GG[CG][ACT][CT][CT]' # cell
    # curStr = 'AGC[ACT]C[CT]' # sw620

    # MERCI
    # curStr = '[ATG][AG][AG][ACTG][ACG]' # nature
    # curStr = '[AT][AC][TG][AT]' # cell
    curStr = '[ACTG][ACTG][CTG][ACTG][ACTG]' # sw620


    print FetchICFromPattern(curStr)

    pass
# end_test

def main():
    func1()
# end_main

if __name__ == "__main__":
    main()
# end_if
