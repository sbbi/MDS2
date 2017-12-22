# coding:utf8

"""
Author          :   Tian Gao
CreationDate    :   2016-10-14 10:52:02
Description     :
"""

import datetime
import pickle
import os
import bitarray
import logging

INFO = logging.getLogger('root').info

def bitarray2int(x):
    """
    convert from bitarray to int 
    """
    return int(bitarray.bitarray.to01(x), 2)
#end_func

def int2bitarray(x):
    """
    convert from int to bitarray
    """
    return bitarray.bitarray(bin(x)[2:])
#end_func

def saveSingleDict(curDict, filenm):
    with open(filenm, 'w') as curFile:
        pickle.dump(curDict, curFile)
#end_def

def loadSingleDict(filenm):
    with open(filenm) as curFile:
        curDict = pickle.load(curFile)
    return curDict
#end_def

def print0(isLogging=False):
    if not isLogging:
        print '=' * 40
    else:
        INFO('=' * 40)
# end_func

def print9(isLogging=False):
    if not isLogging:
        print '==' * 40
    else:
        INFO('==' * 40)
# end_func

def PrintTime():
    now = str(datetime.datetime.now())
    print now[:19]
# end_func

def PrintWithTime(obj, isLogging=False):
    if not isLogging:
        now = str(datetime.datetime.now())
        print "%s - %s" % (now[:19], obj)
    else:
        INFO(obj)
#end_func

def showDict(curDict, isLogging=False):
    if not isLogging:
        for ker, value in curDict.iteritems():
            print "%s : %s" % (ker, value)
    else:
        for ker, value in curDict.iteritems():
            INFO("%s : %s" % (ker, value))
#end_func

def showList(curLis, isLogging=False):
    if not isLogging:
        for item in curLis:
            print item
    else:
        for item in curLis:
            INFO(item)
#end_func

def ShowMatrix(dataMatrix):
    for dataLis in dataMatrix:
        dataStr = '\t'.join(map(lambda x:str(x), dataLis))
        print dataStr
#end_func

def FetchLisMinOffset(lis1, lis2):
    if not lis1 or not lis2: return -1
    minOffset = float('inf')
    for val1 in lis1:
        for val2 in lis2:
            if val2 <= val1: continue
            curOffset = val2 - val1
            minOffset = min(curOffset, minOffset)
    if minOffset == float('inf'): minOffset = -1
    return minOffset
#end_func

def CountAmountInRanges(numLis, boundaryLis=None, isWithLabel=False):
    """
    count how many numbers located in each range(bounded by boundaryLis)
    ([1.1] * 4, [0, 1, 2]) -> [0, 0, 4, 0] since 1 < 1.1 < 2
    """
    if not boundaryLis:
        maxNum = int(max(numLis))
        minNum = int(min(numLis))
        deltaNum = int(max((maxNum - minNum) * 1.0 / 10, 1))
        boundaryLis = range(minNum, maxNum, deltaNum)
    boundaryLis.append(boundaryLis[-1] + boundaryLis[1] - boundaryLis[0])
    boundaryLis.append(boundaryLis[-1] + boundaryLis[1] - boundaryLis[0])
    countLis = [0] * (len(boundaryLis))
    for num in numLis:
        rangeIdx = GetRangeIdx(num, boundaryLis)
        # print num, rangeIdx, boundaryLis
        countLis[rangeIdx] += 1
    return countLis, boundaryLis
#end_func

def GetRangeIdx(num, boundaryLis):
    """
    get the range index of num in boundaryLis
    (-1, [0, 1, 2]) -> 0
    (1.1, [0, 1, 2]) -> 2
    """
    for rangeIdx, boundary in enumerate(boundaryLis):
        if num < boundary:
            return rangeIdx
    return rangeIdx + 1
#end_func

def GetBoundryLis(start, end, interval, epislon):
    """
    :return:
     [('-inf', start - epislon), (start - epislon, start - epislon + interval), (start - epislon + interval, start - epislon + 2 * interval), ..., (end - epislon, 'inf')]
     eg:
     GetBoundryLis(30, 90, 10, 0.1) ->
     [(-inf, 29.9), (29.9, 39.9), (39.9, 49.9), (49.9, 59.9), (59.9, 69.9), (69.9, 79.9), (79.9, 89.9), (89.9, inf)]
    """
    boundaryLis = []
    curFloor = - float('inf')
    curCeil = start - epislon
    while curFloor < end:
        boundaryLis.append((curFloor, curCeil))
        curFloor = curCeil
        curCeil = curCeil + interval
    boundaryLis.append((curFloor, float('inf')))
    return boundaryLis
#end_func

def GetMinAndMax(valueLis, indexLis):
    """
    get min and max value of valueLis, item in this list must be tuple or list
    eg:
    GetMinAndMax([(29.9, 39.9), (39.9, 49.9), (49.9, 59.9), (59.9, 69.9), (69.9, 79.9)], [0, 1]) ->
    [(29.9, 69.9), (39.9, 79.9)]
    """
    minAndMaxValueLis = [(float('inf'), - float('inf'))] * len(indexLis)
    for item in valueLis:
        for i in indexLis:
            value = item[i]
            minValue = minAndMaxValueLis[i][0]
            maxValue = minAndMaxValueLis[i][1]
            minValue = value if value < minValue else minValue
            maxValue = value if value > maxValue else maxValue
            minAndMaxValueLis[i] = (minValue, maxValue)
    return minAndMaxValueLis
#end_func

def FormatIdx2DTo1D(x, y, xRange):
    idx1D = x * xRange + y
    return idx1D
#end_func

def FormatIdx1DTo2D(n, xRange):
    idx2D = (n / xRange, n % xRange)
    return idx2D
#end_func

def floatRange(start, end, step):
    rltLis = []
    curValue = start
    while curValue < end:
        rltLis.append(curValue)
        curValue += step
    return rltLis
#end_func

def CntFreq(lis):
    """
    count the frequency of each item in the inputted list
    [1, 1, 1, 2, 3, 3] -> {1:3, 2:1, 3:2}
    :return:
     frequency dict
    """
    freqDict = {}
    for item in lis:
        freqDict.setdefault(item, 0)
        freqDict[item] += 1
    return freqDict
#end_func

def GroupLisByValue(lis):
    """
    count the frequency of each item in the inputted list
    [1, 1, 1, 2, 3, 3] -> {1:[0, 1, 2], 2:[3], 3:[4, 5]}
    :return:
     group dict
    """
    groupDic = {}
    for idx, item in enumerate(lis):
        groupDic.setdefault(item, [])
        groupDic[item].append(idx)
    return groupDic
#end_func

def UniqCnt(lis):
    """
    count how many elements are unique in a list
    """
    return len(set(lis))
#end_func

def TransposeLis(strLis):
    posStrLis = []
    strLen = len(strLis[0])
    for pos in range(strLen):
        posStrLis.append(''.join([curStr[pos] for curStr in strLis]))
    return posStrLis
#end_func

def FetchPathInLis(dataLis, length=2):
    pathLis = []
    maxStartIdx = len(dataLis) - length + 1
    for startIdx in range(maxStartIdx):
        pathLis.append(dataLis[startIdx: startIdx + length])
    return pathLis
#end_func

def FetchAllPath(dataLis):
    allPathLis = []
    for length in range(len(dataLis))[1:]:
        allPathLis.append(FetchPathInLis(dataLis, length))
    return allPathLis
#end_func

def GetLCM(x, y):
    greater = max(x, y)
    while True:
        if (greater % x == 0) and (greater % y == 0):
            return greater
        greater += 1
#end_func

def GetLisLCM(lis):
    curLCM = lis[0]
    for curNum in lis[1:]:
        curLCM = GetLCM(curLCM, curNum)
    return curLCM
# end_func

def CreateFolder(dirLis):
    """
    create a list of folder
    :param dirLis: 
    :return: 
    """
    for curDir in dirLis:
        if not os.path.exists(curDir): os.makedirs(curDir)
#end_func

def FetchAllComb(strLis):
    pass
#end_func

def sendEmail(toAddrs, title, content):
    import smtplib
    import email.mime.text
    # my test mail
    mail_username = 'gtfish1987@gmail.com'
    mail_password = 'Gaomg2162319'
    from_addr = mail_username
    toAddrs = ['tgaochn@gmail.com', ]

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
    msg = email.mime.text.MIMEText(content)
    msg['From'] = from_addr
    msg['To'] = ';'.join(toAddrs)
    msg['Subject'] = title
    print msg.as_string()
    smtp.sendmail(from_addr, toAddrs, msg.as_string())
    smtp.quit()
#end_func

def rawFormat(fn):
    # format filename so that it can be processed in cmd
    fn = fn.replace('[', '\[').replace(']', '\]').replace(' ', '\ ')
    return fn
#end_func

def test():
    # print GetRangeIdx(3.1, [0,1,2])
    # print CountAmountInRanges([1.1] * 4, [100, 200, 300])

    # print GetBoundryLis(30, 90, 10, 0.1)
    # print GetMinAndMax([(29.9, 39.9), (39.9, 49.9), (49.9, 59.9), (59.9, 69.9), (69.9, 79.9)], [0, 1])
    # a = {'julius1':{'cellphone':'13800000000','tel':'0512-34343534','qq':'354564656'},
    #            'julius2':{'cellphone':'13300000000','tel':'0513-34343534','qq':'454564656'}}
    # filenm = r'D:\project\test\PyTest\1.txt'
    # saveSingleDict(a, filenm)
    #
    # b = loadSingleDict(filenm)
    # print b
    # print GroupByValue([1, 1, 1, 2, 3, 3])
    # FetchPathInLis([1,2,3,4])
    # print FetchAllPath([1,2,3,4])
    # PrintWithTime('123')
    # print GetLisLCM([1,2,3])
    # CreateFolder([r'C:\Users\Administrator\Desktop\3'])
    # print int2bitarray(8) | bitarray.bitarray('0001')
    # print bitarray2int(bitarray.bitarray('1000')) | 1
    print FetchLisMinOffset([5], [4,10])
# end_test

def main():
    test()
# end_main

if __name__ == "__main__":
    main()
# end_if