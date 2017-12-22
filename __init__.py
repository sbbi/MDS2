#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ProjectName     : 
Author          : Tian Gao
CreationDate    : 2017/5/27
Description     :

"""
import sys
import pdb
import BioinfoComm, Comm
from Comm import print0

ENABLE_CPROFILE = False
PRINT_INTO_LOG = False

LOG_FN = r'C:\Users\Administrator\Desktop\log.txt'

class X:
    def __init__(self):
        pass
    #end_func

    def fun2(self):
        pass
    #end_func
#end_class

def fun1():
    print 123
#end_func

def test():
    fun1()
    
    x = X()
    x.fun2()
#end_test

def main():
    if PRINT_INTO_LOG:
        sys.stdout = open(LOG_FN, 'w')
        test()
        sys.stdout.close()
    else:
        test()
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
