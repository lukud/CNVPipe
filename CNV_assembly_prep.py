#!/usr/bin/env python3

"""
a little wrapper to automate assembly
prep for WSSD calling with GEM/mrCaNaVar
"""


import sys,os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description='Prepare assembly for WSSD calling')
    parser.add_argument('-a', '--assembly', type=str, help='Path to assembly')
    parser.add_argument('-m', '--masks', type=str, nargs='+', help='Path to \
                        BED-file(s) with masking coordinates')
    parser.add_argument('-k,', '--ksize', type=int, help='Size of k-mers for masking',\
                        default=36, metavar='K')
    parser.add_argument('-s', '--step', type=int, help='Offset of k-mers for masking',\
                        default=36, metavar='S')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads to run GEM with',\
                        default=1)
    parser.add_argument('--gem', type=str, help='Alternative path to GEM-mapper if not in $PATH',\
                            default='', )
    args = parser.parse_args()

    return args

if __name__=='__main__':
    args=parse_arguments()
    print(args)
