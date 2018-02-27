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
                        BED-file(s) with masking coordinates. All files must be \
                        sorted by scaffold and position (sort -k1,1 -k2,2n)')
    parser.add_argument('-k,', '--ksize', type=int, help='Size of k-mers for masking',\
                        default=36, metavar='K')
    parser.add_argument('-s', '--step', type=int, help='Offset of k-mers for masking',\
                        default=5, metavar='S')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads to run GEM with',\
                        default=1)
    parser.add_argument('-o', '--outfolder', type=str, help='Output folder prefix,\
                        defaults to current directory', default='./')
    parser.add_argument('--gem', type=str, help='Alternative path to GEM-mapper if not in $PATH',\
                            default='', )
    args = parser.parse_args()

    return args

def load_assembly(fasta):
    """load a fasta assembly and returns a dict of scaffold -> sequence"""

    assembly={}

    with open(fasta_filename,"r") as fd:
        sequence = ""
        for line in fd:
            if line.startswith(">"):
                if len(sequence) > 0:
                        assembly[identifier]=list(sequence)
                identifier = line[1:].strip()
                sequence = ""
            else:
                    sequence+=line.strip()
        if len(sequence) > 0:
            assembly[identifier]=list(sequence)

    return assembly

def assembly_stats(assembly):
    pass

def mask_regions(*args, assembly, mask="N"):
    """
    parse a variable number of bed files and mask the correspoding sequences \
    in the asembly dict. I think this is horribly inefficient..but who cares
    """


    mask_count=0
    for bed in args:
        with open(bed, 'r') as b:
            for line in bed:
                scaffold, start, stop = line.split()
                start=int(start)
                stop=int(stop)
                mask_count+=stop-start
                assembly[scaffold][start:stop]=mask*(stop-start)

    return assembly, mask_count

def get_kmers(assembly, k, s, out_folder, out_prefix='assembly_kmers'):
    """
    create kmers for kmer-based maskgin fo the assemblies. Parameters k and s
    determine ksize and step
    """

    outfp="{}/{}_k_{}_s_{}.fasta".format(out_folder, out_prefix, k, s)

    with open(outfp) as kf:
        for scaffold in assembly:
            sequence=assembly[scaffold]
            for i in range(0,len(sequence), s):
                kmer=sequence[i:i+k]
                kf.write(">{}_{}_{}\n{}\n".format(scaffold, i, i+k, kmer))
    return outfp

def map_kmers(assembly, kmers, out_prefix, gem_path, gem_parameter_string=\
            '-q ignore --mismatch-alphabet "ACGT" -m 0.06 -e 0.1 -s 2 -d 1000 -D 1'):

    """
    generate index file and map kmers from the assembly
    """

    pass


if __name__=='__main__':
    args=parse_arguments()
    print(args)
