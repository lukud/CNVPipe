#!/usr/bin/env python3

"""
a little wrapper to automate assembly
prep for WSSD calling with GEM/mrCaNaVar
"""


import sys,os
import argparse
import subprocess
import logging
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Prepare assembly for WSSD calling',\
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser.add_argument('-o', '--out_folder', type=str, help='Output folder prefix,\
                        defaults to current directory', default='CNV_assembly_prep')
    parser.add_argument('--gem', type=str, help='Alternative path to GEM-mapper if not in $PATH',\
                            default='', )
    args = parser.parse_args()

    return args

def load_assembly(fasta_filename):
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

def write_assembly(assembly, filepath):
    """
    Writes a dictionary of name->seq as a fasta file to disk
    """
    with open(filepath, 'w') as f:
        for seq in assembly:
            f.write(">{}\n{}\n".format(seq, "".join(assembly[seq])))

def assembly_stats(assembly):
    pass

def mask_regions(*args, assembly, mask="N"):
    """
    parse a variable number of bed files and mask the correspoding sequences \
    in the asembly dict. I think this is horribly inefficient..but who cares
    """

    mask_count=0
    beds=args[0]

    for bed in beds:
        with open(bed, 'r') as b:
            for line in b:
                scaffold, start, stop = line.split("\t")
                start=int(start)
                stop=int(stop)
                mask_count+=stop-start
                assembly[scaffold][start:stop]=mask*(stop-start)

    return assembly, mask_count

def get_kmers(assembly, k, s, out_folder, out_prefix='asm_kmers'):
    """
    create kmers for kmer-based maskgin fo the assemblies. Parameters k and s
    determine ksize and step
    """

    outfp="{}/{}_k_{}_s_{}.fasta".format(out_folder, out_prefix, k, s)

    with open(outfp, 'w') as kf:
        for scaffold in assembly:
            sequence=assembly[scaffold]
            for i in range(0,len(sequence), s):
                kmer=sequence[i:i+k]
                kmerstr="".join(kmer)
                if len(kmerstr) == k:
                    kf.write(">{}:{}-{}\n{}\n".format(scaffold, i, i+k, kmerstr))
    return outfp

def map_kmers(assembly, kmers, threads, out_prefix, \
            out_folder, gem_path, gem_parameter_string=\
            '-q ignore --mismatch-alphabet "ACGT" -m 0.06 -e 0.1 -s 2 -d 1000 -D 1'):

    """
    generate index file and map kmers from the assembly
    """
    #create the index
    out="{}/{}".format(out_folder, out_prefix)

    gem_index_cmd="{} gem-indexer -t {} -i {} -o {}".format(gem_path, threads, \
                                                    assembly, out)

    print(gem_index_cmd)

    try:
        p = subprocess.Popen([gem_index_cmd], stdout=subprocess.PIPE, \
                            stderr=subprocess.PIPE, shell=True)
        stdout, stderr = p.communicate()
    except:
        logging.critical("Couldn't create the gem-index for basemasking...abort")
        sys.exit(1)

    logging.info('\nGEM-INDEX LOG:\n')
    print(stdout, stderr)
    logging.info('\nGEM-INDEX LOG END\n')

    gem_mapping_cmd = "{} gem-mapper -T {} -I {} -i {} -o {} {}".format(gem_path, threads, \
                                                                    out+'.gem', kmers, out, gem_parameter_string)
    print(gem_mapping_cmd)
    try:
        p = subprocess.Popen([gem_mapping_cmd], stdout=subprocess.PIPE, \
                            stderr=subprocess.PIPE, shell=True)
        stdout, stdserr = p.communicate()
    except:
        logging.critical("Couldn't map kmers with gem...abort")
        sys.exit(1)

    logging.info('\nGEM-MAPPING LOG:\n')
    print(stdout, stderr)
    logging.info('\nGEM-MAPPING LOG END\n')

    return out+".map"

def count_mask_kmer_mappings(map_file, assembly, out_prefix, max_placements=20,\
                            N_threshold=10, mask="N",
                            write_overrepresented_only=True):
    """
    read in the map_file of self mappings, mask out all kmers with more then
    max_placements mappings
    """

    mappings=defaultdict(int)

    with open(map_file, 'r') as m:
        for line in m:
            read_id, sequence, trash, maps=line.rstrip().split('\t')
            #discard mappings with too many N's from the assembly
            if sequence.upper().count('N')>N_threshold:
                continue

            mappings[read_id]+=len(maps.split(','))

    with open(out_prefix+"/asm_basemasking_overrep_kmers.bed", 'w') as kbed,\
        open(out_prefix+"/asm_basemasking_overrep_kmers.histogram", 'w') as hist:
        hist.write("#{}\t{}\n".format('Placements', 'n_Kmers'))
        khist=defaultdict(int)
        mask_count=0

        for mapping in mappings:
            holder=mapping.split(':')
            scaffold=holder[-2]
            start,stop=holder[-1].split('-')
            start=int(start)
            stop=int(stop)
            placements=mappings[mapping]
            #populate the histogram dictionary
            khist[placements]+=1

            if placements > max_placements:
                mask_count+=stop-start
                assembly[scaffold][start:stop]=mask*(stop-start)
                kbed.write('{}\t{}\t{}\t{}'.format(scaffold, start, stop, placements))

        for count in sorted(khist):
            hist.write("{}\t{}\n".format(count,khist[count]))

    return assembly
if __name__=='__main__':
    logging.basicConfig(level=logging.INFO)

    args=parse_arguments()
    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    assembly=load_assembly(args.assembly)

    logging.info("Masking assembly with provided masks...")
    assembly, mask_count=mask_regions(args.masks, assembly=assembly)

    logging.info("Writing base-masked assembly to disk...")
    write_assembly(assembly, args.out_folder+'/asm_basemasking.fa')
    #do some stats on masked bases here

    logging.info('Generating kmers for kmermasking...')
    kmerfp=get_kmers(assembly, args.ksize, args.step, args.out_folder)

    logging.info('Running gem-index and gem-mapper to get overrepresented kmers...')
    map_file=map_kmers(args.out_folder+'/asm_basemasking.fa', kmerfp, args.threads, \
                'asm_basemasking' ,args.out_folder, args.gem)

    logging.info('Count kmers and mask overrepresented ones...')
    assembly=count_mask_kmer_mappings(map_file, assembly, args.out_folder)

    logging.info('Writing assembly with masked overrepresented kmers to disk...')
    write_assembly(assembly, args.out_folder+'/asm_basemasking.kmermasking.fa')
