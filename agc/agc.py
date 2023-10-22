#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "OUADAH Lilia"
__copyright__ = "Universite Paris Cité"
__credits__ = ["OUADAH Lilia"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "OUADAH Lilia"
__email__ = "lilia.ouadah@etu.u-paris.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    
    seq = ""
    with gzip.open(amplicon_file, "rt") as fasta :
        for line in fasta :
            if line[0] == ">" :
                if len(seq)>=minseqlen:
                        yield seq
                seq =""
                continue
            seq+=line.strip()
        if len(seq)>=minseqlen:
            yield seq 

def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    seq_counts = {}
    sequences = []

    for seq in read_fasta(amplicon_file, minseqlen):
        if seq in seq_counts:
            seq_counts[seq] += 1
        else:
            seq_counts[seq] = 1

    for sequence, count in seq_counts.items():
        if count >= mincount:
            sequences.append([sequence, count])

    sequences.sort(key=lambda x: x[1], reverse=True)

    return iter(sequences)
    

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    if len(alignment_list) != 2:
        raise ValueError("Alignment list must contain exactly two sequences.")

    sequence1, sequence2 = alignment_list

    if len(sequence1) != len(sequence2):
        raise ValueError("Sequences must have the same length for identity comparison.")

    identical_nucleotides = sum(c1 == c2 for c1, c2 in zip(sequence1, sequence2))
    alignment_length = len(sequence1)

    # Calculate the identity rate as a percentage
    identity_rate = (identical_nucleotides / alignment_length) * 100.0

    return identity_rate

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """

    sequences = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    sequences.sort(key=lambda x: x[1], reverse=True)
    return sequences


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open (output_file, "w") as file:
        for i in range(len(OTU_list)):
            file.write(f">OTU_{i+1} occurrence:{OTU_list[i][1]}\n")
            file.write(textwrap.fill(OTU_list[i][0]) + "\n")


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, 0, 0)
    write_OTU(OTU_list, args.output_file)



if __name__ == '__main__':
    main()
