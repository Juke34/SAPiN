#!/usr/bin/env python3
from SAPiN.version import __version__
import pysam
import os
import sys
import gzip
import json
import pprint
import time
import shutil
import logging
import argparse
import re

shameless_plug="""
    #############################################################################
    # SAPiN v{str1}                                                                #
    # IRD - French National Research Institute for Sustainable Development      #
    # Author: Jacques Dainat                                                    #
    # Please visit https://github.com/Juke34/SAPiN for more information         #
    #############################################################################
\n""".format(str1=__version__)

TODO="""
TODO: ...
"""

# define constants
SCRIPT_DIR=os.path.dirname(__file__)


##########################
#        MAIN            #
##########################

def main():
    """
    SAPiN aims to summarize SAM/BAM read alignment by pileup or reads at each position in a tabulated way. 
    More convenient as a mpileup format and containing extra information.
    """

    # create a parser for the arguments
    parser = argparse.ArgumentParser( description = __doc__ )

    #positional arguments - none - 

    #single character parameter
    parser.add_argument("-i", "--input", default=None, help="String. Path to the SAM/BAM input file.")
    parser.add_argument("-o", "--output", default=None, help="Output filename.")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="Decrease verbosity.")
    parser.add_argument("-v", "--verbose", action="count", default=2, help="Increase verbosity.")
    parser.add_argument("-z", "--gzip", default=False, action="store_true", help="Gzip output file.")
    parser.add_argument("-s", "--shame", action="store_true", help="Suppress the shameless plug.")
    parser.add_argument("-cf", "--cover_filter", default=None, type=int, help="Integer, filter output to report only site with coverage >= <Integer>.")
    parser.add_argument("-qf", "--quality_filter", default=None, type=int, help="Integer, filter output to report only site with quality >= <Integer>.")

    args = parser.parse_args()

    # take care of the loggin depth
    logging.basicConfig(format='%(asctime)s %(levelname)s %(module)s: %(message)s',
                        level = (5-args.verbose+args.quiet)*10,
                        datefmt="%H:%M:%S")

    # store cover_filter value
    cover_filter = args.cover_filter

    # store quality_filter value
    quality_filter = args.quality_filter

    if args.output:
        outfile = args.output
        if args.gzip:
            if not outfile.endswith(".gz"):
                outfile += ".gz"
            outfile = gzip.open(outfile, "wt")
        else:
            outfile = open(outfile, "w")
    else:
        outfile = sys.stdout

    # Print Info
    if not args.shame:
        sys.stderr.write(shameless_plug)

    logging.info("Reading " + args.input + " file")

    # Read sam/bam file
    samfile = pysam.AlignmentFile(args.input, "rb")
    print(f"SEQID\tPOS\tREF\tQUAL\tA\tT\tG\tC\tN\tINS\tDEL\tIUPAC\tCOV\tREADID") 
    dict_info = {
        "SEQID": "",  
        "POS": "",
        "REF": "",
        "QUAL": "",
        "A": "",
        "T": "",
        "G": "",
        "C": "",
        "N": "",
        "IUPAC": "",
        "COV": ""        
        }
    for pileupcolumn in samfile.pileup(add_indels=True):
      
        dict_info["COV"] = pileupcolumn.n
        dict_info["POS"] = pileupcolumn.pos
        dict_info["SEQID"] = pileupcolumn.reference_name

        dict_nuc = {
        "A": 0,
        "T": 0,
        "G": 0,
        "C": 0,
        "N": 0,
        "IUPAC": 0,
        "INS": 0,
        "DEL": 0,
        "QUAL": "/"       
        }

        nb_qual=0
        for pileupread in pileupcolumn.pileups:
  
            # query position is None if is_del or is_refskip is set. (del is noted D is CIGAR and N for skipped region) This is to distinguish between deletions in exons and large skips due to introns.
            if pileupread.is_del or pileupread.is_refskip:
                logging.info (f"deletion of size {pileupread.indel+1}")
                dict_nuc["DEL"]+=1
            elif pileupread.indel > 0:
                logging.info (f"insertion of size {pileupread.indel+1}")
                dict_nuc["INS"]+=1
            else:
                # take care of quality. Might not exists as for minimap sam style output
                if pileupread.alignment.query_qualities:  
                    qual = pileupread.alignment.query_qualities[pileupread.query_position]
                    nb_qual+=1
                    if dict_nuc["QUAL"] == "/":
                         dict_nuc["QUAL"] = qual
                    else:
                        dict_nuc["QUAL"] += qual
                
                nuc = pileupread.alignment.query_sequence[pileupread.query_position]
                if nuc in dict_nuc:
                    dict_nuc[nuc] += 1
                elif nuc.lower() == "u":
                    dict_nuc["T"] += 1
                elif nuc.lower() == "n":
                    dict_nuc["N"] += 1
                elif nuc in ["R","Y","S","W","K","M","B","D","H","V"]:
                     dict_nuc["IUPAC"] += 1
                else:
                    print(f"{nuc} not expected value in read") 

        #Fill dict_info by dict_nuc
        dict_info["A"] = dict_nuc["A"]
        dict_info["T"] = dict_nuc["T"]
        dict_info["G"] = dict_nuc["G"]
        dict_info["C"] = dict_nuc["C"]
        dict_info["N"] = dict_nuc["N"]
        dict_info["INS"] = dict_nuc["INS"]
        dict_info["DEL"] = dict_nuc["DEL"]
        dict_info["IUPAC"] = dict_nuc["IUPAC"]
        if dict_nuc["QUAL"] == "/":
             dict_info["QUAL"] = "/"
        else:
            dict_info["QUAL"] = round(dict_nuc["QUAL"]/nb_qual,2)

         
        print_line=1
        # apply cover filter
        if cover_filter:
            if dict_info["COV"] < cover_filter:
                print_line=None
        # apply quality filter
        if quality_filter:
            if dict_info["QUAL"] == "/":
                print_line=None
            elif dict_info["QUAL"] < quality_filter:
                 print_line=None

        # print result line for studied position
        if print_line:
            print(dict_info["SEQID"] + "\t" + str(dict_info["POS"]) + "\t" + dict_info["REF"] + "\t" + str(dict_info["QUAL"]) + "\t" + str(dict_info["A"]) + "\t"
                + str(dict_info["T"]) + "\t" + str(dict_info["G"]) + "\t" + str(dict_info["C"]) + "\t" + str(dict_info["N"]) + "\t" + str(dict_info["INS"]) + "\t" 
                + str(dict_info["DEL"]) + "\t" + str(dict_info["IUPAC"]) + "\t" + str(dict_info["COV"]) )

    # close file       
    samfile.close()

    sys.stderr.write( """Work done\n""")