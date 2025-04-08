#!/usr/bin/env python3
import pysam
import gffutils
import matplotlib.pyplot as plt
import numpy as np
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
import importlib.metadata

script_name = "SAPiN"
script_version = importlib.metadata.version(script_name)

shameless_plug="""
    #############################################################################
    # SAPiN v{str1}                                                               #
    # IRD - French National Research Institute for Sustainable Development      #
    # Author: Jacques Dainat                                                    #
    # Please visit https://github.com/Juke34/SAPiN for more information         #
    #############################################################################
\n""".format(str1=script_version)

TODO="""
TODO: ...
"""
# define constants
SCRIPT_DIR=os.path.dirname(__file__)

##########################
#        Functions       #
##########################



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
    parser.add_argument("-a", "--ali", required=True, default=None, help="String. Path to the sorted BAM input file.")
    parser.add_argument("-f", "--fasta", default=None, help="String. Path to the reference fasta file used to align the reads against.")
    parser.add_argument("-g", "--gff", default=None, help="String. Path to the reference gff.")
    parser.add_argument("-o", "--output", default=None, help="Output filename.")
    parser.add_argument("-p", "--plot", action="store_true", help="Boolean, plot level of mutation per position (sapin_plot.svg by default. If outpout provided output.svg).")
    parser.add_argument("-q", "--quiet", action="count", default=0, help="Decrease verbosity.")
    parser.add_argument("-v", "--verbose", action="count", default=2, help="Increase verbosity.")
    parser.add_argument("-z", "--gzip", default=False, action="store_true", help="Gzip output file.")
    parser.add_argument("-s", "--shame", action="store_true", help="Suppress the shameless plug.")
    parser.add_argument("-mf", "--mutation_filter", default=0, type=int, help="Integer, filter output to report only site where percentage of mutation >= <Integer>.")
    parser.add_argument("-cf", "--cover_filter", default=None, type=int, help="Integer, filter output to report only site with coverage >= <Integer>.")
    parser.add_argument("-bqf", "--base_quality_filter", default=0, type=int, help="Integer, filter output to report only site with base quality >= <Integer>.")
    parser.add_argument("-mqf", "--mapping_quality_filter", default=0, type=int, help="Integer, filter output to report only site with mapping quality >= <Integer>.")

    args = parser.parse_args()

    # take care of the loggin depth
    logging.basicConfig(format='%(asctime)s %(levelname)s %(module)s: %(message)s',
                        level = (5-args.verbose+args.quiet)*10,
                        datefmt="%H:%M:%S")

    # store cover_filter value
    cover_filter = args.cover_filter

    # store quality_filter value
    base_quality_filter = args.base_quality_filter

    # store quality_filter value
    mapping_quality_filter = args.mapping_quality_filter

    # handle output
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
    logging.info("Reading " + args.ali + " file")

    if not args.fasta:
        logging.info("No reference fasta file provided!")

    # --------- HANDLE INPUT --------

    # Handle input BAM/SAM
    # get extension, try again if .gz
    input_format = "bam"
    filename, file_extension = os.path.splitext(args.ali)
    if file_extension == ".gz":
        filename, file_extension = os.path.splitext(filename)
    if file_extension == ".sam":
        sys.exit(".sam extension detected! SAPiN does not accept sam file. Please use sorted bam file.")

    # Handle input fasta file 
    if args.fasta:
        genome = pysam.FastaFile(args.fasta)

    # dictionary in case we want to plot
    plot_mut={}

    # Handle annotation
    if args.gff:
        filename, file_extension = os.path.splitext(args.gff)
        # Put in memory
        gffdb = gffutils.create_db(args.gff, ':memory:')
        # add introns
        try:
            gffdb.update(gffdb.create_introns())
        except Exception as e: 
            logging.error('Failed to create intron features: '+ str(e))
        # add splice sites
        try:
            gffdb.update(gffdb.create_splice_sites())
        except Exception as e:
            logging.error('Failed to create splice_sites features: '+ str(e))
        # add intergenic
        try:
            gffdb.update(gffdb.interfeatures(gffdb.features_of_type("gene")))
        except Exception as e:
            logging.error('Failed to create intergenic features: '+ str(e))

    # Output header
    header = "SEQID\tPOS\tREF\tQUAL\tA\tT\tG\tC\tN\tINS\tDEL\tIUPAC\tCOV\tCOV_ATGC\tMUT_RAT\tAPOBEC\tADAR\tREGION"
    if args.gff:
        header+="\tCODON\tNUC\tDESC"
    header+="\n"
    outfile.write(header)

    # --- Read sam/bam file ---
    samfile = pysam.AlignmentFile(args.ali, "rb")
    samfile2 = pysam.AlignmentFile(args.ali, "rb")
    
    dict_info = {
        "SEQID": "",  
        "POS": None,
        "REF": ".",
        "QUAL": ".",
        "A": "",
        "T": "",
        "G": "",
        "C": "",
        "N": "",
        "IUPAC": "",
        "COV": "",
        "COV_ATGC": "",
        "APOBEC": "",
        "ADAR": ""         
    }


    # min_base_quality (int) – Minimum base quality. Bases below the minimum qual- ity will not be output. The default is 13.
    # min_mapping_quality (int) – only use reads above a minimum mapping quality. The default is 0.
    cpt=0
    for pileupcolumn in samfile.pileup( min_base_quality=base_quality_filter, min_mapping_quality=mapping_quality_filter): #, add_indels=True):
        cpt+=1
        dict_info["COV"] = pileupcolumn.nsegments
        dict_info["POS"] = pileupcolumn.reference_pos + 1 # position is 0-based so need +1
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
        "QUAL": "."       
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
                    if dict_nuc["QUAL"] == ".":
                         dict_nuc["QUAL"] = qual
                    else:
                        dict_nuc["QUAL"] += qual
                
                nuc = pileupread.alignment.query_sequence[pileupread.query_position].upper()

                if nuc in dict_nuc:
                    dict_nuc[nuc] += 1
                elif nuc == "U":
                    dict_nuc["T"] += 1
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
        if dict_nuc["QUAL"] == ".":
             dict_info["QUAL"] = "."
        else:
            dict_info["QUAL"] = round(dict_nuc["QUAL"]/nb_qual,2)

        # Check if we print line based on Coverage and Quality threshold
        print_line=1
        # apply cover filter
        if cover_filter:
            if dict_info["COV"] < cover_filter:
                print_line=None

        # Do time consuming things only if we still conisder the line
        if print_line:

            # get ref base
            ref_base = "."
            if args.fasta:
                ref_base = genome.fetch(dict_info["SEQID"], dict_info["POS"]-1, dict_info["POS"]).upper()
            """ 
            else:
                # If bam format because fetch is not possible on sam 
                if input_format == "bam":
                    logging.warning("experimental!")
                    # fetch , fetch a record in a region, specified either by contig, start, and end (which are 0-based, half-open); or alternatively by a samtools region string (which is 1-based inclusive).
                    for alignment in samfile2.fetch(dict_info["SEQID"], dict_info["POS"]-1, dict_info["POS"]+1):
                        ref_base = find_base_in_alignment(alignment, dict_info["POS"])
                        # we continue until we find a base
                        if ref_base:
                            break    # break here because we found the original base
                else: 
                    logging.warning("Not possible to reconstruct the reference sequence from sam. Please provide the fasta file if you need reference information!")
                """
            dict_info["REF"]=ref_base 
            
            # COV_FILTER
            dict_info["COV_ATGC"] = dict_info["A"] + dict_info["T"] + dict_info["G"] + dict_info["C"]
            
            # apply mutation filter
            cov_mut=dict_info["COV_ATGC"]+dict_info["INS"]+dict_info["DEL"] # Need to get coverage for mutations checked
            nb_mut=0
            nb_not_mut=0
            for nuc in ["A","T","G","C","INS","DEL"]:
                if nuc != ref_base:
                    nb_mut+=dict_info[nuc]
            perc_mut = round(nb_mut*100/cov_mut,2)
            if perc_mut <= args.mutation_filter:
                print_line=None
            else:
                dict_info["MUT_RAT"]=perc_mut

            # check to be sure that mutation filter didnt deactivate the line
            if print_line:           
                
                if args.plot:
                    plot_mut[dict_info["POS"]]=perc_mut

                # APOBEC
                dict_info["APOBEC"]="."
                if dict_info["REF"] == "C" and dict_info["COV_ATGC"] > 0:
                    dict_info["APOBEC"] = round(dict_info["T"] / dict_info["COV_ATGC"] * 100, 1)
                if dict_info["REF"] == "G" and dict_info["COV_ATGC"] > 0:
                    dict_info["APOBEC"] = round(dict_info["A"] / dict_info["COV_ATGC"] * 100, 2)

                # ADAR
                dict_info["ADAR"]="."
                if dict_info["REF"] == "A" and dict_info["COV_ATGC"] > 0:
                    dict_info["ADAR"] = round(dict_info["G"] / dict_info["COV_ATGC"] * 100, 2)
                if dict_info["REF"] == "T" and dict_info["COV_ATGC"] > 0:
                        dict_info["ADAR"] = round(dict_info["C"] / dict_info["COV_ATGC"] * 100, 2)

                # get region
                ref_region="."
                if dict_info["ADAR"]!= "." or dict_info["APOBEC"] != "." and args.fasta: 
                    start = dict_info["POS"]-6
                    if start < 0:
                        start = 0
                    ref_region = genome.fetch(dict_info["SEQID"], start, dict_info["POS"]+5).upper()

                # get annotation
                if args.gff:
                    # get all features that overlap at the position
                    features = gffdb.region(seqid=dict_info["SEQID"], start=dict_info["POS"], end=dict_info["POS"])
                    # gather feature type and attribute information in a single line 
                    # e.g: FeatureType:Attribute-key:Attribute-values;Attribute-key:Attribute-values@@FeatureType:Attribute-key:Attribute-values;Attribute-key:Attribute-values
                    # example gene:ID=gene-1;Name=E6@@mRNA:ID=nbis-rna-1;Parent=gene-1;Name=E6@@exon:ID=nbis-exon-1;Parent=nbis-rna-1;Name=E6@@CDS:ID=cds-1;Parent=nbis-rna-1;Name=E6
                    info_list=[]
                    AA = "."
                    phase = "."
                    for feature in features:
                        featuretype = feature.featuretype

                        # Get AA if CDS
                        if(args.fasta):
                            if featuretype.upper() == "CDS":                 
                                if feature.strand == "+" or feature.strand == "-":
                                    if feature.frame == "0" or feature.frame == "1" or feature.frame == "2":
                                        cds_piece_length = dict_info["POS"] - feature.start + 1
                                        cds_mod = cds_piece_length % 3
                                        codon_shift = 3 - cds_mod
                                        AA = genome.fetch(dict_info["SEQID"], dict_info["POS"] - codon_shift + int(feature.frame), dict_info["POS"] + 3 - codon_shift + int(feature.frame) ).upper()
                                        phase = codon_shift
                                        #print("feature.start" + str(feature.start) +  "cds_mod " + str(cds_mod)  + " cds_piece_length " + str(cds_piece_length) + " AA="+AA )
                        
                        att_key=featuretype+":"
                        info_att_list = []
                        for attribute in feature.attributes.items():  
                            info_att =  attribute[0]+"="
                            info_att +=  ",".join(attribute[1])
                            info_att_list.append(info_att)
                        att_val = ";".join(info_att_list)
                        info_list.append(att_key+att_val)
                    info = "@@".join(info_list)

            if print_line:
                outfile.write(dict_info["SEQID"] + "\t" + str(dict_info["POS"]) + "\t" + dict_info["REF"] + "\t" + str(dict_info["QUAL"]) + "\t" + str(dict_info["A"]) + "\t"
                    + str(dict_info["T"]) + "\t" + str(dict_info["G"]) + "\t" + str(dict_info["C"]) + "\t" + str(dict_info["N"]) + "\t" + str(dict_info["INS"]) + "\t" 
                    + str(dict_info["DEL"]) + "\t" + str(dict_info["IUPAC"]) + "\t" + str(dict_info["COV"]) + "\t" + str(dict_info["COV_ATGC"]) + "\t" + str(dict_info["MUT_RAT"]) + "\t" 
                    + str(dict_info["APOBEC"]) + "\t" +  str(dict_info["ADAR"]) + "\t" + ref_region)
                
                if args.gff:
                    outfile.write("\t" + AA + "\t" + str(phase) + "\t" + info)

                outfile.write("\n")
            
        #if cpt == 680:
         #   sys.exit()
 
    # close files     
    samfile.close()
    samfile2.close()
    outfile.close()

    if args.plot:
        try:
            lists = sorted(plot_mut.items()) # sorted by key, return a list of tuples
            x, y = zip(*lists) # unpack a list of pairs into two tuples
            plt.plot(x, y)
            plt.ylabel("Percentage mutation")
            plt.xlabel("Reference position")
            image_format = 'svg' # e.g .png, .svg, etc. 
            image_name = 'sapin_plot.svg'
            if args.output:
                filename, file_extension = os.path.splitext(args.output)
                image_name=filename+".svg"

            plt.savefig(image_name, format=image_format, dpi=1200)

        except Exception as e:
            logging.error('Failed to make the plot: '+ str(e))
    
    sys.stderr.write( """Work done\n""")