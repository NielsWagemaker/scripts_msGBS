#!/usr/bin/env python
from __future__ import print_function
import pandas as pd
import os
import sys
import tempfile
import subprocess
import argparse
from Bio import SeqIO
from Bio import Restriction

###verwerkt blast against NR file###

###Produces three types of files : 1 : Name lists 2: loci lists 3: reference files
###
### Plant genus names came fomr :
### Fungi database from : http://www.mycobank.org/mb/283905

pd.set_option('display.max_columns', 300)
pd.options.mode.chained_assignment = None # getting rid of som error messages

def parse_args():
    """Pass command line arguments"""
    parser = argparse.ArgumentParser(description='Remove PCR duplicates from bam file')
    #input files
    parser.add_argument('-i', '--input',
                        help='outputblast_kingdoms_nt_original_tax.txt input file')
    parser.add_argument('-r', '--ref',
                        help='ref.fa input file')
    parser.add_argument('-F', '--fungi',
                        help='All fungi genus names txt file')
    parser.add_argument('-dir', '--directory',
                        help='output directory')
    parser.add_argument('-f', '--flower',
                        help='All flowering plants genus names txt file')
    parser.add_argument('-b', '--Bryophytes',
                        help='All Bryophytes genus names txt file')
    parser.add_argument('-G', '--Gymnosperms',
                        help='All Gymnosperms genus names txt file')
    parser.add_argument('-P', '--Pteridophytes',
                        help='All Pteridophytes genus names txt file')
    args = parser.parse_args()

    return args

def make_separate_lists(args):
    ### Here i parse the blast results (done at command line) and make a list of contigs that occur in more than one species ###
    print("started making Eukaryote and other list")
    df = pd.read_csv(args.input, sep= '\t', names = ["contig", "hit", "pident", "E_value", "bp_hit", "kindoms",
                                                                                                         "Fullname", "alignment_length", "start", "stop"],index_col=None)
    print("start computing")
    contig =''
    Bacteria_count = 0
    Eukaryota_count = 0
    other_list_count = 0
    AM_count =0
    fungi_count =0
    Archaea_count =0
    Virus_count =0

    number = int(df['contig'].count())
    print("number of lines to process : %s" % number)
    Eukaryota_list = []
    Flower_plants_output_list = []
    Bacteria_list = []
    AM_list = []
    other_fungi_list = []
    Archaea_list = []
    Virus_list = []
    Name_Eukaryota_list = []
    Name_Bacteria_list = []
    Name_nan_list = []
    Name_AM_list = []
    Name_other_fungi_list = []
    Name_Archaea_list = []
    Name_Virus_list = []


    #The lists below are at Genus level
    arbuscular_mycorrhiza_list = ['Funneliformis', 'Claroideoglomus', 'Rhizophagus', 'Redeckera', 'Acaulospora',
                                  'Scutellospora', 'Gigaspora', 'Septoglomus', 'Paraglomus', 'Rhizophagus', 'Dominikia',
                                  'Glomus', 'Kamienskia', 'Oehlia', 'Sclerocystis', 'Ambispora', 'Archaeospora',
                                  'Horodyskia', 'Geosiphon', 'Diskagma', 'Entrophospora', 'Diversispora', 'Cetraspora',
                                  'Dentiscutata', 'Quatunica', 'Racocetra', 'Pacispora', 'Rhizoglomus', 'Oehlia',
                                  'Halonatospora']

    ### make fungi list from txt (fungi at genus level) list ###
    ALL_fungi_list = []
    with open(args.fungi, 'r') as f:
        lines = f.readlines()
        for line in lines:
            ALL_fungi_list.append((line.rstrip()))

            ### make flowering_plants list from txt (fungi at genus level) list ###
    Flowering_plants_list = []
    with open(args.flower, 'r') as f:
        lines = f.readlines()
        for line in lines:
            Flowering_plants_list.append((line.rstrip()))

    ### make Bryophytes_plants list from txt (fungi at genus level) list ###
    Bryophytes_plants_list = []
    with open(args.Bryophytes, 'r') as f:
        lines = f.readlines()
        for line in lines:
            Bryophytes_plants_list.append((line.rstrip()))

    ### make Gymnosperms_plants list from txt (fungi at genus level) list ###
    Gymnosperms_plants_list = []
    with open(args.Gymnosperms, 'r') as f:
        lines = f.readlines()
        for line in lines:
            Gymnosperms_plants_list.append((line.rstrip()))

    ### make Pterophytes_plants list from txt (fungi at genus level) list ###
    Pteridophytes_plants_list = []
    with open(args.Pteridophytes, 'r') as f:
        lines = f.readlines()
        for line in lines:
            Pteridophytes_plants_list.append((line.rstrip()))



    Eukaryota = 0
    Flowering_plants = 0
    Bacteria = 0
    Archaea = 0
    Viruses = 0
    Fungi =0
    AM_Fungi =0
    Bryophytes =0
    Gymnosperms =0
    Pteridophytes =0
    NA = 0

    for i, s in df.iterrows():
        #print i
        if s[5] == "Bacteria":
            if s[7] > 40 :
                if s[3] < 1e-20:
                    Bacteria_list.append(str(s[0]))
                    Name_Bacteria_list.append(s[6])
                    Bacteria_count +=1
                    Bacteria +=1
                    pass
        elif s[5] == "Eukaryota":
            Genera = s[6].split()
            Genera2 = Genera[0]
            if Genera2 in arbuscular_mycorrhiza_list:
                if s[7] > 40 :
                    if s[3] < 1e-20:
                        Genera = Genera[:2]
                        Genera = " ".join(Genera)
                        AM_list.append(str(s[0]))
                        Name_AM_list.append(Genera)
                        AM_count +=1
                        AM_Fungi += 1
                        Fungi += 1
                        Eukaryota += 1
                        pass
            elif Genera2 in ALL_fungi_list:
                if s[7] > 40 :
                    if s[3] < 1e-20:
                        Genera = Genera[:2]
                        Genera = " ".join(Genera)
                        other_fungi_list.append(str(s[0]))
                        Name_other_fungi_list.append(Genera)
                        fungi_count +=1
                        Fungi += 1
                        Eukaryota += 1
                        pass
            else:
                Eukaryota_list.append(str(s[0]))
                Name_Eukaryota_list.append(s[6])
                Eukaryota_count +=1
                Genera = s[6].split()
                Genera2 = Genera[0]
                if Genera2 in Flowering_plants_list:
                    Flowering_plants += 1
                    Eukaryota += 1
                    Flower_plants_output_list.append(str(s[0]))
                elif Genera2 in Bryophytes_plants_list:
                    Bryophytes += 1
                    Eukaryota += 1
                elif Genera2 in Gymnosperms_plants_list:
                    Gymnosperms += 1
                    Eukaryota += 1
                elif Genera2 in Pteridophytes_plants_list:
                    Pteridophytes += 1
                    Eukaryota += 1
                else:
                    Eukaryota += 1
                    pass
        elif s[5] == "Archaea":
            if s[7] > 40 :
                if s[3] < 1e-20:
                    Archaea_list.append(str(s[0]))
                    Name_Archaea_list.append(s[6])
                    Archaea_count +=1
                    Archaea +=1
                    pass
        elif s[5] == "Viruses":
            if s[7] > 40 :
                if s[3] < 1e-20:
                    Virus_list.append(str(s[0]))
                    Name_Virus_list.append(s[6])
                    Virus_count +=1
                    Viruses +=1
                    pass
        elif s[5] == "N/A":
            if s[7] > 40 :
                if s[3] < 1e-20:
                    Eukaryota_list.append(str(s[0]))
                    Name_nan_list.append(s[6])
                    other_list_count +=1
                    NA +=1
                    pass
        else:
            if s[7] > 40 :
                if s[3] < 1e-20:
                    Eukaryota_list.append(str(s[0]))
                    Name_nan_list.append(s[6])
                    other_list_count +=1
                    NA +=1
                    pass


    Flower_plants_output_list = set(Flower_plants_output_list)
    Bacteria_list = set(Bacteria_list)
    Eukaryota_list = set(Eukaryota_list)
    AM_list = set(AM_list)
    other_fungi_list =set(other_fungi_list)
    Archaea_list =set(Archaea_list)
    Virus_list = set(Virus_list)

    Name_Eukaryota_list = set(Name_Eukaryota_list)
    Name_Bacteria_list = set(Name_Bacteria_list)
    Name_AM_list = set(Name_AM_list)
    Name_other_fungi_list = set(Name_other_fungi_list)
    Name_Archaea_list = set(Name_Archaea_list)
    Name_Virus_list = set(Name_Virus_list)

    print('hit to Eukaryota                           :%s'% Eukaryota)
    print('       of which Flowering plants                    :%s'% Flowering_plants)
    print('       of which Fungi                               :%s'% Fungi)
    print('       of which AM Fungi                            :%s'% AM_Fungi)
    print('       of which Bryophytes                          :%s'% Bryophytes)
    print('       of which Gymnosperms                         :%s'% Gymnosperms)
    print('       of which Pteridophytes                       :%s'% Pteridophytes)
    print('       no Plant Genera available or Other Eukaryote :%s'% (Eukaryota-(Flowering_plants+Fungi+AM_Fungi+Bryophytes+Gymnosperms+Pteridophytes)))
    print('hit to Bacteria                            :%s'% Bacteria)
    print('hit to Archaea                             :%s'% Archaea)
    print('hit to Viruses                             :%s'% Viruses)

    print('hit to NA or nan                           :%s'% NA)

    print("wait")

    bact_out_handle = open(os.path.join(args.directory, 'output_blast/bact_list.txt'), 'w')
    for item in Bacteria_list:
        item = str(item)
        bact_out_handle.write(item + '\n')
    bact_out_handle.close()

    euk_out_handle = open(os.path.join(args.directory, 'output_blast/euk_list.txt'), 'w')
    for item in Eukaryota_list:
        item = str(item)
        euk_out_handle.write(item + '\n')
    euk_out_handle.close()

    AM_fungi_out_handle = open(os.path.join(args.directory, 'output_blast/AM_fungi_list.txt'), 'w')
    for item in AM_list:
        item = str(item)
        AM_fungi_out_handle.write(item + '\n')
    AM_fungi_out_handle.close()

    other_fungi_out_handle = open(os.path.join(args.directory, 'output_blast/other_fungi_list.txt'), 'w')
    for item in other_fungi_list:
        item = str(item)
        other_fungi_out_handle.write(item + '\n')
    other_fungi_out_handle.close()

    Archaea_out_handle = open(os.path.join(args.directory, 'output_blast/Archaea_list.txt'), 'w')
    for item in Archaea_list:
        item = str(item)
        Archaea_out_handle.write(item + '\n')
    Archaea_out_handle.close()

    Virus_out_handle = open(os.path.join(args.directory, 'output_blast/Virus_list.txt'), 'w')
    for item in Virus_list:
        item = str(item)
        Virus_out_handle.write(item + '\n')
    Virus_out_handle.close()

    ### here i parse the ref.fa and sort out the false (overspecies duplicate contigs) and true contigs ###
    print("")
    print("started producing Eukaryota AND Bactria reference seq lists")

    with open(args.ref, 'r') as f:
        lines = f.readlines()
    number = int(len(lines))
    print("number of lines to process : %s" % number)
    teller =0
    flowerdatalist = []
    bactdatalist = []
    eukdatalist = []
    AMdatalist = []
    otherFUNGIlist = []
    Archaealist = []
    Viruslist =[]
    for line in lines:
        teller += 1
        if not teller % 1000000:
            print('%s lines processed of %s)' % (teller, number))
        if line.startswith('>'):
            string = line[1:]
            string = string.rstrip()
            string =str(string)
            if string in Flower_plants_output_list:
                flowerdatalist.append(line)
                eukdatalist.append(line)
                true =1
            elif string in Eukaryota_list:
                eukdatalist.append(line)
                true =2
            elif string in Bacteria_list:
                bactdatalist.append(line)
                true =3
            elif string in AM_list:
                AMdatalist.append(line)
                true =4
            elif string in other_fungi_list:
                otherFUNGIlist.append(line)
                true =5
            elif string in Archaea_list:
                Archaealist.append(line)
                true =6
            elif string in Virus_list:
                Viruslist.append(line)
                true =7
            else:
                eukdatalist.append(line)
                true =2
        else:
            #TODO: make sure all nnnnnnnn are NNNNNNNN ?
            if true == 1:
                flowerdatalist.append(line)
                eukdatalist.append(line)
            elif true == 2:
                eukdatalist.append(line)
            elif true == 3:
                bactdatalist.append(line)
            elif true == 4:
                AMdatalist.append(line)
            elif true == 5:
                otherFUNGIlist.append(line)
            elif true == 6:
                Archaealist.append(line)
            elif true == 7:
                Viruslist.append(line)

    print("")
    print("started writing Eukaryota reference list")
    out_handle = open('%s/output_denovo/refBlasted.fa' % args.directory, 'w')
    teller =0
    number = int(len(eukdatalist))
    print("number of items to process : %s" % number)
    for item in eukdatalist:
        teller +=1
        out_handle.write(item)
        if not teller % 100000:
            print('%s lines processed of %s)' % (teller, number))
    out_handle.close()

    print("")
    print("started writing flowering plants reference list")
    out_handle = open('%s/output_blast/Ref_flowering_plants.txt' % args.directory, 'w')
    teller =0
    number = int(len(flowerdatalist))
    print("number of items to process : %s" % number)
    for item in flowerdatalist:
        teller +=1
        out_handle.write(item)
        if not teller % 100000:
            print('%s lines processed of %s)' % (teller, number))
    out_handle.close()

    print("")
    print("started writing Bacteria reference list")
    out_handle = open('%s/output_blast/Ref_Bactria.txt' % args.directory, 'w')
    teller =0
    number = int(len(bactdatalist))
    print("number of items to process : %s" % number)
    for item in bactdatalist:
        teller +=1
        out_handle.write(item)
        if not teller % 100000:
            print('%s lines processed of %s)' % (teller, number))
    out_handle.close()

    print("")
    print("started writing AM fungi reference list")
    out_handle = open('%s/output_blast/Ref_AM_Fungi.txt' % args.directory, 'w')
    teller =0
    number = int(len(AMdatalist))
    print("number of items to process : %s" % number)
    for item in AMdatalist:
        teller +=1
        out_handle.write(item)
        if not teller % 100000:
            print('%s lines processed of %s)' % (teller, number))
    out_handle.close()

    print("")
    print("started writing Other fungi reference list")
    out_handle = open('%s/output_blast/Ref_other_Fungi.txt' % args.directory, 'w')
    teller =0
    number = int(len(otherFUNGIlist))
    print("number of items to process : %s" % number)
    for item in otherFUNGIlist:
        teller +=1
        out_handle.write(item)
        if not teller % 100000:
            print('%s lines processed of %s)' % (teller, number))
    out_handle.close()

    print("")
    print("started writing Archaea reference list")
    out_handle = open('%s/output_blast/Ref_Archaea.txt' % args.directory, 'w')
    teller =0
    number = int(len(Archaealist))
    print("number of items to process : %s" % number)
    for item in Archaealist:
        teller +=1
        out_handle.write(item)
        if not teller % 100000:
            print('%s lines processed of %s)' % (teller, number))
    out_handle.close()

    print("")
    print("started writing Virus reference list")
    out_handle = open('%s/output_blast/Ref_Virus.txt' % args.directory, 'w')
    teller =0
    number = int(len(Viruslist))
    print("number of items to process : %s" % number)
    for item in Viruslist:
        teller +=1
        out_handle.write(item)
        if not teller % 100000:
            print('%s lines processed of %s)' % (teller, number))
    out_handle.close()


    ### Procuding NAME txt Documents

    print("")
    print("started writing Eukaryota  Name list")
    out_handle = open('%s/output_blast/EUKARYOTA_NAMES.txt' % args.directory, 'w')
    teller =0
    number = int(len(Name_Eukaryota_list))
    print("number of items to process : %s" % number)
    for item in Name_Eukaryota_list:
        teller +=1
        out_handle.write(item+ '\n')
    out_handle.close()

    print("started writing Bacteria Name list")
    out_handle = open('%s/output_blast/BACTERIA_NAMES.txt' % args.directory, 'w')
    teller =0
    number = int(len(Name_Bacteria_list))
    print("number of items to process : %s" % number)
    for item in Name_Bacteria_list:
        teller +=1
        out_handle.write(item+ '\n')
    out_handle.close()

    print("started writing AM fungi Name list")
    out_handle = open('%s/output_blast/AM_FUNGI_NAMES.txt' % args.directory, 'w')
    teller =0
    number = int(len(Name_AM_list))
    print("number of items to process : %s" % number)
    for item in Name_AM_list:
        teller +=1
        out_handle.write(item+ '\n')
    out_handle.close()

    print("started writing other fungi Name list")
    out_handle = open('%s/output_blast/OTHER_FUNGI_NAMES.txt' % args.directory, 'w')
    teller =0
    number = int(len(Name_other_fungi_list))
    print("number of items to process : %s" % number)
    for item in Name_other_fungi_list:
        teller +=1
        out_handle.write(item+ '\n')
    out_handle.close()

    print("started writing Archaea Name list")
    out_handle = open('%s/output_blast/Archaea_NAMES.txt' % args.directory, 'w')
    teller =0
    number = int(len(Name_Archaea_list))
    print("number of items to process : %s" % number)
    for item in Name_Archaea_list:
        teller +=1
        out_handle.write(item+ '\n')
    out_handle.close()

    print("started writing Virus Name list")
    out_handle = open('%s/output_blast/Virus_NAMES.txt' % args.directory, 'w')
    teller =0
    number = int(len(Name_Virus_list))
    print("number of items to process : %s" % number)
    for item in Name_Virus_list:
        teller +=1
        out_handle.write(item+ '\n')
    out_handle.close()


def main():
    """main function loop"""
    args = parse_args()
    #args = remove_hashtaq(args)
    out_lists = make_separate_lists(args)


if __name__ == '__main__':
    main()
