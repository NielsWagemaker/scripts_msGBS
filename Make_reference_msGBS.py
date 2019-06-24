#!/usr/bin/env pypy
__author__ = 'thomasvangurp'
# Date created: 20/11/2014 (europe date)
# Function: Pipeline for creation of reference
#Python version: 2.7.3
#External dependencies: vsearch,samtools,vcfutils.pl,rename_fast.py
#Known bugs: None
#Modifications: None
import argparse
import subprocess
import tempfile
import os
import operator
import gzip

origWD = os.getcwd()
os.chdir(origWD)


dependencies = {}
dependencies['samtools'] = '<0.1.18'         #ok
dependencies['vcfutils.pl'] = ''             #ok
dependencies['usearch'] = 'usearch_8.0.1409' #staat in PATH
dependencies['seqtk'] = '1.0-r31'            #staat in path
#dependencies['pear'] = 'v0.9.7'              #zou moeten werken staat in /usr/local/bin
#dependencies['NGmerge'] = vUnknown           #werkt op server
dependencies['pigz'] = ''                    #staat in path

usearch = "usearch"#_8.0.1409_i86osx32"
vsearch = "vsearch"
seqtk = "seqtk"
pear = "pear"
create_consensus = "create_consensus.py"
vcfutils = "vcfutils.pl"

def parse_args():
    "Pass command line arguments"
    parser = argparse.ArgumentParser(description='Process input files')
    #input files
    parser.add_argument('-s','--sequences',
                        help='number of sequences to take for testing, useful for debugging')
    # parser.add_argument('--forward_in',
    #                     help='forward input reads fastq')
    parser.add_argument('--forward',
                        help='forward reads fastq')
    # parser.add_argument('--reverse_in',
    #                     help='reverse input reads fastq')
    #    parser.add_argument('--reverse_tussen',
    #                        help='reverse temp tussen reads fastq')
    parser.add_argument('--reverse',
                        help='reverse reads fastq')
    parser.add_argument('-o','--output',
                        help='merged reads fastq')
    parser.add_argument('--barcodes',
                        help='max barcode length used to trim joined reads')
    parser.add_argument('--cycles',default='150',
                        help='Number of sequencing cycles / read length')
    parser.add_argument('--min_unique_size',default="2",
                        help='Minimum unique cluster size') # TODO: min unique cluster sie was 2 ???? i changed this
    parser.add_argument('-q','--qual_profile',
                        help='qual_profile txt input file')
    parser.add_argument('-t','--tmpdir',
                        help='tmp directory',default='/tmp')
    parser.add_argument('-n','--NGthreads',
                        help='Number of threads to used where multithreading is possible')
    parser.add_argument('-u','--unassembled_output',
                        help='Unassembled_reads - merged.unassembled')
    parser.add_argument('--threads',
                        help='Number of threads to used where multithreading is possible')
    parser.add_argument('--outputdir',
                        help='Optional: output directory')
    parser.add_argument('--log',
                        help='log of output operation')
    parser.add_argument('--samout',#TODO: describe what is the purpose of this file
                        help='sam output of clustering process')
    parser.add_argument('--consensus',#TODO: describe what is the purpose of this file
                        help='consensus output')
    parser.add_argument('--consensus_cluster',#TODO: describe what is the purpose of this file
                        help='consensus clustering output')
    args = parser.parse_args()
    if args.outputdir:
        if not os.path.exists(args.outputdir):
            try:
                os.mkdir(args.outputdir)
            except OSError:
                raise
        args.log = os.path.join(args.outputdir,'make_reference.log')
        args.samout = os.path.join(args.outputdir,'clustering.bam')
        args.consensus = os.path.join(args.outputdir,'consensus.fa')
        args.consensus_cluster = os.path.join(args.outputdir,'consensus_cluster.fa')
    return args

def run_subprocess(cmd,args,log_message):
    "Run subprocess under standardized settings"
    #force the cmds to be a string.
    if len(cmd) != 1:
        cmd = [" ".join(cmd)]
    with open(args.log,'a') as log:
        log.write("now starting:\t%s\n"%log_message)
        log.write('running:\t%s\n'%(' '.join(cmd)))
        log.flush()
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True,executable='bash')
        stdout, stderr = p.communicate()
        stdout = stdout.replace('\r','\n')
        stderr = stderr.replace('\r','\n')
        if stdout:
            log.write('stdout:\n%s\n'%stdout)
        if stderr:
            log.write('stderr:\n%s\n'%stderr)
        log.write('finished:\t%s\n\n'%log_message)
    return 0

def merge_reads(args):
    "Merged reads using NGmerge"
    out_files = {}

    cmd = ['NGmerge']
    #set reads
    cmd+=['-1',args.forward]
    cmd+=['-2',args.reverse]
    #set output directory
    cmd+=['-o',args.output]
    #set number of threads for NGmerge
    cmd+=['-n',args.NGthreads]
    #Set unassemble_output_prefix
    cmd+=['-f',args.unassembled_output]
    #Set location of qual_profile_matrix
    cmd+=['-w',args.qual_profile]
    #set FASTQ quality offset
    cmd+=['-q','33'] # was o ; 140 was niet goed?
    # Maximum input quality score (0-based; def. 40)
    cmd+=['-u','41']
    # Option to gzip (-z) or not (-y) FASTQ output(s)
    cmd+=['-z']
    cmd+=['-j','%sMergedLog.txt'%args.outputdir]
    log = "run NGmerge for merging reads"

    log = "way to go!"
    run_subprocess(cmd,args,log)

    #TODO: Delete input files and output file name that are no longer needed??
    # append output files as dictionary
    # pear put reads in unassembled in REVCOMP (to change this add option -k
    between_files = {'merged':'%s/merged'%args.outputdir + ".fastq.txt.gz",
                     'single_R1':'%s/merged'%args.outputdir + ".unassembled_1.fastq.gz",
                     'single_R2':'%s/merged'%args.outputdir + ".unassembled_2.fastq.gz"}

    out_files = {'merged':'%s/merged'%args.outputdir + ".assembled_correct_header.fastq.gz",
                 'single_R1':'%s/merged'%args.outputdir + ".unassembled_1_correct_header.fastq.gz",
                 'single_R2':'%s/merged'%args.outputdir + ".unassembled_2_correct_header.fastq.gz"}

    cmd = ['zcat %s|sed "s/#BC:Z:/\tBC:Z:/g;s/#RG:Z:/\tRG:Z:/g;s/#ST:Z:/\tST:Z:/g;s/#RN:Z:/\tRN:Z:/g" |pigz -p %s >%s' % (between_files['merged'], args.NGthreads, out_files['merged'])]
    log = "R1 change # to tab :!"
    run_subprocess(cmd,args,log)

    cmd = ['zcat %s|sed "s/#BC:Z:/\tBC:Z:/g;s/#RG:Z:/\tRG:Z:/g;s/#ST:Z:/\tST:Z:/g;s/#RN:Z:/\tRN:Z:/g" |pigz -p %s >%s' % (between_files['single_R1'], args.NGthreads, out_files['single_R1'])]
    log = "R1 change # to tab :!"
    run_subprocess(cmd,args,log)

    cmd = ['zcat %s|sed "s/#BC:Z:/\tBC:Z:/g;s/#RG:Z:/\tRG:Z:/g;s/#ST:Z:/\tST:Z:/g;s/#RN:Z:/\tRN:Z:/g" |pigz -p %s >%s' % (between_files['single_R2'], args.NGthreads, out_files['single_R2'])]
    log = "R1 change # to tab :!"
    run_subprocess(cmd,args,log)

    return out_files


def join_fastq(r1,r2,outfile,args):
    """join fastq files with 'NNNN' between forward and reverse complemented reverse read"""
    #get max length of forward and reverse barcodes
    if args.barcodes:
        with open(args.barcodes) as bc_handle:
            header = bc_handle.readline()[:-1].split('\t')
            barcode_1_index = header.index('Barcode_R1')
            barcode_2_index = header.index('Barcode_R2')
            try:
                wobble_R1_index = header.index('Wobble_R1')
                wobble_R2_index = header.index('Wobble_R2')
            except ValueError:
                wobble_R1_index = None
                wobble_R2_index = None
            barcode_1_max_len = 0
            barcode_2_max_len = 0
            for line in bc_handle:
                split_line = line.rstrip('\n').split('\t')
                try:
                    #TODO: make control nucleotide explicit option in barcode file, now harccoded!
                    wobble_R1_len = int(split_line[wobble_R1_index]) + 1
                    wobble_R2_len = int(split_line[wobble_R2_index]) + 1
                except TypeError:
                    wobble_R1_len = 0
                    wobble_R2_len = 0
                if len(split_line[barcode_1_index]) > barcode_1_max_len:
                    barcode_1_max_len = len(split_line[barcode_1_index])
                if len(split_line[barcode_2_index]) > barcode_2_max_len:
                    barcode_2_max_len = len(split_line[barcode_2_index])
        max_len_R1 = int(args.cycles) - barcode_1_max_len - wobble_R1_len
        max_len_R2 = int(args.cycles) - barcode_2_max_len - wobble_R2_len
    else:
        #no trimming required
        max_len_R1 = 200
        max_len_R2 = 200
    #Trim the reads up to the min expected length to improve de novo reference creation for joined reads
    #Ik heb de -r uit de tweede regel verwijderd: was : "<(seqtk seq -r -A %s |cut -c1-%s)|cut -f1-5" % (r2, max_len_R2)+
    #En weer terug geplaatst; blijkbaar is het toch nodig om deze onlogische orientatie te gebruiken
    #WEL weer aangepast : nu eerst -rA dan cut later -r
    #TODO: Checken of dit nu goed is ??? staat het tweede haakje goed (na seqtk seq -r ))
    cmd = ["paste <(seqtk seq -A %s | cut -c1-%s) " % (r1, max_len_R1) +
           "<(seqtk seq -rA %s |cut -c1-%s|seqtk seq -r )|cut -f1-5" % (r2, max_len_R2)+
           "|sed '/^>/!s/\t/NNNNNNNN/g' |pigz -p %s -c > %s" % (args.threads, outfile)]
    log = "Combine joined fastq file into single fasta file"
    if not os.path.exists(outfile):
        run_subprocess(cmd,args,log)
    return True

def join_non_overlapping(in_files,mapping_dict, args): #ok
    """join non overlapping PE reads"""
    for key in mapping_dict:
        key = key.replace(' ','_')
        if 'rootmono' not in key and key != 'all':
            continue
        left = os.path.join(args.outputdir, '%s.R1.fq.gz' % key)
        right = os.path.join(args.outputdir, '%s.R2.fq.gz' % key)
        joined = os.path.join(args.outputdir, '%s.joined.fq.gz' % key)
        if not os.path.exists(joined):
            join_fastq(left, right, joined, args)
    return in_files

def get_mapping_dict(args):
    """get mapping dict for mono's"""
    mapping_dict = {'all': None}
    with open(args.barcodes) as barcode_handle:
        header = barcode_handle.readline()[:-1].split('\t') # \t = tab delimited
        for line in barcode_handle:
            split_line = line[:-1].split('\t')
            try:
                id = split_line[header.index('Sample')]
            except KeyError:
                raise KeyError('barcode file does not contain column : Sample, please revise your file %s' % args.barcodes)
            mapping_dict[id] = {}
            for k, v in zip(header, split_line):
                if k != 'Sample':
                    mapping_dict[id][k] = v
    return mapping_dict

def trim_split_and_zip(in_files, mapping_dict, args):
    """Trim , split and zip fastq files for mono species"""
    in_files['trimmed'] = {}

    log = 'Zip and split forward reads'
    file_in = in_files['single_R1']
    cmd = seqtk + ' seq %s |tee >(pigz -c > %s)' % (file_in, os.path.join(args.outputdir,'all.R1.fq.gz'))
    #TODO : 'mono' veranderd in "rootmono' : Dit werkte goed
    mono_list = [k for k in mapping_dict.keys() if 'rootmono' in k.lower()]
    #TODO : Hier mogelijk set monolist invoeren om dubbele eruit te halen wanneer er duplo's voor de mono's inzitten
    for k in mono_list[:-1]:
        cmd += "|tee >(grep $'%s\t' -A3|sed '/^--$/d' |pigz -c > %s)" % (k, os.path.join(args.outputdir,'%s.R1.fq.gz'%(k.replace(' ','_'))))
    #add the latest key to mono 1.
    k = mono_list[-1]
    cmd += "|grep $'%s\t' -A3|sed '/^--$/d' |pigz -c > %s" % (k, os.path.join(args.outputdir,'%s.R1.fq.gz'%k.replace(' ','_')))
    if not os.path.exists(os.path.join(args.outputdir,'all.R1.fq.gz')):
        run_subprocess([cmd],args,log)

    log = 'Zip and split reverse reads'
    file_in = in_files['single_R2']
    cmd = seqtk + ' seq -r %s |tee >(pigz -c > %s)' % (file_in, os.path.join(args.outputdir, 'all.R2.fq.gz'))
    mono_list = [k for k in mapping_dict.keys() if 'rootmono' in k.lower()]
    # TODO : Hier mogelijk set monolist invoeren om dubbele eruit te halen wanneer er duplo's voor de mono's inzitten
    for k in mono_list[:-1]:
        cmd += "|tee >(grep $'%s\t' -A3|sed '/^--$/d' |pigz -c > %s)" % (
            k, os.path.join(args.outputdir, '%s.R2.fq.gz' % (k.replace(' ', '_'))))
    # add the latest key to mono 1.
    k = mono_list[-1]
    cmd += "|grep $'%s\t' -A3|sed '/^--$/d' |pigz -c > %s" % (
        k, os.path.join(args.outputdir, '%s.R2.fq.gz' % k.replace(' ', '_')))
    if not os.path.exists(os.path.join(args.outputdir, 'all.R2.fq.gz')):
        run_subprocess([cmd], args, log)

    #Process merged files
    log = 'Zip and split merged reads'
    file_in = in_files['merged']
    cmd = seqtk + ' seq %s |tee >(pigz -c > %s)' % (file_in, os.path.join(args.outputdir, 'all.merged.fq.gz'))
    mono_list = [k for k in mapping_dict.keys() if 'rootmono' in k.lower()]
    # TODO : Hier mogelijk set monolist invoeren om dubbele eruit te halen wanneer er duplo's voor de mono's inzitten
    for k in mono_list[:-1]:
        cmd += "|tee >(grep $'%s\t' -A3|sed '/^--$/d' |pigz -c > %s)" % (
            k, os.path.join(args.outputdir, '%s.merged.fq.gz' % (k.replace(' ', '_'))))
    # add the latest key to mono 1.
    k = mono_list[-1]
    cmd += "|grep $'%s\t' -A3|sed '/^--$/d' |pigz -c > %s" % (
        k, os.path.join(args.outputdir, '%s.merged.fq.gz' % k.replace(' ', '_')))
    if not os.path.exists(os.path.join(args.outputdir, 'all.merged.fq.gz')):
        run_subprocess([cmd], args, log)

    return in_files



def dereplicate_reads(in_files, mapping_dict, args):
    """dereplicate reads using vsearch"""
    in_files['derep'] = {}
    for type in ['joined','merged']:
        for name in [k for k in mapping_dict.keys() if 'rootmono' in k]:
            name = name.replace(' ','_')
            #Voor de zekerheid sorteer ik de file hier nog even op lengte
            file_in_temp = os.path.join(args.outputdir, 'temp%s.%s.fq.gz' % (name, type))
            file_in = os.path.join(args.outputdir, '%s.%s.fq.gz' % (name, type))
            file_out = '.'.join(file_in.split('.')[:-2]) + '.derep.fa'
            cmd = ['vsearch -sortbylength %s --output %s' % (file_in_temp, file_in)]
            log = 'Sort input %s sequences by length.' % name
            if not os.path.exists(file_in):
                run_subprocess(cmd, args, log)
            cmd = [vsearch +' -derep_fulllength %s -sizeout -minuniquesize %s -output %s'%(file_in, args.min_unique_size, file_out)]
            #HIER ZOU -minuniquesize VERHOGEN KUNNEN LEIDEN TOT MINDER CONTIGS 3 gedaan
            log = "Dereplicate full_length of %s using vsearch"%(name)
            if not os.path.exists(file_out):
                run_subprocess(cmd, args, log)
            in_files['derep']['%s_%s_derep' % (name, type)] = file_out
    return in_files


def cluster_consensus(in_files, mapping_dict, args):
    "Cluster concensus with preset id"
    #concatenate and order derep output
    for name in [k for k in mapping_dict.keys() if 'rootmono' in k]:
        #concatenate and order
        name = name.replace(' ','_')
        joined = os.path.join(args.outputdir,'%s.joined.derep.fa' % name)
        merged = os.path.join(args.outputdir,'%s.merged.derep.fa' % name)
        combined = os.path.join(args.outputdir,'%s.combined.derep.fa' % name)
        combined_ordered = os.path.join(args.outputdir,'%s.combined.ordered.derep.fa' % name)
        output_derep = os.path.join(args.outputdir,'%s.clustered.fa' % name)
        output_derep_renamed = os.path.join(args.outputdir,'%s.clustered.renamed.fa' % name)
        cmd = ['cat %s %s > %s' % (joined, merged, combined)]
        log = "combine joined and merged dereplicates reads of %s for sorting" % name
        if not os.path.exists(combined):
            run_subprocess(cmd, args, log)


        cmd = ['vsearch -sortbylength %s --output %s' % (combined, combined_ordered)]
        log = 'Sort derep sequences by length.'
        if not os.path.exists(combined_ordered):
            run_subprocess(cmd, args, log)

        # TODO: NIELS: check if joined and merged clusters from same loci (in the borderzone length of ~280bp )are merged properly
        cmd = ["vsearch -cluster_smallmem %s -id 0.95 -centroids %s -sizeout -strand both"% (combined_ordered,output_derep)]
        # TODO: Variate clustering (was 0.80)
        log = "Clustering consensus with 95% identity"
        if not os.path.exists(output_derep):
            run_subprocess(cmd,args,log)
        # in_files['consensus']['consensus_clustered'] = args.consensus_cluster

        cmd = ['cat %s |rename_fast.py -g %s -n > %s'% (output_derep, name, output_derep_renamed)]
        log = "rename resulting clusters with number and origin name"
        if not os.path.exists(output_derep_renamed):
            run_subprocess(cmd,args,log)
            cmd = ["samtools faidx %s"%output_derep_renamed]
            log = "faidx index %s" % output_derep_renamed
            run_subprocess(cmd, args, log)
    #concatenate all renamed clusters
    ref = os.path.join(args.outputdir, 'ref.fa')
    cmd = ['cat %s/*.renamed.fa > %s' % (args.outputdir, ref)]
    log = 'concatenate all renamed clusters'
    run_subprocess(cmd, args, log)
    #index this new reference sequence
    cmd = ["samtools faidx %s" % ref]
    log = "faidx index %s" % ref
    # TODO : dit moet ik na blast nogmaals doen met de hand...?
    run_subprocess(cmd, args, log)
    return in_files

def check_dependencies():
    """check for presence of dependencies and if not present say where they can be installed"""

    return 0

def clear_tmp(file_dict):
    """clear tmp files"""
    purge_list = []
    for v in file_dict.keys():
        for key,value in file_dict[v].items():
            try:
                if value.startswith('/tmp'):
                    purge_list.append(value)
            except AttributeError:
                if type(value) == type([]):
                    purge_list.append(value[0])
    for item in purge_list:
        print "removing %s" % item
        os.remove(item)
    return 0


def main():
    "Main function loop"
    #Check if os is windows, if so, throw error and tell the user that the software cannot run on windows.
    check_dependencies()
    args = parse_args()
    #Make sure log is empty at start
    if os.path.isfile(args.log):
        os.remove(args.log)
    #Step 1: get mapping dict for mono's
    mapping_dict = get_mapping_dict(args)
    #Step 2 use seqtk to trim merged and joined reads from enzyme recognition site
    files = merge_reads(args)
    #Step 3: Dereplicate all watson and crick reads
    files = trim_split_and_zip(files, mapping_dict, args)
    #Step 4: join the non overlapping PE reads from watson and crick using vsearch
    files = join_non_overlapping(files, mapping_dict, args)
    files = dereplicate_reads(files, mapping_dict, args)
    #step 5: Cluster consensus
    files = cluster_consensus(files, mapping_dict, args)

if __name__ == '__main__':
    main()
