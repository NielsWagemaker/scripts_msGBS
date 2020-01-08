import argparse
import subprocess
import os
parser = argparse.ArgumentParser(description='Process input files')
# input files
parser.add_argument('-r1',
                    help='forward reads fastq')
parser.add_argument('-r2',
                    help='reverse reads fastq')
parser.add_argument('-o', '--output',
                    help='merged reads fastq')
parser.add_argument('--barcodes',
                    help='max barcode length used to trim joined reads')
parser.add_argument('--cycles', default='150',
                    help='Number of sequencing cycles / read length')
args = parser.parse_args()


def run_subprocess(cmd, args, log_message):
    "Run subprocess under standardized settings"
    # force the cmds to be a string.
    if len(cmd) != 1:
        cmd = [" ".join(cmd)]
    with open(".snakemake/test.log", 'w') as log:
        log.write("now starting:\t%s\n" % log_message)
        log.write('running:\t%s\n' % (' '.join(cmd)))
        log.flush()
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, executable='bash')
        stdout, stderr = p.communicate()
        stdout = stdout.decode().replace('\r', '\n')
        stderr = stderr.decode().replace('\r', '\n')
        if stdout:
            log.write('stdout:\n%s\n' % stdout)
        if stderr:
            log.write('stderr:\n%s\n' % stderr)
        log.write('finished:\t%s\n\n' % log_message)
    return 0

"""join fastq files with 'NNNN' between forward and reverse complemented reverse read"""
# get max length of forward and reverse barcodes
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
            if line == '\n':
                continue
            split_line = line.rstrip('\n').split('\t')
            try:
                # TODO: make control nucleotide explicit option in barcode file, now harccoded!
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






# Trim the reads up to the min expected length to improve de novo reference creation for joined reads
cmd = ["paste <(seqtk seq %s | cut -c1-%s) " % (args.r1, max_len_R1) +
       "<(seqtk seq -r %s |cut -c1-%s|seqtk seq -r -)|cut -f1-5" % (args.r2, max_len_R2) +
       "|sed '/^@/!s/\t/NNNNNNNN/g'| sed s/+NNNNNNNN+/+/g| sed 's/ /\t/' | cut -f1,2 |  pigz -p 1 -c > %s" % (args.output)]
log = "Combine joined fastq file into single fasta file"
print(cmd)
run_subprocess(cmd, args, log)