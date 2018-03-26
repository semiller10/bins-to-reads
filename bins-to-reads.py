'''
Recover all reads associated with bins in Anvi'o database
'''

import argparse
import os.path
import pandas as pd
import subprocess
import sys

def main():

    args = get_args()

    bin_table_hdr = ['bin', 'collection', 'profile_db', 'contig_db', 'bam']
    table = pd.read_csv(args.table, sep='\t', header=None, names=bin_table_hdr)

    out_dir = os.path.dirname(args.out)
    out_name = os.path.splitext(os.path.basename(args.out))[0]
    for i, row in table.iterrows():
        tmp_out = os.path.join(
            out_dir, out_name + '.bin' + str(i) + '.' + row['bin'] + '.tmp.fasta'
        )
        subprocess.call([
            'anvi-get-short-reads-from-bam', 
            '-p', row['profile_db'], 
            '-c', row['contig_db'], 
            '-C', row['collection'], 
            '-b', row['bin'], 
            '-o', tmp_out, 
            row['bam']
        ])

        tmp_out_1 = os.path.join(
            os.path.dirname(tmp_out), 
            os.path.splitext(os.path.basename(tmp_out))[0] + '.sorted.fasta'
        )
        with open(tmp_out_1, 'w') as handle:
            subprocess.call(
                ['seqkit', 'sort', '-w', '0', tmp_out], 
                stdout=handle
            )
        subprocess.call(['mv', tmp_out_1, tmp_out])

        contig_db_name = os.path.splitext(os.path.basename(row['contig_db']))[0]
        hdr_suffix = '_' + contig_db_name + '\n'
        # WHEN ONLY ONE PAIRED-END READ MAPS TO THE CONTIG, 
        # RECOVER THE OTHER READ
        # CURRENTLY IN DEVELOPMENT: REQUIRES READS FILE
        # with open(tmp_out) as handle:
        #     read_hdrs = [hdr.replace(hdr_suffix, '') for hdr in handle.readlines()[::2]]        
        # missing_read_hdrs = determine_missing_reads(read_hdrs)

        # PAIRED-END READS CAN ALSO BE DISCARDED 
        # WHEN ONE READ DOES NOT MAP TO A CONTIG
        # TO ENSURE OUTPUT ENTIRELY OF P-E READS
        with open(tmp_out) as handle:
            mapped_reads = []
            for line in handle.readlines():
                if line[0] == '>':
                    mapped_reads.append(line.replace(hdr_suffix, ''))
                else:
                    mapped_reads.append(line.rstrip())
        pe_reads = remove_unpaired_reads(mapped_reads)
        with open(args.out, 'a') as handle:
            for line in pe_reads:
                if line[0] == '>':
                    handle.write(line + hdr_suffix)
                else:
                    handle.write(line + '\n')

    subprocess.call(['rm', tmp_out])

    return

def remove_unpaired_reads(reads):
    '''
    Remove unpaired reads from a sorted fasta file
    '''

    pe_reads = []
    paired = False
    prev_hdr = reads[0]
    prev_basename = prev_hdr[:-2]
    prev_seq = reads[1]
    for line in reads[2:]:
        if line[0] == '>':
            if paired:
                prev_hdr = line
                prev_basename = prev_hdr[:-2]
                paired = False
            else:
                if line[:-2] == prev_basename:
                    paired = True
                    curr_hdr = line
                else:
                    prev_hdr = line
                    prev_basename = prev_hdr[:-2]
        else:
            if paired:
                pe_reads.append(prev_hdr)
                pe_reads.append(prev_seq)
                pe_reads.append(curr_hdr)
                pe_reads.append(line)
                paired = False
            else:
                prev_seq = line

    return pe_reads

def determine_missing_reads(read_hdrs):
    '''
    Find any missing members of paired reads
    '''
    
    missing_read_hdrs = []
    unpaired = True
    previous_hdr = read_hdrs[0]
    previous_basename = previous_hdr[:-2]
    previous_member = previous_hdr[-1]
    for hdr in read_hdrs[1:]:
        if unpaired:
            if hdr[:-2] == previous_basename:
                unpaired = False
            else:
                if previous_member == '1':
                    missing_read_hdrs.append(previous_basename + '.2\n')
                else:
                    missing_read_hdrs.append(previous_basename + '.1\n')
        else:
            unpaired = True
            previous_hdr = read_hdrs[0]
            previous_basename = hdr[:-2]
            previous_member = hdr[-1]

    return missing_read_hdrs

def get_args():
    '''
    Get command line arguments
    '''

    parser = argparse.ArgumentParser(
        description=(
            'Recover all paired reads associated with bin(s) in an Anvi\'o database.\n'
            'When multiple bins are specified, '
            'reads mapping to the different bins will be output together.\n'
            'Make sure that the Anvi\'o environment is loaded.\n'
            'The following programs must be in the path:\n'
            'anvi-get-short-reads-from-bam\n'
            'seqkit\n'
            'samtools'
        )
    )

    parser.add_argument(
        '-t', 
        '--table', 
        help=(
            'Path to bin table (tsv) with following column format (no headers):\n'
            '1. Bin name from Anvi\'o collection\n'
            '2. Anvi\'o collection name\n'
            '3. Path to Anvi\'o profile db\n'
            '4. Path to Anvi\'o contig db\n'
            '5. Path to bam file of mapped paired end reads'
        )
    )
    parser.add_argument('-o', '--out', help='Path to fasta file output')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()