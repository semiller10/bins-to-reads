'''
Recover all reads associated with bins in Anvi'o database
'''

import argparse
import os.path
import pandas as pd
import subprocess

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
            '-o', tmp_out
        ])

        tmp_out_1 = os.path.join(
            os.path.dirname(tmp_out), 
            os.path.splitext(os.path.basename(tmp_out))[0] + '.sorted.fasta'
        )
        subprocess.call(
            ['seqkit', 'sort', tmp_out], 
            stdout=tmp_out_1
        )
        subprocess.call(['mv', tmp_out_1, tmp_out])

    return

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