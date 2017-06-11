'''
Split protein sequences onto individual lines
'''

import argparse
import os

from Bio import SeqIO

def main(fin_path, fout_path):

    with open(fin_path, 'r') as fin_handle:
        with open(fout_path, 'w') as fout_handle:
            for record in SeqIO.parse(fin_handle, 'fasta'):
                protein_seqs = record.seq._data.strip('*').split('*')
                split_count = 0
                for s in protein_seqs:
                    if len(s) > 5:
                        fout_handle.write('>' + record.id + '_' + str(split_count) + '\n')
                        fout_handle.write(s + '\n')
                    split_count += 1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__file__.__doc__)
    parser.add_argument('--fin_path', default='test.faa', help='Path to input faa file')
    args = parser.parse_args()
    fin_path = os.path.abspath(args.fin_path)
    fout_path = os.path.join(
        os.path.dirname(fin_path),
        os.path.splitext(os.path.basename(fin_path))[0] + '.split.faa'
        )
    main(fin_path, fout_path)