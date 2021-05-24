#!/usr/bin/env python
import os
import subprocess

db_path = '/bioseq/biodb/BLAST/Proteins2/nr'
output_path = '/dorotheeh/zivlang/output/blast/blastx/'

# files = ['btcaA1.fa', 'btcaB1.fa']
input_path = '/dorotheeh/zivlang/output/fasta_files/'


# for file in files:
def blastx(file):
    command = f'blastx -query {input_path}{file} -db {db_path} -max_hsps 1 -max_target_seqs 50 -num_threads 10 ' \
              f'-evalue 1e-5 -out {output_path + file[:-2]}txt -outfmt'.split()
    command.append("6 qseqid sseqid pident staxids sskingdoms qstart qend qlen length sstart send slen evalue mismatch "
                   "gapopen bitscore stitle")
    print(command)
    subprocess.run(command, shell=False)
    # subprocess.run(command, shell=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('file_name')
    args = parser.parse_args()
    file_name = args.file_name
    blastx(file_name)
