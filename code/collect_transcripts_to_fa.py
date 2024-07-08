import sys
import os
import gzip
from Bio import SeqIO

if __name__ == "__main__":

    ids_file = "id_list.txt"
    ids = {}
    with open(ids_file) as f:
        for g in f:
            organism, gid = g.rstrip().split("\t")
            ids[gid] = organism

    top_transcript_dir = "../transcripts"

    for path, subdirs, files, in os.walk(top_transcript_dir):
        for f in files:
            if f[-6:] == '.fa.gz':
                fa_path = f"{path}/{f}"
                with gzip.open(fa_path, 'rt') as fa_file:
                    
                    fasta_sequences = SeqIO.parse(fa_file,'fasta')
                    for s in fasta_sequences:
                        print(s.id)
                        #for gid in ids:
                        #    if s.id == gid:
                        #        print(s)
