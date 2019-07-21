import skbio
import os
from Bio import SeqIO
import pandas as pd
from collections import defaultdict
from plumbum import local

class UnifracPipe:
    def __init__(self, discard_unamed=True):
        self.workdir = os.path.abspath(os.path.dirname(__file__))
        self.outdir = os.path.join(self.workdir, "output", "unifrac_testing")
        os.makedirs(self.outdir, exist_ok=True)
        self.all_seqs_input_path = os.path.join(self.workdir, "all_seqs.fasta")
        self.all_seqs_aligned_output_path = os.path.join(self.outdir, 'all_seqs.aligned.fasta')
        self.discard_unnamed = discard_unamed
        self.mafft = local['mafft']
        self.aligned_individual_clade_fasta_paths = []
        self.clade_collection_of_seq_objects_dd = self._populate_clade_collection_of_seq_objects_dd()
        # self.sorted_clade_abund_tup = self._populate_sorted_clade_abund_tup()
        self._align_individual_clade_seq_collections()
        self._create_all_seqs_alignment()

    def _populate_clade_collection_of_seq_objects_dd(self):
        print('seperating sequences by clade for alignment')
        clade_collection_of_seq_objects_dd = defaultdict(list)
        for seq in SeqIO.parse(self.all_seqs_input_path, "fasta"):
            if "_" in seq.id:
                if self.discard_unnamed:
                    continue
                else:
                    clade = seq.id.split("_")[1]
            else:
                clade = seq.id[0]

            clade_collection_of_seq_objects_dd[clade].append(seq)
        return clade_collection_of_seq_objects_dd

    # def _populate_sorted_clade_abund_tup(self):
    #     print('counting the number of sequences for each clade')
    #     clade_sizes = {}
    #     for clade, seq_list in self.clade_collection_of_seq_objects_dd.items():
    #         clade_sizes[clade] = len(seq_list)
    #         print(f'{clade} {len(seq_list)}')
    #
    #     return sorted(clade_sizes.items(), key=lambda x:x[1], reverse=True)

    def _align_individual_clade_seq_collections(self):
        """For each clade separated collection of sequences
        write out the collection as a fasta and then align them using Mafft"""
        for clade in list('ABCDEFGHI'):
            if clade not in self.clade_collection_of_seq_objects_dd:
                continue

            seq_collection = self.clade_collection_of_seq_objects_dd[clade]
            nseqs = len(seq_collection)

            unaligned_fasta_write_out_path = os.path.join(self.outdir, f'its2_clade_{clade}.fasta')

            print(f'Writing {nseqs} clade {clade} seq(s) to {unaligned_fasta_write_out_path}')

            # write out the seq collection
            SeqIO.write(seq_collection, unaligned_fasta_write_out_path, "fasta-2line")

            # if more than one seq in collection, align
            if nseqs > 1:
                print(f'Aligning clade {clade}:')
                aligned_fasta_write_out_path = os.path.join(self.outdir, f'its2_clade_{clade}.aligned.fasta')
                self.aligned_individual_clade_fasta_paths.append(aligned_fasta_write_out_path)
                (self.mafft['--thread', '4','--maxiterate', '1000', '--ep', '0', '--genafpair', unaligned_fasta_write_out_path] > aligned_fasta_write_out_path)()
            else:
                print('only 1 sequence in collection. nothing to align.')
                aligned_fasta_write_out_path = os.path.join(self.outdir, f'its2_clade_{clade}.aligned.fasta')
                self.aligned_individual_clade_fasta_paths.append(aligned_fasta_write_out_path)

    def _create_all_seqs_alignment(self):
        # Create the seeds string that will be passed to the mafft program
        seeds = ""
        for path in self.aligned_individual_clade_fasta_paths:
            seeds += f'--seed {path}'

        # Make the all seqs alignment
        (self.mafft[
             '--thread', '4',
             '--maxiterate', '1000',
             '--ep', '0',
             '--genafpair',
             seeds] > self.all_seqs_aligned_output_path)()




















































up = UnifracPipe()