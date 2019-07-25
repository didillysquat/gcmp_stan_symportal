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
        self.iqtree = local['iqtree']
        # This is the list that will hold the path to the aligned fastas of each given clade
        self.aligned_individual_clade_fasta_paths = []
        # This list will hold a list of seq collections where each seq collection contains only one sequence
        # of a given clade. This will be used when running mafft.
        self.isolated_clade_collections = []
        self.clade_collection_of_seq_objects_dd = self._populate_clade_collection_of_seq_objects_dd()
        # self.sorted_clade_abund_tup = self._populate_sorted_clade_abund_tup()
        self._align_individual_clade_seq_collections()
        # For the time being I am going to move away from the all seqs alignment.
        # I don't think that it is a good idea to work across clades.
        # self._create_all_seqs_alignment()
        # self._generate_tree_for_each_alignment()
        self.type_profile_output_path = '/Users/humebc/Google_Drive/symportal_outputs/ryan/data_analysis/72_DBV_20190609_2019-06-10_01-30-10.537772.profiles.relative.txt'
        self.profile_to_clade_dict = {}
        self.profile_seq_abundance_df = self._process_profile_output()
        self.do_clade_unifrac()

    def do_clade_unifrac(self):
        # first get a df for each clade
        for clade in list('A'):
            clade_df = self.profile_seq_abundance_df.loc[[tup[0] for tup in self.profile_to_clade_dict.items() if tup[1] == clade]]
            clade_df = clade_df.loc[:,(clade_df != 0).any(axis=0)]
            divTree = skbio.TreeNode.read('/Users/humebc/Google_Drive/projects/mcMinds/gcmp_stan_symportal/phylogenetics/output/unifrac_testing/its2_clade_A.aligned.fasta.treefile').root_at_midpoint()
            divTree.write('/Users/humebc/Google_Drive/projects/mcMinds/gcmp_stan_symportal/phylogenetics/output/unifrac_testing/its2_clade_A.aligned.fasta.rooted.treefile')
            wu = skbio.diversity.beta_diversity('weighted_unifrac', clade_df.to_numpy(), ids=list(clade_df.index), tree=divTree, otu_ids=list(clade_df.columns))
            fooo = 'asdf'

    def _process_profile_output(self):
        import re
        df = pd.read_csv(self.type_profile_output_path, sep='\t', index_col=0)
        df.drop(columns=[df.columns[0]], inplace=True)

        self.profile_to_clade_dict = {prof_uid: clade for prof_uid, clade in zip(list(df), df.loc['Clade'].values.tolist())}
        abunds = {}
        for string, profdef, names in zip(
                df.loc['Average defining sequence proportions and [stdev]'], df.loc['ITS2 type profile'], df.columns):
            abunds[names] = {}
            for substr, div in zip(string.split('-'), profdef.replace('/', '-').split('-')):
                abunds[names][div] = float(re.sub('\[.*\]', '', substr)) * 1000  # need integer values for some reason! provided decimal values were at precision of 1e-3
        df = pd.DataFrame.from_dict(abunds, orient='index')
        df[pd.isna(df)] = 0
        return df

    # def _generate_tree_for_each_alignment(self):
    #     for alignment_path in self.aligned_individual_clade_fasta_paths:
    #         (self.iqtree[
    #              '-s', alignment_path, '-nt', 'AUTO', '-ntmax', '8',] > aligned_fasta_write_out_path)()


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
                if not os.path.exists(aligned_fasta_write_out_path):
                    (self.mafft['--thread', '4','--maxiterate', '1000', '--ep', '0', '--genafpair', unaligned_fasta_write_out_path] > aligned_fasta_write_out_path)()
            else:
                print('only 1 sequence in collection. nothing to align.')
                self.isolated_clade_collections.append(seq_collection)


    def _create_all_seqs_alignment(self, redo=False):
        if not os.path.exists(self.all_seqs_aligned_output_path) or redo:
            # If we are only working with --seed arguments then we need to add this '/dev/null' string at the end
            # However, if we are appending a fasta that contains the isolated sequences, then we do not need it
            terminating_str = '/dev/null'

            # Create the seeds string that will be passed to the mafft program
            seeds = []
            for path in self.aligned_individual_clade_fasta_paths:
                seeds.extend(['--seed', f'{path}'])

            # Now check to see if there are any isolated seq collections
            # if there are then these need to be concatenated together and written out
            if self.isolated_clade_collections:
                if len(self.isolated_clade_collections) == 1:
                    self.isolated_clade_collections = self.isolated_clade_collections[0]
                print('Writing isolated clades to a file')
                terminating_str = os.path.join(self.outdir, 'its2_isolated_concat.fasta')
                SeqIO.write(self.isolated_clade_collections, terminating_str, 'fasta-2line')

            # Make the all seqs alignment
            (self.mafft[
                 '--thread', '4',
                 '--maxiterate', '1000',
                 '--ep', '0',
                 '--genafpair',
                 seeds, terminating_str] > self.all_seqs_aligned_output_path)()





















































up = UnifracPipe()