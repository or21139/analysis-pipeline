# importing libraries
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo

# requires aligned, translated sequences

# an array to use for keys in the dictionary
sample_letters = ["A", "B", "C", "D"]

# fetching the aligned sequences
alignmentA = AlignIO.read("sampleA_aligned.fasta", "fasta")
alignmentB = AlignIO.read("sampleB_aligned.fasta", "fasta")
alignmentC = AlignIO.read("sampleC_aligned.fasta", "fasta")
alignmentD = AlignIO.read("sampleD_aligned.fasta", "fasta")

# creating a dictionary for cleaning
unclean_alignments = {"A": alignmentA, "B": alignmentB,
                      "C": alignmentC, "D": alignmentD}

# clean and format the data as SeqRecords
alignments = {}
for key in sample_letters:
    alignments[key] = []
    for record in unclean_alignments[key]:
        cleaned_seq = record.seq
        parts = record.description.split()
        latin_name = " ".join(parts[1:3])
        genus_species = latin_name.replace(" ", "_")
        if genus_species == "":
            genus_species = f"sample{key}"
        cleaned_record = SeqRecord(
            cleaned_seq,
            id=genus_species,  # id and name are switched for the benefit of the distance matrix
            name=record.id,
            description=record.description
        )
        alignments[key].append(cleaned_record)

# turning cleaned data into an MSA object
for key in sample_letters:
    alignments[key] = MultipleSeqAlignment(alignments[key])

# creating the distance matrix for each sample
calculator = DistanceCalculator("blosum62")
dist_matrices = {}
for key in sample_letters:
    dist_matrices[key] = calculator.get_distance(alignments[key])

# fixing the name and id for the seqrecords
for key in sample_letters:
    for record in alignments[key]:
        old_name = record.name
        old_id = record.id
        record.id = old_name
        record.name = old_id

# defining constructor for future commands
constructor = DistanceTreeConstructor()

# defining outgroup dictionary for trees
outgroups = {"A": "Balaenoptera_musculus", "B": "Stenella_attenuata",
             "C": "Auxis_rochei", "D": "Lepidochelys_kempii"}

# function to build a tree for each sample and output the each tree as a nexus file


def tree_builder(dm):
    for key in sample_letters:
        tree = constructor.nj(dm[key])
        outgroup = tree.find_any(name=outgroups[key])
        tree.root_with_outgroup(outgroup)
        ingroup_terms = [
            t for t in tree.get_terminals()
            if t.name != outgroup
        ]
        ingroup_clade = tree.common_ancestor(ingroup_terms)
        ingroup_clade.name = "Ingroup"
        X = ingroup_clade.branch_length
        ingroup_clade.branch_length = X/2
        outgroup.branch_length = X/2
        Phylo.write(tree, f"tree_{key}.nex", "nexus")


# calling function to construct trees
tree_builder(dist_matrices)
