# importing libraries
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# requires aligned sequences from protein for each sample

# an array to use for keys in the dictionary
sample_letters = ["A", "B", "C", "D"]

# loading sequences from aligned fasta files
sequencesA = list(SeqIO.parse("sampleA_sequences.fasta", "fasta"))
sequencesB = list(SeqIO.parse("sampleB_sequences.fasta", "fasta"))
sequencesC = list(SeqIO.parse("sampleC_sequences.fasta", "fasta"))
sequencesD = list(SeqIO.parse("sampleD_sequences.fasta", "fasta"))
sample_sequences = {"A": sequencesA, "B": sequencesB,
                    "C": sequencesC, "D": sequencesD}

# defining ambiguous base and aa sets
ambiguous_bases = set("NRYSWKMBDHV")
ambiguous_aas = set("XBZJUO")

# function to calculate percentage of ambiguous bases in sequences


def ambiguity_dna_percent(sequences):
    ambig_dna_percent = []
    for sample_seq in sequences:
        num_ambiguous = 0
        for base in sample_seq:
            if base in ambiguous_bases:
                num_ambiguous += 1
        ambig_dna_percent.append(num_ambiguous/len(sample_seq)*100)
    return ambig_dna_percent


# adding ambiguous base percentages to a list
ambig_dna = []
for key in sample_letters:
    ambig_dna.append(ambiguity_dna_percent(sample_sequences[key]))

# function to translate nucleotide to amino acid sequences


def translate_sequences(sequences):
    aa_sequences = []

    # define variables to hold final seq
    for record in sequences:
        nucleo_seq = record.seq
        best_orf = ""
        best_frame = None

        # cycling through each reading frame, find the longest orf
        for frame in range(3):
            protein_sequence = nucleo_seq[frame:].translate(to_stop=False)
            fragments = str(protein_sequence).split("*")
            longest_in_frame = max(fragments, key=len)

            if len(longest_in_frame) > len(best_orf):
                best_orf = longest_in_frame
                best_frame = frame

        # extract the name of the sequence
        parts = record.description.split()
        latin_name = " ".join(parts[1:3])
        genus_species = latin_name.replace(" ", "_")

        # add as a seqrecord
        aa_record = SeqRecord(
            Seq(best_orf),
            id=record.id,
            name=genus_species,
            description=record.description,
            annotations={
                "type": "longest ORF",
                "frame": "best_frame"
            }
        )
        aa_sequences.append(aa_record)

    # output translateed sequences to a file called aa_*sample name*.fasta
    SeqIO.write(aa_sequences, f"aa_{aa_sequences[0].id}.fasta", "fasta")


# calling function to translate sequences
for key in sample_letters:
    translate_sequences(sample_sequences[key])

# function to calculate ambiguous amino acid percentage


def ambiguity_aa_percent(sequences):
    ambig_aa_percent = []
    for sample_seq in sequences:
        num_ambiguous = 0
        for aa in sample_seq:
            if aa in ambiguous_aas:
                num_ambiguous += 1
        ambig_aa_percent.append(num_ambiguous/len(sample_seq)*100)
    return ambig_aa_percent


# defining protein sequences and adding to a dictionary
proteinsA = list(SeqIO.parse("aa_sampleA.fasta", "fasta"))
proteinsB = list(SeqIO.parse("aa_sampleB.fasta", "fasta"))
proteinsC = list(SeqIO.parse("aa_sampleC.fasta", "fasta"))
proteinsD = list(SeqIO.parse("aa_sampleD.fasta", "fasta"))
sample_proteins = {"A": proteinsA, "B": proteinsB,
                   "C": proteinsC, "D": proteinsD}

# calling aa ambiguity function and creating an array
ambig_aa = []
for key in sample_letters:
    ambig_aa.append(ambiguity_aa_percent(sample_proteins[key]))

print(ambig_aa)
print(ambig_dna)
