from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from pathlib import Path
from subprocess import PIPE
import os, sys, argparse, subprocess, re

##  Location of input fasta files
input_dir = "output"

##  Output directory
output_dir = "gene_alignments" 

parser = argparse.ArgumentParser(description = "Align against reference genome (gb) and extract corresponding features")
parser.add_argument("--reference", type = str, help = "File to be used as a reference. Genbank format", required = True)
parser.add_argument("--type", type = str, nargs = '?', help = "Feature type to extract.  Default=product", default = "product")
parser.add_argument("--conserved", help = "Only output sequences conserved in all samples.  Default=False", action = 'store_true')
parser.add_argument("--min_cov", type = float, nargs = '?', help = "Minimum fraction of gene present to accept sequence.  Default=0.95", default = 0.95)
parser.add_argument("--genome", help = "Output whole genome sequence alignment as \"wholegenome.fas\".  Requires extract_sequences.py to have been run with \"--genome\".  Ignores --conserved switch for this file.  Default=False", action = 'store_true')


args = parser.parse_args()

##  Parse genbank file, only retaining first record
ref_record = next(SeqIO.parse (args.reference, "genbank"))

## Get target list and expected length
targets = {}
for feature in ref_record.features:
    if args.type in feature.qualifiers:
        feature_desc = re.sub("\W+", "", str(feature.qualifiers[args.type]))
        feature_length = len(feature.extract(ref_record.seq))
        targets[feature_desc] = feature_length

##  Get all sequence files
dir = Path(input_dir)
fasta_files = list(dir.glob("*.fas"))

##  Open output handles
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

file_handles = {}
for target in targets:
    output_file = target + ".fas"
    output_file_path = os.path.join (output_dir, output_file)
    f = open(output_file_path,'w')
    file_handles[target] = f

##  Extract whole genome alignment if requested
if args.genome == True:
    output_file = "wholegenome.fas"
    output_file_path = os.path.join (output_dir, output_file)
    f = open(output_file_path,'w')
    file_handles["wholegenome"] = f

## Extract sequences
for file in fasta_files:
    sequences = {}
    id = ""
    for record in SeqIO.parse (file, "fasta"):
        ##  Check that description conformed to expected form
        matches = re.match("^(\S+) (\S+)", record.description)
        if matches:
            id = matches.group(1)
            description = matches.group(2)
            ##  Output whole genome alignment if requested
            if description == "wholegenome" and args.genome == True:
                sequence = record.seq.upper()
                print (f'>{id}\n{sequence}', file = file_handles["wholegenome"])
            ##  Make sure sequence is in our desired set
            if description in targets:
                sequence = record.seq.upper()
                ##  Knock out gaps and ambiguous sites to determine high quality sequence
                trimmed_sequence = re.sub("[^ATCG]", "", str(sequence))
                if len(trimmed_sequence) / targets[description] >= args.min_cov:
                    sequences[description]=sequence
                else:
                    print (f'\tRejected {description} from {id}: Match too short')
    ##  Check that file contained all desired targets
    missing_targets = []
    if args.conserved:
        for target in targets:
            if target not in sequences:
                missing_targets.append(target)
    if len(missing_targets) == 0:
        for target in targets:
            if target in sequences:
                seq = str(sequences[target])
                print (f'>{id}\n{seq}', file = file_handles[target]) 
    else :
        output_string = ", ".join(missing_targets)
        print (f'Rejected {id}: Missing {output_string}')


# #  Align query against reference
# if not os.path.exists(alignment_dir):
    # os.makedirs(alignment_dir)
# alignment_file = args.reference + ".aln"
# alignment_file_path = os.path.join (alignment_dir, alignment_file)
# aln_fh = open (alignment_file_path, "w")
# result = subprocess.run(
    # [ "mafft",
      # "--keeplength",
      # "--anysymbol", 
      # "--adjustdirection",
      # "--quiet",
      # "--thread", args.threads,
      # "--add", args.query,
      # temp_file_path, 
    # ],
    # stdout = aln_fh,
    # universal_newlines=True
    # )

# ##  Read alignment
# alignment = SeqIO.parse (alignment_file_path, "fasta")

# ##  Discard reference genome
# next (alignment)

# ##  Extract features
# extracted_records = []
# for record in alignment:
    # if args.genome == True:
        # genome_record = SeqRecord( record.seq,
                                # id = record.id,
                                # description = "wholegenome"
                              # )
        # extracted_records.append(genome_record)
    # for feature in ref_record.features:
        # if args.type in feature.qualifiers:
            # feature_desc = re.sub("\W+", "", str(feature.qualifiers[args.type]))
            # feature_locus_tag = str(feature.qualifiers["locus_tag"])
            # extracted_sequence = feature.extract(record.seq)
            # new_record = SeqRecord( extracted_sequence,
                                    # id = record.id,
                                    # name = feature_locus_tag,
                                    # description = feature_desc
                                  # )
            # extracted_records.append(new_record)

# ## Write output
# if not os.path.exists(output_dir):
    # os.makedirs(output_dir)
# output_file = args.output + ".fas"
# output_file_path = os.path.join (output_dir, output_file)

# SeqIO.write(extracted_records, output_file_path, "fasta-2line")

# ## Clean up temp files unless "retain" is specified
# os.remove(temp_file_path)
# if args.retain == False:
    # os.remove(alignment_file_path)