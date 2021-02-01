from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from pathlib import Path
from subprocess import PIPE
import os, sys, argparse, subprocess, re

temp_dir = "temp"
output_dir = "output"
alignment_dir = "alignments"
# genome_dir = "genomes"

parser = argparse.ArgumentParser(description = "Align against reference genome (gb) and extract corresponding features")
parser.add_argument("--reference", type = str, help = "File to be used as a reference. Genbank format", required = True)
parser.add_argument("--query", type = str, help = "File to extract features from. FASTA format", required = True)
parser.add_argument("--name", type = str, help = "File name for output files.", required = True)
parser.add_argument("--type", type = str, nargs = '?', help = "Feature type to extract.  Default=product", default = "product")
parser.add_argument("--retain", help = "Retain intermediate files.  Default=False", action = 'store_true')
parser.add_argument("--threads", type = str, nargs = '?', help = "Threads to use for mafft.  Default=1", default = "1")
parser.add_argument("--genome", help = "Output whole genome sequence.  Default=False", action = 'store_true')
parser.add_argument("--overwrite", help = "Overwrite existing output files.  Default=False", action = 'store_true')

args = parser.parse_args()

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
output_file = args.name + ".fas"
output_file_path = os.path.join (output_dir, output_file)

##  Do not run and overwrite existing output file unless overwrite == True
if args.overwrite == False and os.path.exists(output_file_path):
    print (f"Output file {output_file_path} already exists.  If you wish to overwrite, use \"--overwrite\".")
    sys.exit(0)

##  Refuse input files with more or less than one genome
query_records = list(SeqIO.parse(args.query, "fasta"))
record_count = len(query_records)
if record_count > 1:
    print (f"Query {args.query} contains more than one record.  Terminating execution.", file=sys.stderr)
    sys.exit(1)
if record_count < 1:
    print (f"Query {args.query} does not contain a record.  Check that input file contains a single record in fasta format.  Terminating execution.", file=sys.stderr)
    sys.exit(1)
    
##  Parse genbank file, only retaining first record
ref_record = next(SeqIO.parse (args.reference, "genbank"))

##  Write out a fasta file with reference sequence
if not os.path.exists(temp_dir):
    os.makedirs(temp_dir)
temp_file = args.reference + ".temp"
temp_file_path = os.path.join (temp_dir, temp_file)
if not os.path.exists(temp_file_path):
    SeqIO.write(ref_record, temp_file_path, "fasta-2line")

#  Align query against reference
if not os.path.exists(alignment_dir):
    os.makedirs(alignment_dir)
alignment_file = args.name + ".aln"
alignment_file_path = os.path.join (alignment_dir, alignment_file)
aln_fh = open (alignment_file_path, "w")
result = subprocess.run(
    [ "mafft",
      "--keeplength",
      "--anysymbol", 
      "--adjustdirection",
      "--quiet",
      "--thread", args.threads,
      "--add", args.query,
      temp_file_path, 
    ],
    stdout = aln_fh,
    universal_newlines=True
    )

##  Read alignment
alignment = SeqIO.parse (alignment_file_path, "fasta")

##  Discard reference genome
next (alignment)

##  Extract features
extracted_records = []
for record in alignment:
    if args.genome == True:
        genome_record = SeqRecord( record.seq,
                                id = record.id,
                                description = "wholegenome"
                              )
        extracted_records.append(genome_record)
    for feature in ref_record.features:
        if args.type in feature.qualifiers:
            feature_desc = re.sub("\W+", "", str(feature.qualifiers[args.type]))
            feature_locus_tag = str(feature.qualifiers["locus_tag"])
            extracted_sequence = feature.extract(record.seq)
            new_record = SeqRecord( extracted_sequence,
                                    id = record.id,
                                    name = feature_locus_tag,
                                    description = feature_desc
                                  )
            extracted_records.append(new_record)

## Write output
SeqIO.write(extracted_records, output_file_path, "fasta-2line")

# Clean up temp files unless "retain" is specified
if args.retain == False:
    os.remove(alignment_file_path)
