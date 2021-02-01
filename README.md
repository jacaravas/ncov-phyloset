#ncov-phyloset

##SYNOPSIS

This tool set is designed to roughly annotate highly conserved genomes and generate pseudo-alignments of the called gene sequences. Additionally it includes the option to generate conserved core gene alignments containing only isolates present in all analyzed genomes.

While numerous annotation tools and pre-annotated genomes already exist for many taxa, this tool was necessitated by the unusual nature of SARS-CoV-2 genome. Individual annotation of the nsp peptides within the ORF1AB transcript is not typically performed in public repositories. The highly conserved nature of the genome allows for reasonably accurate annotation using a simple pairwise alignment against a reference genome with these features annotated. With the rapid addition of new genome sequences of varying quality, this tool is also designed with filtering capabilities and a to rapidly add new genome sequences to your alignment files.

There are two main components to this pipeline:
1. `extract_sequences.py`, which creates an output file for each genome with the requested sequences.  This step can be parallelized by the user using their own code or the provided Sun Grid Engine subission scripts `run_jobs.pl` and/or `extract_seqs_array.sh`
2. `build_gene_files.py`, which uses fasta files produced by the previous step to build custom alignment sets.

##USAGE
####`extract_sequences.py`
```
usage: extract_sequences.py [-h] --reference REFERENCE --query QUERY --name
                            NAME [--type [TYPE]] [--retain]
                            [--threads [THREADS]] [--genome] [--overwrite]

Align against reference genome (gb) and extract corresponding features

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE File to be used as a reference. Genbank format
  --query QUERY         File to extract features from. FASTA format
  --name NAME           File name for output files.
  --type [TYPE]         Feature type to extract. Default=product
  --retain              Retain intermediate files. Default=False
  --threads [THREADS]   Threads to use for mafft. Default=1
  --genome              Output whole genome sequence. Default=False
  --overwrite           Overwrite existing output files. Default=False
```
Output files will be written as fasta files to a directory paradoxically named "input"
The `--type` argument specifies which genbank tag to look for to include.  Only the presence of the tag is evaluated, not the value.  The default type "product" has been tested, but others have not.  The name of a sequence in the resulting fasta files will be the value associated with the provided tag.  Spaces, punctuation, and most non-alphanumeric chracters are removed from the tag value.
Example tags are highlighted in bold below:
```
     CDS             join(266..13468,13468..21555)
                     /**gene**="ORF1ab"
                     /**locus_tag**="GU280_gp01"
                     /**ribosomal_slippage**
                     /**note**="pp1ab; translated by -1 ribosomal frameshift"
                     /**codon_start**=1
                     /**product**="ORF1ab polyprotein"
                     /**protein_id**="YP_009724389.1"
                     /**db_xref**="GI:1796318597"
                     /**db_xref**="GeneID:437405780"
```

####`build_gene_files.py`
```
usage: build_gene_files.py [-h] --reference REFERENCE [--type [TYPE]]
                           [--conserved] [--min_cov [MIN_COV]] [--genome]

Align against reference genome (gb) and extract corresponding features

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE
                        File to be used as a reference. Genbank format
  --type [TYPE]         Feature type to extract. Default=product
  --conserved           Only output sequences conserved in all samples.
                        Default=False
  --min_cov [MIN_COV]   Minimum fraction of gene present to accept sequence.
                        Default=0.95
  --genome              Output whole genome sequence alignment as
                        "wholegenome.fas". Requires extract_sequences.py to
                        have been run with "--genome". Ignores --conserved
                        switch for this file. Default=False
```
Output files are written to the "gene_alignments" directory.
This function always overwrites existing gene files, so be sure to copy the files another directory if you wish to preserve existing files.
The `--conserved` flag specifies whether all sequences that meet the `--min_cov` value are included (default) or if only genomes which contained all targets at the desired cutoff are included (flag activated).
The `--min_cov` value controls how much of the sequence must be present to consider the sequence present. This value is the fraction of non-ambiguous base calls in the query sequence to the length annotated on the provided reference.

##DATA SOURCES
This program was designed to work with nucleotide FASTA sequence data downloaded from the [NCBI SARS-CoV-2 Resources](https://www.ncbi.nlm.nih.gov/sars-cov-2/) page.  It should be compatible with other single genome fasta format files.
The NCBI reference genome `NC_045512.2.gb` is included in the repository for your convenience.

##CAVEATS
All sites which would result in an insertion relative to the reference sequence are ommitted from the result.
The resulting alignment files are not true multiple alignments.  They are all aligned relative to the reference genome only.


##DEPENDENCIES
Python3 with Biopython.
mafft aligner in path.

##Disclaimer
This software is provided for research purposes only.