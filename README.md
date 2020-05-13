# ViralMSA
ViralMSA is a tool to perform reference-guided multiple sequence alignment of viral genomes. ViralMSA wraps around existing read mapping tools such as [Minimap2](https://doi.org/10.1093/bioinformatics/bty191), and as such, it can natively improve as methods of read mapping evolve. Importantly, this approach scales linearly with the number of sequences and can be massively parallelized. However, insertions with respect to the reference genome will be thrown away. This is fair for many viral analyses (e.g. phylogenetic inference, as insertions with respect to the reference likely lack phylogenetic information), but it may not be appropriate for all contexts.

Importantly, ViralMSA differs significantly from [VIRULIGN](https://doi.org/10.1093/bioinformatics/bty851), a codon-correct reference-guided alignment tool designed for viruses, in three key ways:

1. VIRULIGN requires a thorough annotation of the reference genome, which may be difficult to obtain (especially towards the beginning of a novel outbreak) and does not provide the user a simple way to utilize different reference genomes for different viral strains. ViralMSA, on the other hand, only requires the reference genome assembly's GenBank accession number and can build any required index files on-the-fly
2. ViralMSA in its default settings is orders of magnitude faster than VIRULIGN, which is critical for rapidly-growing epidemics
3. VIRULIGN is codon-aware, making it appropriate for the alignment of coding sequences, whereas ViralMSA is appropriate for the alignment of full viral genomes

## Installation
ViralMSA is written in Python 3 and depends on [BioPython](https://biopython.org/). You can simply download [ViralMSA.py](ViralMSA.py) to your machine and make it executable.

ViralMSA also requires at least one of the following tools to perform the alignment:

* [Minimap2](https://github.com/lh3/minimap2) (used by default)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [HISAT2](http://daehwankimlab.github.io/hisat2/)
* [STAR](https://github.com/alexdobin/STAR)

## Usage
ViralMSA can be used as follows:

```
usage: ViralMSA.py [-h] -s SEQUENCES -r REFERENCE -e EMAIL -o OUTPUT [-a ALIGNER] [-t THREADS] [-l] [--omit_ref] [--viralmsa_dir VIRALMSA_DIR]

  -h, --help                            show this help message and exit
  -s SEQUENCES, --sequences SEQUENCES   Input Sequences (FASTA format) (default: None)
  -r REFERENCE, --reference REFERENCE   Reference (default: None)
  -e EMAIL, --email EMAIL               Email Address (for Entrez) (default: None)
  -o OUTPUT, --output OUTPUT            Output Directory (default: None)
  -a ALIGNER, --aligner ALIGNER         Aligner (default: Minimap2)
  -t THREADS, --threads THREADS         Number of Threads (default: max)
  -l, --list_references                 List all reference sequences (default: False)
  --omit_ref                            Omit reference sequence from output alignment (default: False)
  --viralmsa_dir VIRALMSA_DIR           ViralMSA Cache Directory (default: ~/.viralmsa)
```

For the reference, you can provide a GenBank accession number, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r MT072688
```

For specific viruses of interest, you can simply use their name, and we have provided what we believe would be a good choice of reference genome, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r SARS-CoV-2
```

# Citing ViralMSA
If you use ViralMSA in your work, please cite:

> **Moshiri N** (2020). "ViralMSA: Massively scalable reference-guided multiple sequence alignment of viral genomes." *bioRxiv*. [doi:10.1101/2020.04.20.052068](https://doi.org/10.1101/2020.04.20.052068)
