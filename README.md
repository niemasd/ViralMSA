# ViralMSA
ViralMSA is a tool to perform reference-guided multiple sequence alignment of viral genomes. ViralMSA wraps around existing read mapping tools such as [Minimap2](https://doi.org/10.1093/bioinformatics/bty191), and as such, it can natively improve as methods of read mapping evolve. Importantly, this approach scales linearly with the number of sequences and can be massively parallelized. However, insertions with respect to the reference genome will be thrown away. This is fair for many viral analyses (e.g. phylogenetic inference, as insertions with respect to the reference likely lack phylogenetic information), but it may not be appropriate for all contexts.

Importantly, ViralMSA differs significantly from [VIRULIGN](https://doi.org/10.1093/bioinformatics/bty851), a codon-correct reference-guided alignment tool designed for viruses, in three key ways:

1. ViralMSA in its default settings is orders of magnitude faster than VIRULIGN, which is critical for rapidly-growing epidemics
2. ViralMSA in its default settings uses significantly less memory than does VIRULIGN, which makes it feasible for running analyses on a personal machine
    * On a SARS-CoV-2 whole-genome dataset containing 100 sequences, VIRULIGN consumed roughly 10 GB of memory, whereas ViralMSA consumed hundreds of MB
3. VIRULIGN is codon-aware, making it appropriate for the alignment of coding sequences, whereas ViralMSA is appropriate for the alignment of full viral genomes

## Installation
ViralMSA is written in Python 3 and depends on [BioPython](https://biopython.org/). You can simply download [ViralMSA.py](ViralMSA.py) to your machine and make it executable:

```bash
wget "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py"
chmod a+x ViralMSA.py
sudo mv ViralMSA.py /usr/local/bin/ViralMSA.py # optional step to install globally
```

ViralMSA also requires at least one of the following tools to perform the alignment:

* [Minimap2](https://github.com/lh3/minimap2) (used by default; strongly recommended)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [HISAT2](http://daehwankimlab.github.io/hisat2/)
* [STAR](https://github.com/alexdobin/STAR)

## Usage
ViralMSA can be used as follows:

```
usage: ViralMSA.py [-h] -s SEQUENCES -r REFERENCE -e EMAIL -o OUTPUT [-a ALIGNER] [-t THREADS] [-l] [--omit_ref] [--viralmsa_dir VIRALMSA_DIR]

  -h, --help                                  show this help message and exit
  -s SEQUENCES, --sequences SEQUENCES         Input Sequences (FASTA format)
  -r REFERENCE, --reference REFERENCE         Reference
  -e EMAIL, --email EMAIL                     Email Address (for Entrez)
  -o OUTPUT, --output OUTPUT                  Output Directory
  -a ALIGNER, --aligner ALIGNER               Aligner (default: Minimap2)
  -t THREADS, --threads THREADS               Number of Threads (default: max)
  -b BUFFER_SIZE, --buffer_size BUFFER_SIZE   File Stream Buffer Size (bytes) (default: 1048576)
  -l, --list_references                       List all reference sequences (default: False)
  --omit_ref                                  Omit reference sequence from output alignment (default: False)
  --viralmsa_dir VIRALMSA_DIR                 ViralMSA Cache Directory (default: ~/.viralmsa)
  -u, --update                                Update ViralMSA (default: False)
```

For the reference, you can provide a GenBank accession number, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r MT072688
```

For specific viruses of interest, you can simply use their name, and we have provided what we believe would be a good choice of reference genome, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r SARS-CoV-2
```

If you have a local reference genome you would like to use, you can provide the path to a FASTA file with a single sequence, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r my_reference.fas
```

# Citing ViralMSA
If you use ViralMSA in your work, please cite:

> **Moshiri N** (2020). "ViralMSA: Massively scalable reference-guided multiple sequence alignment of viral genomes." *Bioinformatics*. btaa743. [doi:10.1093/bioinformatics/btaa743](https://doi.org/10.1093/bioinformatics/btaa743)

Please also cite the read mapper you selected. The citation for Minimap2 (the default selection) is the following:

> Li H (2018). "Minimap2: pairwise alignment for nucleotide sequences." *Bioinformatics*. 34(18):3094–3100. [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)

The citation for bowtie2 is the following:

> Langmead B, Salzberg SL (2012). "Fast gapped-read alignment with Bowtie 2." *Nat Methods*. 9(4):357-359. [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)

The citation for HISAT2 is the following:

> Kim D, Paggi JM, Park C, Bennett C, Salzberg SL (2019). "Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype." *Nat Biotechnol*. 37:907-915. [doi:10.1038/s41587-019-0201-4](https://doi.org/10.1038/s41587-019-0201-4)

The citation for STAR is the following:

> Dobin A, Davis CA, Schlesinger F, Drehkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013). "STAR: ultrafast universal RNA-seq aligner." *Bioinformatics*. 29(1):15-21. [doi:10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

# Common Issues
## Weird ViralMSA output on sequences with many `N`s
It seems as though, in some cases in which an input viral sequence has many `N`s within the sequence, Minimap2 splits the input sequence at each long consecutive chain of `N`s and produces an alignment for each fragment of the input sequence, with only one of these alignments (probably the longest one?) being labeled as the primary alignment (flag 0 in the SAM file) and all others being labeled as supplementary alignments (flag 2048 in the SAM file). I'm hoping to find a way to fix this permanently by somehow merging the primary and supplementary alignments into a single alignment (see [this GitHub Issue](https://github.com/lh3/minimap2/issues/720)), but for now, a simple fix that seems to work well is to simply delete all `N`s from the sequence before running ViralMSA.

**TL;DR:** If you are getting a weird ViralMSA output on sequences with many `N`s, try deleting the `N`s (not replacing with `-`: just deleting) before running ViralMSA.
