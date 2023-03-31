# ViralMSA
ViralMSA is a tool to perform reference-guided multiple sequence alignment of viral genomes. ViralMSA wraps around existing read mapping tools such as [Minimap2](https://doi.org/10.1093/bioinformatics/bty191), and as such, it can natively improve as methods of read mapping evolve. Importantly, this approach scales linearly with the number of sequences and can be massively parallelized. However, insertions with respect to the reference genome will be thrown away. This is fair for many viral analyses (e.g. phylogenetic inference, as insertions with respect to the reference likely lack phylogenetic information), but it may not be appropriate for all contexts.

To run ViralMSA, you can either install the command-line tool (instructions below), or you can try out the [web app](https://niema.net/ViralMSA) created by my student, [Daniel Ji](https://www.linkedin.com/in/danielji26), which is a complete WebAssembly port of ViralMSA (meaning it runs fully client-side in your own web browser!). The web app works well for reasonably small datasets (e.g. a few thousand full genomes), but for larger datasets, you will want to use the command-line tool.

## Installation
ViralMSA is written in Python 3 and depends on [BioPython](https://biopython.org/). You can simply download [ViralMSA.py](ViralMSA.py) to your machine and make it executable:

```bash
wget "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py"
chmod a+x ViralMSA.py
sudo mv ViralMSA.py /usr/local/bin/ViralMSA.py # optional step to install globally
```

ViralMSA also requires at least one of the following tools to perform the alignment:

* **[Minimap2](https://github.com/lh3/minimap2) (used by default; strongly recommended)**
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [DRAGMAP](https://github.com/Illumina/DRAGMAP)
* [HISAT2](http://daehwankimlab.github.io/hisat2)
* [LRA](https://github.com/ChaissonLab/LRA)
* [mm2-fast](https://github.com/bwa-mem2/mm2-fast)
* [NGMLR](https://github.com/philres/ngmlr)
* [STAR](https://github.com/alexdobin/STAR)
* [Unimap](https://github.com/lh3/unimap)
* [wfmash](https://github.com/ekg/wfmash)
* [Winnowmap](https://github.com/marbl/Winnowmap)

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

For the reference genome, you can provide a GenBank accession number, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r NC_045512
```

For specific viruses of interest, you can simply use their name, and we have provided what we believe would be a good choice of reference genome, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r SARS-CoV-2
```

If you have a local reference genome you would like to use, you can provide the path to a FASTA file with a single sequence, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r my_reference.fas
```

## Citing ViralMSA
If you use ViralMSA in your work, please cite:

> Moshiri N (2021). "ViralMSA: Massively scalable reference-guided multiple sequence alignment of viral genomes." *Bioinformatics*. 37(5):714–716. [doi:10.1093/bioinformatics/btaa743](https://doi.org/10.1093/bioinformatics/btaa743)

Please also cite the read mapper you selected.

### **Minimap2 (the default selection)**

> Li H (2018). "Minimap2: pairwise alignment for nucleotide sequences." *Bioinformatics*. 34(18):3094–3100. [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)

### bowtie2

> Langmead B, Salzberg SL (2012). "Fast gapped-read alignment with Bowtie 2." *Nat Methods*. 9(4):357-359. [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)

### DRAGMAP

> https://github.com/Illumina/DRAGMAP

### HISAT2

> Kim D, Paggi JM, Park C, Bennett C, Salzberg SL (2019). "Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype." *Nat Biotechnol*. 37:907-915. [doi:10.1038/s41587-019-0201-4](https://doi.org/10.1038/s41587-019-0201-4)

### LRA

> Ren J, Chaisson MJP (2021). "lra: A long read aligner for sequences and contigs." *PLoS Comput Biol*. 17(6):e1009078. [doi:10.1371/journal.pcbi.1009078](https://doi.org/10.1371/journal.pcbi.1009078)

### mm2-fast

> Kalikar S, Jain C, Vasimuddin M, Misra S (2022). "Accelerating minimap2 for long-read sequencing applications on modern CPUs." *Nat Comput Sci*. 2:78-83. [doi:10.1038/s43588-022-00201-8](https://doi.org/10.1038/s43588-022-00201-8)

### NGMLR

> Sedlazeck FJ, Rescheneder P, Smolka M, Fang H, Nattestad M, von Haeseler A, Schatz MC (2018). "Accurate detection of complex structural variations using single-molecule sequencing." *Nat Methods*. 15:461-468. [doi:10.1038/s41592-018-0001-7](https://doi.org/10.1038/s41592-018-0001-7)

### STAR

> Dobin A, Davis CA, Schlesinger F, Drehkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013). "STAR: ultrafast universal RNA-seq aligner." *Bioinformatics*. 29(1):15-21. [doi:10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

### Unimap

> Li H (2021). "Unimap: A fork of minimap2 optimized for assembly-to-reference alignment." https://github.com/lh3/unimap

### wfmash

> Jain C, Koren S, Dilthey A, Phillippy AM, Aluru S (2018). "A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps". *Bioinformatics*. 34(17):i748-i756. [doi:10.1093/bioinformatics/bty597](https://doi.org/10.1093/bioinformatics/bty597)

> Marco-Sola S, Moure JC, Moreto M, Espinosa A (2021). "Fast gap-affine pairwise alignment using the wavefront algorithm". *Bioinformatics*. 37(4):456-463. [doi:10.1093/bioinformatics/btaa777](https://doi.org/10.1093/bioinformatics/btaa777)

### Winnowmap

> Jain C, Rhie A, Hansen NF, Koren S, Phillippy AM (2022). "Long-read mapping to repetitive reference sequences using Winnowmap2". *Nature Methods*. 19:705-710. [doi:10.1038/s41592-022-01457-8](https://doi.org/10.1038/s41592-022-01457-8)

# Common Issues
## Weird ViralMSA output on sequences with many `N`s
It seems as though, in some cases in which an input viral sequence has many `N`s within the sequence, Minimap2 splits the input sequence at each long consecutive chain of `N`s and produces an alignment for each fragment of the input sequence, with only one of these alignments (probably the longest one?) being labeled as the primary alignment (flag 0 in the SAM file) and all others being labeled as supplementary alignments (flag 2048 in the SAM file). This issue should be fixed in [ViralMSA 1.1.12](https://github.com/niemasd/ViralMSA/releases/tag/1.1.12), but if you run into this issue, a simple fix that seems to work well is to simply shorten the streches of `N`s in the sequence to be at most ~100 before running ViralMSA.
