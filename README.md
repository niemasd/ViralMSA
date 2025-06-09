# ViralMSA
ViralMSA is a tool to perform reference-guided multiple sequence alignment of viral genomes. ViralMSA wraps around existing read mapping tools such as [Minimap2](https://doi.org/10.1093/bioinformatics/bty191), and as such, it can natively improve as methods of read mapping evolve. Importantly, this approach scales linearly with the number of sequences and can be massively parallelized. However, insertions with respect to the reference genome will be thrown away. This is fair for many viral analyses (e.g. phylogenetic inference, as insertions with respect to the reference likely lack phylogenetic information), but it may not be appropriate for all contexts.

To run ViralMSA, you can either install the command-line tool (instructions below), or you can use the [ViralWasm-Epi web app](https://niema-lab.github.io/ViralWasm-Epi) created by my student, [Daniel Ji](https://www.linkedin.com/in/danielji26), which is a complete WebAssembly port of ViralMSA and other downstream analyses, e.g. phylogenetic inference and molecular clustering, which runs fully client-side in your own web browser. The web app works well for reasonably small datasets (e.g. a few thousand full genomes), but for larger datasets, you will want to use the command-line tool.

## Installation
ViralMSA is written in Python 3 and depends on [BioPython](https://biopython.org/). You can simply download [ViralMSA.py](https://github.com/niemasd/ViralMSA/releases/latest/download/ViralMSA.py) to your machine and make it executable:

```bash
wget "https://github.com/niemasd/ViralMSA/releases/latest/download/ViralMSA.py"
chmod a+x ViralMSA.py
sudo mv ViralMSA.py /usr/local/bin/ViralMSA.py # optional step to install globally
```

If you already have ViralMSA installed, you can update it to the newest release version easily:

```bash
ViralMSA.py -u
```

ViralMSA also requires at least one of the following tools to perform the alignment:

* **[Minimap2](https://github.com/lh3/minimap2) (used by default; strongly recommended)**
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [BWA](https://bio-bwa.sourceforge.net/)
* [DRAGMAP](https://github.com/Illumina/DRAGMAP)
* [HISAT2](http://daehwankimlab.github.io/hisat2)
* [LRA](https://github.com/ChaissonLab/LRA)
* [mm2-fast](https://github.com/bwa-mem2/mm2-fast)
* [NGMLR](https://github.com/philres/ngmlr)
* [seq-align](https://github.com/noporpoise/seq-align)
* [STAR](https://github.com/alexdobin/STAR)
* [Unimap](https://github.com/lh3/unimap)
* [wfmash](https://github.com/ekg/wfmash)
* [Winnowmap](https://github.com/marbl/Winnowmap)

We also provide a Docker image with all dependencies installed: [niemasd/viralmsa](https://hub.docker.com/r/niemasd/viralmsa)

## Usage
ViralMSA can be used as follows:

```
usage: ViralMSA.py [-h] -s SEQUENCES -r REFERENCE -e EMAIL -o OUTPUT [-a ALIGNER] [-t THREADS] [-l] [--omit_ref] [--viralmsa_dir VIRALMSA_DIR]

  -h, --help                                  show this help message and exit
  -s SEQUENCES, --sequences SEQUENCES         Input Sequences (FASTA format)
  -r REFERENCE, --reference REFERENCE         Reference
  -o OUTPUT, --output OUTPUT                  Output Directory
  -e EMAIL, --email EMAIL                     Email Address (for Entrez)
  -a ALIGNER, --aligner ALIGNER               Aligner (default: Minimap2)
  -t THREADS, --threads THREADS               Number of Threads (default: max)
  -b BUFFER_SIZE, --buffer_size BUFFER_SIZE   File Stream Buffer Size (bytes) (default: 1048576)
  -l, --list_references                       List all reference sequences (default: False)
  --omit_ref                                  Omit reference sequence from output alignment (default: False)
  --viralmsa_dir VIRALMSA_DIR                 ViralMSA Cache Directory (default: ~/.viralmsa)
  -u, --update                                Update ViralMSA (default: False)
```

## Reference Genome Selection
ViralMSA provides multiple options for selecting a reference genome against which reference-guided multiple sequence alignment will be performed.

### GenBank Accession Number
For the reference genome, you can provide a GenBank accession number, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r NC_045512
```

### Preselected Reference Genomes
For specific viruses of interest, you can simply use their name, and we have provided what we believe would be a good choice of reference genome, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r SARS-CoV-2
```

Our preselected viral reference genomes can be found in the following GitHub repository:

https://github.com/Niema-Lab/Reference-Genomes

If you would like to contribute a new or updated viral reference genome to our collection, please feel free to submit a [GitHub Issue](https://github.com/Niema-Lab/Reference-Genomes/issues/new) or a [Pull Request](https://github.com/Niema-Lab/Reference-Genomes/pulls) with all relevant information.

### Local File
If you have a local reference genome you would like to use, you can provide the path to a FASTA file with a single sequence, such as the following:

```
ViralMSA.py -e email@address.com -s sequences.fas -o output -r my_reference.fas
```

## Citing ViralMSA
If you use ViralMSA in your work, please cite:

> Moshiri N (2021). "ViralMSA: Massively scalable reference-guided multiple sequence alignment of viral genomes." *Bioinformatics*. 37(5):714–716. [doi:10.1093/bioinformatics/btaa743](https://doi.org/10.1093/bioinformatics/btaa743)

If you use ViralMSA via the [ViralWasm-Epi web application](https://niema-lab.github.io/ViralWasm-Epi) (rather than the ViralMSA command-line tool), please *also* cite:

> Ji D, Aboukhalil R, Moshiri N (2023). "ViralWasm: a client-side user-friendly web application suite for viral genomics." *Bioinformatics*. btae018. [doi:10.1093/bioinformatics/btae018](https://doi.org/10.1093/bioinformatics/btae018)

Please also cite the read mapper you selected.

### **Minimap2 (default selection; only option for web app)**

> Li H (2018). "Minimap2: pairwise alignment for nucleotide sequences." *Bioinformatics*. 34(18):3094–3100. [doi:10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)

### bowtie2

> Langmead B, Salzberg SL (2012). "Fast gapped-read alignment with Bowtie 2." *Nat Methods*. 9(4):357-359. [doi:10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)

### BWA

> Li H, Durbin R (2009). "Fast and accurate short read alignment with Burrows–Wheeler transform." *Bioinformatics*. 25(14):1754-1760. [doi:10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324)

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

### seq-align
> Turner I (2015). "seq-align". https://github.com/noporpoise/seq-align

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

## Missing letters from the beginning/end of some sequences
If you notice that some sequences are missing letters at the beginning/end, even though they don't appear as insertions with respect to the reference genome, it could be that the read mapper (e.g. Minimap2) is soft-clipping those bases from the alignment (if you open the `.sam` file in the output folder, those entries should have soft-clipping, denoted by `S`, at the beginning/end of the CIGAR string). A simple hacky fix is to add a bunch of matching letters (e.g. `AAAAA...`) to the beginnings/ends of all of your sequences (including your reference genome) to force a match at the beginning/end of the pairwise alignment done by the read mapper, thus preventing soft-clipping.

## Missing Sequences in Output
Because ViralMSA relies on read mappers to compute each pairwise alignment to the reference genome, and because read mappers can fail to align a sequence if it deviates too significantly from the reference genome, sequences from your input can be omitted from ViralMSA's output if they did not appear in the output of the underlying read mapper you used, and you'll be shown the following warning:

> WARNING: Some sequences from the input are missing from the output. Perhaps try a different aligner or reference genome?

As is mentioned in the warning, to remedy this, you can either try a different reference genome (if you believe there is one that be more similar to *all* of your sequences than the reference you're currently using), or you can try a different aligner (using the `-a` argument). I will briefly mention that there is typically a trade-off between accuracy/sensitivity and speed for different read mappers, so feel free to review the literature to compare/contrast the different read mappers supported by ViralMSA.

As a last resort, ViralMSA also supports [seq-align](https://github.com/noporpoise/seq-align), which performs complete pairwise global alignments. This method is the most sensitive (it's guaranteed to align every sequence in your input), but it will likely be ***much*** slower than the other options. To use it, you can specify `-a seq-align`.
