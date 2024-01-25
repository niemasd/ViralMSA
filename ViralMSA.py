#! /usr/bin/env python3
'''
ViralMSA: Reference-guided multiple sequence alignment of viral genomes
'''

# standard imports
from datetime import datetime
from gzip import open as gopen
from hashlib import md5
from io import BufferedReader, TextIOWrapper
from json import load as jload
from math import log2
from multiprocessing import cpu_count, Pool
from os import chdir, getcwd, makedirs, remove
from os.path import abspath, expanduser, isdir, isfile, split
from shutil import copy, move
from urllib.request import urlopen
import argparse
import subprocess
import sys

# useful constants
VERSION = '1.1.42'
RELEASES_URL = 'https://api.github.com/repos/niemasd/ViralMSA/tags'
CIGAR_LETTERS = {'M','D','I','S','H','=','X'}
DEFAULT_BUFSIZE = 1048576 # 1 MB #8192 # 8 KB
DEFAULT_ALIGNER = 'minimap2'
DEFAULT_THREADS = cpu_count()
global QUIET; QUIET = False
global LOGFILE; LOGFILE = None

# mapper-specific options
WINNOWMAP_K = 15 # using Minimap2's default of k=15
WINNOWMAP_DISTINCT = 0.9998

# read mappers that output FASTA
ALIGNERS_FAS = {
    'seq-align',
}

# read mappers that output PAF
ALIGNERS_PAF = {
    'minigraph',
}

# citations
CITATION = {
    'bowtie2':   'Bowtie2: Langmead B, Salzberg SL (2012). "Fast gapped-read alignment with Bowtie 2." Nat Methods. 9(4):357-359. doi:10.1038/nmeth.1923',
    'bwa':       'BWA: Li H, Durbin R (2009). "Fast and accurate short read alignment with Burrows–Wheeler transform." Bioinformatics. 25(14):1754-1760. doi:10.1093/bioinformatics/btp324',
    'dragmap':   'DRAGMAP: https://github.com/Illumina/DRAGMAP',
    'hisat2':    'HISAT2: Kim D, Paggi JM, Park C, Bennett C, Salzberg SL (2019). "Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype." Nat Biotechnol. 37:907-915. doi:10.1038/s41587-019-0201-4',
    'lra':       'LRA: Ren J, Chaisson MJP (2021). "lra: A long read aligner for sequences and contigs." PLoS Comput Biol. 17(6):e1009078. doi:10.1371/journal.pcbi.1009078',
    'minigraph': 'Minigraph: https://github.com/lh3/minigraph',
    'minimap2':  'Minimap2: Li H (2018). "Minimap2: pairwise alignment for nucleotide sequences." Bioinformatics. 34(18):3094-3100. doi:10.1093/bioinformatics/bty191',
    'mm2-fast':  'mm2-fast: Kalikar S, Jain C, Vasimuddin M, Misra S (2022). "Accelerating minimap2 for long-read sequencing applications on modern CPUs." Nat Comput Sci. 2:78-83. doi:10.1038/s43588-022-00201-8',
    'ngmlr':     'NGMLR: Sedlazeck FJ, Rescheneder P, Smolka M, Fang H, Nattestad M, von Haeseler A, Schatz MC (2018). "Accurate detection of complex structural variations using single-molecule sequencing." Nat Methods. 15:461-468. doi:10.1038/s41592-018-0001-7',
    'seq-align': 'seq-align: Turner I (2015). "seq-align". https://github.com/noporpoise/seq-align',
    'star':      'STAR: Dobin A, Davis CA, Schlesinger F, Drehkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013). "STAR: ultrafast universal RNA-seq aligner." Bioinformatics. 29(1):15-21. doi:10.1093/bioinformatics/bts635',
    'unimap':    'Unimap: Li H (2021). "Unimap: A fork of minimap2 optimized for assembly-to-reference alignment." https://github.com/lh3/unimap',
    'viralmsa':  'ViralMSA: Moshiri N (2021). "ViralMSA: Massively scalable reference-guided multiple sequence alignment of viral genomes." Bioinformatics. 37(5):714–716. doi:10.1093/bioinformatics/btaa743',
    'wfmash':    'wfmash: Jain C, Koren S, Dilthey A, Phillippy AM, Aluru S (2018). "A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps". Bioinformatics. 34(17):i748-i756. doi:10.1093/bioinformatics/bty597. Marco-Sola S, Moure JC, Moreto M, Espinosa A (2021). "Fast gap-affine pairwise alignment using the wavefront algorithm". Bioinformatics. 37(4):456-463. doi:10.1093/bioinformatics/btaa777',
    'winnowmap': 'Winnowmap2: Jain C, Rhie A, Hansen NF, Koren S, Phillippy AM (2022). "Long-read mapping to repetitive reference sequences using Winnowmap2". Nature Methods. 19:705-710. doi:10.1038/s41592-022-01457-8',
}

# reference genomes for common viruses
REFS = {
    'bombalivirus':    'NC_039345', # Bombali Virus (Bombali ebolavirus)
    'bundibugyovirus': 'NC_014373', # Bundibugyo Virus (Bundibugyo ebolavirus)
    'chikv':           'NC_004162', # Chikungunya virus
    'denv1':           'NC_001477', # Dengue Virus 1
    'denv2':           'NC_001474', # Dengue Virus 2
    'denv3':           'NC_001475', # Dengue Virus 3
    'denv4':           'NC_002640', # Dengue Virus 4
    'ebolavirus':      'NC_002549', # Ebola Virus (Zaire ebolavirus)
    'hcv1':            'NC_004102', # HCV genotype 1
    'hcv1h77':         'NC_038882', # HCV genotpye 1 (isolate H77)
    'hcv2':            'NC_009823', # HCV genotype 2
    'hcv3':            'NC_009824', # HCV genotype 3
    'hcv4':            'NC_009825', # HCV genotype 4
    'hcv5':            'NC_009826', # HCV genotype 5
    'hcv6':            'NC_009827', # HCV genotype 6
    'hcv7':            'NC_030791', # HCV genotype 7
    'hiv1':            'NC_001802', # HIV-1
    'hiv2':            'NC_001722', # HIV-2
    'measles':         'NC_001498', # Measles Virus (Measles morbillivirus)
    'mpox':            'NC_063383', # Mpox Virus
    'restonvirus':     'NC_004161', # Reston Virus (Reston ebolavirus)
    'sarscov2':        'NC_045512', # SARS-CoV-2 (COVID-19)
    'sudanvirus':      'NC_006432', # Sudan Virus (Sudan ebolavirus)
    'taiforestvirus':  'NC_014372', # Tai Forest Virus (Tai Forest ebolavirus, Cote d'Ivoire ebolavirus)
}

# common names of viruses (for -l listing)
REF_NAMES = {
    'CHIKV': {
        'chikv': 'Chikungunya virus',
    },
    'DENV': {
        'denv1':           'Dengue Virus 1',
        'denv2':           'Dengue Virus 2',
        'denv3':           'Dengue Virus 3',
        'denv4':           'Dengue Virus 4',
    },

    'Ebola': {
        'bombalivirus':    'Bombali Virus (Bombali ebolavirus)',
        'bundibugyovirus': 'Bundibugyo Virus (Bundibugyo ebolavirus)',
        'ebolavirus':      'Ebola Virus (Zaire ebolavirus)',
        'restonvirus':     'Reston Virus (Reston ebolavirus)',
        'sudanvirus':      'Sudan Virus (Sudan ebolavirus)',
        'taiforestvirus':  'Tai Forest Virus (Tai Forest ebolavirus, Cote d\'Ivoire ebolavirus)',
    },

    'HCV': {
        'hcv1':            'HCV genotype 1',
        'hcv1h77':         'HCV genotpye 1 (isolate H77)',
        'hcv2':            'HCV genotype 2',
        'hcv3':            'HCV genotype 3',
        'hcv4':            'HCV genotype 4',
        'hcv5':            'HCV genotype 5',
        'hcv6':            'HCV genotype 6',
        'hcv7':            'HCV genotype 7',
    },
    
    'HIV': {
        'hiv1':            'HIV-1',
        'hiv2':            'HIV-2',
    },

    'Measles': {
        'measles':         'Measles Virus',
    },

    'Mpox': {
        'monkeypox':       'Mpox Virus',
        'mpox':            'Mpox Virus',
    },
    
    'SARS-CoV-2': {
        'covid19':         'SARS-CoV-2 (COVID-19)',
        'sarscov2':        'SARS-CoV-2 (COVID-19)',
        'sc2':             'SARS-CoV-2 (COVID-19)',
    }
}

# aliases for user-friendliness
REF_ALIASES = {
    'covid19':   'sarscov2',
    'monkeypox': 'mpox',
    'sc2':       'sarscov2',
}
for k,v in REF_ALIASES.items():
    REFS[k] = REFS[v]
    for virus in REF_NAMES:
        if v in REF_NAMES[virus]:
            REF_NAMES[virus][k] = REF_NAMES[virus][v]
            break

# check for validity in reference names
tmp = set(REFS.keys()) - {k2 for k1 in REF_NAMES for k2 in REF_NAMES[k1]}
assert len(tmp) == 0, "Value(s) in REFS missing in REF_NAMES: %s" % str(tmp)

# print to log (prefixed by current time)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    if not QUIET:
        print(tmp, file=sys.stderr, end=end); sys.stderr.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()

# convert a ViralMSA version string to a tuple of integers
def parse_version(s):
    return tuple(int(v) for v in s.split('.'))

# update ViralMSA to the newest version
def update_viralmsa():
    tags = jload(urlopen(RELEASES_URL))
    newest = max(tags, key=lambda x: parse_version(x['name']))
    if parse_version(newest['name']) <= parse_version(VERSION):
        print("ViralMSA is already at the newest version (%s)" % VERSION); exit(0)
    old_version = VERSION; new_version = newest['name']
    url = 'https://raw.githubusercontent.com/niemasd/ViralMSA/%s/ViralMSA.py' % newest['commit']['sha']
    content = urlopen(url).read()
    try:
        with open(__file__, 'wb') as f:
            f.write(content)
    except PermissionError:
        print("ERROR: Received a permission error when updating ViralMSA. Perhaps try running as root?", file=sys.stderr); exit(1)
    print("Successfully updated ViralMSA %s --> %s" % (old_version, new_version)); exit(0)

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# count the number of IDs in a FASTA file
def count_IDs_fasta(fn, bufsize=DEFAULT_BUFSIZE):
    return sum(l.startswith('>') for l in open(fn, buffering=bufsize))

# count the number of IDs in a SAM file
def count_IDs_sam(fn, bufsize=DEFAULT_BUFSIZE):
    return len({l.split('\t')[0].strip() for l in open(fn, buffering=bufsize) if not l.startswith('@')})

# parse a CIGAR string
def parse_cigar(s):
    out = list(); ind = len(s)-1
    while ind >= 0:
        let = s[ind]; ind -= 1; num = ''
        while s[ind] not in CIGAR_LETTERS:
            num += s[ind]; ind -= 1
        out.append((let, int(num[::-1])))
    return out[::-1]

# iterate over a FASTA file as (header, seq) tuples
def iter_fasta(fa_path, bufsize=DEFAULT_BUFSIZE):
    if fa_path.lower().endswith('.gz'):
        fa_file = TextIOWrapper(BufferedReader(gopen(fa_path, 'rb'), buffer_size=bufsize))
    else:
        fa_file = open(fa_path, 'r', buffering=bufsize)
    header = None; seq = None; line = None
    for line in fa_file:
        l = line.strip()
        if l.startswith('>'):
            if seq is not None:
                if len(seq) == 0:
                    raise ValueError("Invalid FASTA file: %s" % fa_path)
                else:
                    yield (header, seq)
            header = l[1:]; seq = ''
        elif len(l) != 0:
            seq += l.replace(' ','').replace('\t','')
    if len(seq) == 0:
        raise ValueError("Invalid FASTA file: %s" % fa_path)
    yield (header, seq)

# convert FASTA to FASTQ
def fasta2fastq(fa_path, fq_path, qual='~', bufsize=DEFAULT_BUFSIZE):
    if fq_path.lower().endswith('.gz'):
        gzip_out = True; fq_file = gopen(fq_path, 'wb', 9)
    else:
        gzip_out = False; fq_file = open(fq_path, 'w', buffering=bufsize)
    for fa_header, fa_seq in iter_fasta(fa_path):
        curr = '@%s\n%s\n+\n%s\n' % (fa_header, fa_seq, qual*len(fa_seq))
        if gzip_out:
            fq_file.write(curr.encode())
        else:
            fq_file.write(curr)

# check bowtie2
def check_bowtie2():
    try:
        o = subprocess.check_output(['bowtie2', '-h'])
    except:
        o = None
    if o is None or 'Bowtie 2 version' not in o.decode():
        print("ERROR: bowtie2 is not runnable in your PATH", file=sys.stderr); exit(1)
    try:
        o = subprocess.check_output(['bowtie2-build', '-h'])
    except:
        o = None
    if o is None or 'Bowtie 2 version' not in o.decode():
        print("ERROR: bowtie2-build is not runnable in your PATH", file=sys.stderr); exit(1)

# check BWA
def check_bwa():
    try:
        o = subprocess.run(['bwa'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
    except:
        o = None
    if o is None or 'Program: bwa' not in o.decode():
        print("ERROR: BWA is not runnable in your PATH", file=sys.stderr); exit(1)

# check DRAGMAP
def check_dragmap():
    try:
        o = subprocess.check_output(['dragen-os', '-h'])
    except:
        o = None
    if o is None or 'dragenos -r <reference> -b <base calls> [optional arguments]' not in o.decode():
        print("ERROR: dragen-os is not runnable in your PATH", file=sys.stderr); exit(1)

# check HISAT2
def check_hisat2():
    try:
        o = subprocess.check_output(['hisat2', '-h'])
    except:
        o = None
    if o is None or 'HISAT2 version' not in o.decode():
        print("ERROR: hisat2 is not runnable in your PATH", file=sys.stderr); exit(1)
    try:
        o = subprocess.check_output(['hisat2-build', '-h'])
    except:
        o = None
    if o is None or 'HISAT2 version' not in o.decode():
        print("ERROR: hisat2-build is not runnable in your PATH", file=sys.stderr); exit(1)

# check LRA
def check_lra():
    try:
        o = subprocess.check_output(['lra', '-h'])
    except subprocess.CalledProcessError as cpe:
        o = cpe.output
    except:
        o = None
    if o is None or 'lra (long sequence alignment)' not in o.decode():
        print(o.decode())
        print("ERROR: LRA is not runnable in your PATH", file=sys.stderr); exit(1)

# check minigraph
def check_minigraph():
    try:
        o = subprocess.run(['minigraph'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
    except:
        o = None
    if o is None or 'Usage: minigraph' not in o.decode():
        print("ERROR: Minigraph is not runnable in your PATH", file=sys.stderr); exit(1)

# check minimap2
def check_minimap2():
    try:
        o = subprocess.check_output(['minimap2', '-h'])
    except:
        o = None
    if o is None or 'Usage: minimap2' not in o.decode():
        print("ERROR: Minimap2 is not runnable in your PATH", file=sys.stderr); exit(1)

# check mm2-fast
def check_mm2fast():
    try:
        o = subprocess.check_output(['mm2-fast', '-h'])
    except:
        o = None
    if o is None or 'Usage: ' not in o.decode():
        print("ERROR: mm2-fast is not runnable in your PATH", file=sys.stderr); exit(1)

# check NGMLR
def check_ngmlr():
    try:
        o = subprocess.check_output(['ngmlr', '-h'], stderr=subprocess.STDOUT)
    except:
        o = None
    if o is None or 'Usage: ngmlr' not in o.decode():
        print("ERROR: NGMLR is not runnable in your PATH", file=sys.stderr); exit(1)

# check seq-align
def check_seqalign():
    try:
        o = subprocess.run(['needleman_wunsch', '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
    except:
        o = None
    if o is None or 'usage: needleman_wunsch' not in o.decode():
        print("ERROR: seq-align (needleman_wunsch) is not runnable in your PATH", file=sys.stderr); exit(1)

# check STAR
def check_star():
    try:
        o = subprocess.check_output(['STAR', '-h'])
    except:
        o = None
    if o is None or 'Usage: STAR' not in o.decode():
        print("ERROR: STAR is not runnable in your PATH", file=sys.stderr); exit(1)

# check unimap
def check_unimap():
    try:
        o = subprocess.check_output(['unimap', '-h'])
    except:
        o = None
    if o is None or 'Usage: unimap' not in o.decode():
        print("ERROR: Unimap is not runnable in your PATH", file=sys.stderr); exit(1)

# check wfmash
def check_wfmash():
    try:
        o = subprocess.check_output(['wfmash', '-h'])
    except:
        o = None
    if o is None or 'wfmash [target] [queries...] {OPTIONS}' not in o.decode():
        print("ERROR: wfmash is not runnable in your PATH", file=sys.stderr); exit(1)

# check Winnowmap
def check_winnowmap():
    try:
        o = subprocess.check_output(['winnowmap', '-h'])
    except:
        o = None
    if o is None or 'Usage: winnowmap' not in o.decode():
        print("ERROR: Winnowmap is not runnable in your PATH", file=sys.stderr); exit(1)
    try:
        o = subprocess.run(['meryl'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
    except:
        o = None
    if o is None or 'usage: meryl' not in o.decode():
        print("ERROR: meryl is not runnable in your PATH", file=sys.stderr); exit(1)

# build FAIDX index (multiple tools might need this)
def build_index_faidx(ref_genome_path, threads, verbose=True):
    index_path = '%s.fai' % ref_genome_path
    if isfile(index_path):
        if verbose:
            print_log("FAIDX index found: %s" % index_path)
        return
    command = ['samtools', 'faidx', ref_genome_path]
    if verbose:
        print_log("Building FAIDX index: %s" % ' '.join(command))
    log = open('%s.log' % index_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("FAIDX index built: %s" % index_path)

# build bowtie2 index
def build_index_bowtie2(ref_genome_path, threads, verbose=True):
    exts = ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2']
    all_found = True
    for ext in exts:
        if not isfile('%s.bowtie2.%s' % (ref_genome_path,ext)):
            all_found = False; break
    if all_found:
        if verbose:
            print_log("bowtie2 index found: %s.bowtie2.*.bt2" % ref_genome_path)
        return
    for ext in exts:
        if isfile('%s.bowtie2.%s' % (ref_genome_path,ext)):
            remove('%s.bowtie2.%s' % (ref_genome_path,ext))
    ref_genome_dir, ref_genome_fn = split(ref_genome_path)
    orig_dir = getcwd()
    chdir(ref_genome_dir)
    command = ['bowtie2-build', '--threads', str(threads), ref_genome_path, '%s.bowtie2' % ref_genome_fn]
    if verbose:
        print_log("Building bowtie2 index: %s" % ' '.join(command))
    log = open('%s.bowtie2.log' % ref_genome_path, 'w'); subprocess.call(command, stdout=log, stderr=subprocess.STDOUT); log.close()
    if verbose:
        print_log("bowtie2 index built: %s.bowtie2.*.bt2" % ref_genome_path)
    chdir(orig_dir)

# build BWA index
def build_index_bwa(ref_genome_path, threads, verbose=True):
    exts = ['amb', 'ann', 'bwt', 'pac', 'sa']
    all_found = True
    for ext in exts:
        if not isfile('%s.%s' % (ref_genome_path,ext)):
            all_found = False; break
    if all_found:
        if verbose:
            print_log("BWA index found: %s.[amb, ann, bwt, pac, sa]" % ref_genome_path)
        return
    for ext in exts:
        if isfile('%s.%s' % (ref_genome_path,ext)):
            remove('%s.%s' % (ref_genome_path,ext))
    ref_genome_dir, ref_genome_fn = split(ref_genome_path)
    orig_dir = getcwd()
    chdir(ref_genome_dir)
    command = ['bwa', 'index', ref_genome_path]
    if verbose:
        print_log("Building BWA index: %s" % ' '.join(command))
    log = open('%s.bwa.log' % ref_genome_path, 'w'); subprocess.call(command, stdout=log, stderr=subprocess.STDOUT); log.close()
    if verbose:
        print_log("BWA index built: %s.[amb, ann, bwt, pac, sa]" % ref_genome_path)
    chdir(orig_dir)

# build DRAGMAP index
def build_index_dragmap(ref_genome_path, threads, verbose=True):
    index_path = '%s.DRAGMAP' % ref_genome_path
    if isdir(index_path):
        if verbose:
            print_log("DRAGMAP index found: %s" % index_path)
        return
    makedirs(index_path, exist_ok=False)
    command = ['dragen-os', '--build-hash-table', 'true', '--ht-reference', ref_genome_path, '--ht-num-threads', str(threads), '--output-directory', index_path]
    if verbose:
        print_log("Building DRAGMAP index: %s" % ' '.join(command))
    log = open('%s/index.log' % index_path, 'w'); subprocess.call(command, stdout=log, stderr=subprocess.STDOUT); log.close()
    if verbose:
        print_log("DRAGMAP index built: %s" % index_path)

# build HISAT2 index
def build_index_hisat2(ref_genome_path, threads, verbose=True):
    exts = [('%d.ht2'%i) for i in range(1,9)]
    all_found = True
    for ext in exts:
        if not isfile('%s.hisat2.%s' % (ref_genome_path,ext)):
            all_found = False; break
    if all_found:
        if verbose:
            print_log("HISAT2 index found: %s.hisat2.*.ht2" % ref_genome_path)
        return
    for ext in exts:
        if isfile('%s.hisat2.%s' % (ref_genome_path,ext)):
            remove('%s.hisat2.%s' % (ref_genome_path,ext))
    ref_genome_dir, ref_genome_fn = split(ref_genome_path)
    orig_dir = getcwd()
    chdir(ref_genome_dir)
    command = ['hisat2-build', '--threads', str(threads), ref_genome_path, '%s.hisat2' % ref_genome_fn]
    if verbose:
        print_log("Building HISAT2 index: %s" % ' '.join(command))
    log = open('%s.hisat2.log' % ref_genome_path, 'w'); subprocess.call(command, stdout=log, stderr=subprocess.STDOUT); log.close()
    if verbose:
        print_log("HISAT2 index built: %s.hisat2.*.ht2" % ref_genome_path)
    chdir(orig_dir)

# build LRA index
def build_index_lra(ref_genome_path, threads, verbose=True):
    lra_ref_genome_path = '%s.lra' % ref_genome_path
    gli_index_path = '%s.gli' % lra_ref_genome_path
    mmi_index_path = '%s.mmi' % lra_ref_genome_path
    mms_index_path = '%s.mms' % lra_ref_genome_path
    if isfile(lra_ref_genome_path) and isfile(gli_index_path) and (isfile(mmi_index_path) or isfile(mms_index_path)):
        if verbose:
            if isfile(mmi_index_path):
                print_log("LRA index files found: %s and %s" % (gli_index_path, mmi_index_path))
            else:
                print_log("LRA index files found: %s and %s" % (gli_index_path, mms_index_path))
        return
    elif isfile(lra_ref_genome_path):
        raise RuntimeError("Corrupt LRA index. Please delete the following and try again: %s" % lra_ref_genome_path)
    elif isfile(gli_index_path):
        raise RuntimeError("Corrupt LRA index. Please delete the following and try again: %s" % gli_index_path)
    elif isfile(mmi_index_path):
        raise RuntimeError("Corrupt LRA index. Please delete the following and try again: %s" % mmi_index_path)
    elif isfile(mms_index_path):
        raise RuntimeError("Corrupt LRA index. Please delete the following and try again: %s" % mms_index_path)
    copy(ref_genome_path, lra_ref_genome_path)
    command = ['lra', 'index', '-CONTIG', lra_ref_genome_path]
    if verbose:
        print_log("Building LRA index: %s" % ' '.join(command))
    log = open('%s.log' % lra_ref_genome_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        if isfile(mmi_index_path):
            print_log("LRA index built: %s and %s" % (gli_index_path, mmi_index_path))
        else:
            print_log("LRA index built: %s and %s" % (gli_index_path, mms_index_path))

# build minigraph index
def build_index_minigraph(ref_genome_path, threads, verbose=True):
    pass # Minigraph doesn't seem to do anything with a ref genome with 1 sequence

# build minimap2 index
def build_index_minimap2(ref_genome_path, threads, verbose=True):
    index_path = '%s.mmi' % ref_genome_path
    if isfile(index_path):
        if verbose:
            print_log("Minimap2 index found: %s" % index_path)
        return
    command = ['minimap2', '-t', str(threads), '-d', index_path, ref_genome_path]
    if verbose:
        print_log("Building Minimap2 index: %s" % ' '.join(command))
    log = open('%s.log' % index_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("Minimap2 index built: %s" % index_path)

# build mm2-fast index
def build_index_mm2fast(ref_genome_path, threads, verbose=True):
    index_path = '%s.mm2fast.mmi' % ref_genome_path
    if isfile(index_path):
        if verbose:
            print_log("mm2-fast index found: %s" % index_path)
        return
    command = ['mm2-fast', '-t', str(threads), '-d', index_path, ref_genome_path]
    if verbose:
        print_log("Building mm2-fast index: %s" % ' '.join(command))
    log = open('%s.log' % index_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("mm2-fast index built: %s" % index_path)

# build NGMLR index
def build_index_ngmlr(ref_genome_path, threads, verbose=True):
    enc_index_path = '%s-enc.2.ngm' % ref_genome_path
    ht_index_path = '%s-ht-13-2.2.ngm' % ref_genome_path
    if isfile(enc_index_path) and isfile(ht_index_path):
        if verbose:
            print_log("NGMLR index files found: %s and %s" % (enc_index_path, ht_index_path))
        return
    elif isfile(enc_index_path):
        raise RuntimeError("Corrupt NGMLR index. Please delete the following and try again: %s" % enc_index_path)
    elif isfile(ht_index_path):
        raise RuntimeError("Corrupt NGMLR index. Please delete the following and try again: %s" % ht_index_path)
    command = ['ngmlr', '-x', 'pacbio', '-i', '0', '--no-smallinv', '-t', str(threads), '-r', ref_genome_path]
    if verbose:
        print_log("Building NGMLR index: %s" % ' '.join(command))
    log = open('%s.NGMLR.log' % ref_genome_path, 'w'); subprocess.call(command, stdin=subprocess.DEVNULL, stderr=log, stdout=DEVNULL); log.close()
    if verbose:
        print_log("NGMLR index built: %s and %s" % (enc_index_path, ht_index_path))

# build seq-align index
def build_index_seqalign(ref_genome_path, threads, verbose=True):
    pass # no indexing needed

# build STAR index
def build_index_star(ref_genome_path, threads, verbose=True):
    delete_log = True # STAR by default creates Log.out in the running directory
    if isfile('Log.out'):
        delete_log = False # don't delete it (it existed before running ViralMSA)
    index_path = '%s.STAR' % ref_genome_path
    if isdir(index_path):
        if verbose:
            print_log("STAR index found: %s" % index_path)
        return
    genome_length = sum(len(l.strip()) for l in open(ref_genome_path) if not l.startswith('>'))
    genomeSAindexNbases = min(14, int(log2(genome_length)/2)-1)
    makedirs(index_path, exist_ok=False)
    command = ['STAR', '--runMode', 'genomeGenerate', '--runThreadN', str(threads), '--genomeDir', index_path, '--genomeFastaFiles', ref_genome_path, '--genomeSAindexNbases', str(genomeSAindexNbases)]
    if verbose:
        print_log("Building STAR index: %s" % ' '.join(command))
    log = open('%s/index.log' % index_path, 'w'); subprocess.call(command, stdout=log, stderr=subprocess.STDOUT); log.close()
    if verbose:
        print_log("STAR index built: %s" % index_path)
    if delete_log and isfile('Log.out'):
        remove('Log.out')

# build unimap index
def build_index_unimap(ref_genome_path, threads, verbose=True):
    index_path = '%s.umi' % ref_genome_path
    if isfile(index_path):
        if verbose:
            print_log("Unimap index found: %s" % index_path)
        return
    command = ['unimap', '-t', str(threads), '-d', index_path, ref_genome_path]
    if verbose:
        print_log("Building Unimap index: %s" % ' '.join(command))
    log = open('%s.log' % index_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("Unimap index built: %s" % index_path)

# build wfmash index
def build_index_wfmash(ref_genome_path, threads, verbose=True):
    build_index_faidx(ref_genome_path, threads, verbose=True)

# build Winnowmap index
def build_index_winnowmap(ref_genome_path, threads, verbose=True):
    db_path = '%s.meryl.db' % ref_genome_path
    index_path = '%s.meryl.k%d.txt' % (ref_genome_path, WINNOWMAP_K)
    if isdir(db_path):
        if verbose:
            print_log("Meryl DB found: %s" % db_path)
    else:
        command = ['meryl', 'count', 'k=%d' % WINNOWMAP_K, 'output', db_path, ref_genome_path]
        if verbose:
            print_log("Building Meryl DB: %s" % ' '.join(command))
        log = open('%s.log' % db_path, 'w'); subprocess.call(command, stderr=log); log.close()
        if verbose:
            print_log("Meryl DB built: %s" % db_path)
    if isfile(index_path):
        if verbose:
            print_log("Winnowmap index (Meryl %d-mers) found: %s" % (WINNOWMAP_K, index_path))
        return
    command = ['meryl', 'print', 'greater-than', 'distinct=%s' % WINNOWMAP_DISTINCT, db_path]
    if verbose:
        print_log("Building Winnowmap index (Meryl %d-mers): %s" % (WINNOWMAP_K, ' '.join(command)))
    out = open(index_path, 'w'); log = open('%s.log' % index_path, 'w'); subprocess.call(command, stdout=out, stderr=log); out.close(); log.close()
    if verbose:
        print_log("Winnowmap index (Meryl %d-mers) built: %s" % (WINNOWMAP_K, index_path))

# align genomes using bowtie2
def align_bowtie2(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    command = ['bowtie2', '--very-sensitive', '-p', str(threads), '-f', '-x', '%s.bowtie2' % ref_genome_path, '-U', seqs_path, '-S', out_aln_path]
    if verbose:
        print_log("Aligning using bowtie2: %s" % ' '.join(command))
    log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("bowtie2 alignment complete: %s" % out_aln_path)

# align genomes using BWA MEM
def align_bwa(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    command = ['bwa', 'mem', '-t', str(threads), ref_genome_path, seqs_path]
    if verbose:
        print_log("Aligning using BWA: %s" % ' '.join(command))
    out = open(out_aln_path, 'w'); log = open('%s.log' % out_aln_path, 'w')
    subprocess.call(command, stdout=out, stderr=log); out.close(); log.close()
    if verbose:
        print_log("BWA alignment complete: %s" % out_aln_path)

# align genomes using DRAGMAP
def align_dragmap(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    tmp_fq_path = '%s.fastq.gz' % out_aln_path
    fasta2fastq(seqs_path, tmp_fq_path)
    command = ['dragen-os', '--num-threads', str(threads), '--Aligner.sw-all', '1', '-r', '%s.DRAGMAP' % ref_genome_path, '-1', tmp_fq_path]
    if verbose:
        print_log("Aligning using DRAGMAP: %s" % ' '.join(command))
    out = open(out_aln_path, 'w'); log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stdout=out, stderr=log); out.close(); log.close()
    if verbose:
        print_log("DRAGMAP alignment complete: %s" % out_aln_path)

# align genomes using HISAT2
def align_hisat2(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    command = ['hisat2', '--very-sensitive', '-p', str(threads), '-f', '-x', '%s.hisat2' % ref_genome_path, '-U', seqs_path, '-S', out_aln_path]
    if verbose:
        print_log("Aligning using HISAT2: %s" % ' '.join(command))
    log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("HISAT2 alignment complete: %s" % out_aln_path)

# align genomes using LRA
def align_lra(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    lra_ref_genome_path = '%s.lra' % ref_genome_path
    command = ['lra', 'align', '-t', str(threads), '-CONTIG', '-p', 's', lra_ref_genome_path, seqs_path]
    if verbose:
        print_log("Aligning using LRA: %s" % ' '.join(command))
    out_sam_file = open(out_aln_path, 'w'); log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stdout=out_sam_file, stderr=log); out_sam_file.close(); log.close()
    if verbose:
        print_log("LRA alignment complete: %s" % out_aln_path)

# align genomes using minigraph
def align_minigraph(seqs_path, out_paf_path, ref_genome_path, threads, verbose=True):
    command = ['minigraph', '-c', '-t', str(threads), '--secondary=no', '-l', '0', '-d', '0', '-L', '0', '-o', out_paf_path, ref_genome_path, seqs_path]
    if verbose:
        print_log("Aligning using Minigraph: %s" % ' '.join(command))
    log = open('%s.log' % out_paf_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("Minigraph alignment complete: %s" % out_paf_path)

# align genomes using minimap2
def align_minimap2(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    index_path = '%s.mmi' % ref_genome_path
    command = ['minimap2', '-t', str(threads), '--score-N=0', '--secondary=no', '--sam-hit-only', '-a', '-o', out_aln_path, index_path, seqs_path]
    if verbose:
        print_log("Aligning using Minimap2: %s" % ' '.join(command))
    log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("Minimap2 alignment complete: %s" % out_aln_path)

# align genomes using mm2-fast
def align_mm2fast(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    index_path = '%s.mm2fast.mmi' % ref_genome_path
    command = ['mm2-fast', '-t', str(threads), '--score-N=0', '--secondary=no', '--sam-hit-only', '-a', '-o', out_aln_path, index_path, seqs_path]
    if verbose:
        print_log("Aligning using mm2-fast: %s" % ' '.join(command))
    log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("mm2-fast alignment complete: %s" % out_aln_path)

# align genomes using NGMLR
def align_ngmlr(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    command = ['ngmlr', '--skip-write', '-x', 'pacbio', '-i', '0', '--no-smallinv', '-t', str(threads), '-r', ref_genome_path, '-q', seqs_path, '-o', out_aln_path]
    if verbose:
        print_log("Aligning using NGMLR: %s" % ' '.join(command))
    log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("NGMLR alignment complete: %s" % out_aln_path)

# helper function to perform individual pairwise alignments
def run_needleman_wunsch(x): # x is (command, ref_seq, query_header, query_seq) tuple
    ref_aln, query_aln = subprocess.check_output(x[0] + [x[1], x[3]]).decode().strip().splitlines()
    start_ind_in = 0; end_ind_ex = len(ref_aln)
    for c in ref_aln:
        if c == '-':
            start_ind_in += 1
        else:
            break
    for c in ref_aln[::-1]:
        if c == '-':
            end_ind_ex -= 1
        else:
            break
    return (x[2], query_aln[start_ind_in:end_ind_ex])

# align genomes using seq-align
def align_seqalign(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    # set things up
    command = ['needleman_wunsch', '--printfasta', '--nogapsin1']
    if verbose:
        print_log("Aligning using seq-align: %s" % ' '.join(command))
    log = open('%s.log' % out_aln_path, 'w'); out_aln = open(out_aln_path, 'w', buffering=bufsize)
    ref_header, ref_seq = list(iter_fasta(ref_genome_path))[0]

    # perform alignment
    pool = Pool(processes=threads)
    inputs = ((command, ref_seq, curr_header, curr_seq) for curr_header, curr_seq in iter_fasta(seqs_path))
    for curr in pool.imap_unordered(run_needleman_wunsch, inputs, chunksize=1):
        out_aln.write('>%s\n%s\n' % curr)
        log.write("Finished %s\n" % curr[0])
    log.close(); out_aln.close()
    if verbose:
        print_log("seq-align complete: %s" % out_aln_path)

# align genomes using STAR
def align_star(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    delete_log = True # STAR by default creates Log.out in the running directory
    if isfile('Log.out'):
        delete_log = False # don't delete it (it existed before running ViralMSA)
    index_path = '%s.STAR' % ref_genome_path
    out_sam_dir, out_sam_fn = split(out_aln_path)
    out_file_prefix = '%s.' % '.'.join(out_aln_path.split('.')[:-1])
    command = ['STAR', '--runThreadN', str(threads), '--genomeDir', index_path, '--readFilesIn', seqs_path, '--outFileNamePrefix', out_file_prefix, '--outFilterMismatchNmax', '9999999999']
    if verbose:
        print_log("Aligning using STAR: %s" % ' '.join(command))
    log = open('%s/STAR.log' % out_sam_dir, 'w'); subprocess.call(command, stdout=log); log.close()
    move('%sAligned.out.sam' % out_file_prefix, out_aln_path)
    if verbose:
        print_log("STAR alignment complete: %s" % out_sam_dir)
    if delete_log and isfile('Log.out'):
        remove('Log.out')

# align genomes using unimap
def align_unimap(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    index_path = '%s.umi' % ref_genome_path
    command = ['unimap', '-t', str(threads), '--score-N=0', '--secondary=no', '--sam-hit-only', '-a', '--cs', '-o', out_aln_path, index_path, seqs_path]
    if verbose:
        print_log("Aligning using Unimap: %s" % ' '.join(command))
    log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("Unimap alignment complete: %s" % out_aln_path)

# align genomes using wfmash
def align_wfmash(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    ref_genome_length = sum(len(l.strip()) for l in open(ref_genome_path) if not l.startswith('>'))
    command = ['wfmash', '-N', '--sam-format', '--threads=%d' % threads, ref_genome_path, seqs_path]
    if verbose:
        print_log("Aligning using wfmash: %s" % ' '.join(command))
    sam = open(out_aln_path, 'w'); log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stdout=sam, stderr=log); sam.close(); log.close()
    if verbose:
        print_log("wfmash alignment complete: %s" % out_aln_path)

# align genomes using Winnowmap
def align_winnowmap(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    index_path = '%s.meryl.k%d.txt' % (ref_genome_path, WINNOWMAP_K)
    command = ['winnowmap', '-k', str(WINNOWMAP_K), '-W', index_path, '-t', str(threads), '--score-N=0', '--secondary=no', '--sam-hit-only', '-a', '-o', out_aln_path, ref_genome_path, seqs_path]
    if verbose:
        print_log("Aligning using Winnowmap: %s" % ' '.join(command))
    log = open('%s.log' % out_aln_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("Winnowmap alignment complete: %s" % out_aln_path)

# aligners
ALIGNERS = {
    'bowtie2': {
        'check':       check_bowtie2,
        'build_index': build_index_bowtie2,
        'align':       align_bowtie2,
    },

    'bwa': {
        'check':       check_bwa,
        'build_index': build_index_bwa,
        'align':       align_bwa,
    },

    'dragmap': {
        'check':       check_dragmap,
        'build_index': build_index_dragmap,
        'align':       align_dragmap,
    },

    'hisat2': {
        'check':       check_hisat2,
        'build_index': build_index_hisat2,
        'align':       align_hisat2,
    },

    'lra': {
        'check':       check_lra,
        'build_index': build_index_lra,
        'align':       align_lra,
    },

    # Minigraph still doesn't work (only outputs PAF without sequences)
    #'minigraph': {
    #    'check':       check_minigraph,
    #    'build_index': build_index_minigraph,
    #    'align':       align_minigraph,
    #},

    'minimap2': {
        'check':       check_minimap2,
        'build_index': build_index_minimap2,
        'align':       align_minimap2,
    },

    'mm2-fast': {
        'check':       check_mm2fast,
        'build_index': build_index_mm2fast,
        'align':       align_mm2fast,
    },

    'ngmlr': {
        'check':       check_ngmlr,
        'build_index': build_index_ngmlr,
        'align':       align_ngmlr,
    },

    'seq-align': {
        'check':       check_seqalign,
        'build_index': build_index_seqalign,
        'align':       align_seqalign,
    },

    'star': {
        'check':       check_star,
        'build_index': build_index_star,
        'align':       align_star,
    },

    'unimap': {
        'check':       check_unimap,
        'build_index': build_index_unimap,
        'align':       align_unimap,
    },

    'wfmash': {
        'check':       check_wfmash,
        'build_index': build_index_wfmash,
        'align':       align_wfmash,
    },

    'winnowmap': {
        'check':       check_winnowmap,
        'build_index': build_index_winnowmap,
        'align':       align_winnowmap,
    },
}

# handle GUI (updates argv)
def run_gui():
    def clear_argv():
        tmp = sys.argv[0]; sys.argv.clear(); sys.argv.append(tmp)
    try:
        # imports
        from tkinter import Button, Checkbutton, END, Entry, Frame, IntVar, Label, OptionMenu, StringVar, Tk
        from tkinter.filedialog import askdirectory, askopenfilename

        # helper function to make a popup
        def gui_popup(message, title=None):
            popup = Tk()
            if title:
                popup.wm_title(title)
            label = Label(popup, text=message)
            label.pack()
            button_close = Button(popup, text="Close", command=popup.destroy)
            button_close.pack(padx=3, pady=3)
            popup.mainloop()

        # create applet
        root = Tk()
        root.geometry("600x400")
        #root.configure(background='white')

        # set up main frame
        frame = Frame(root)
        frame.pack()

        # add header
        header = Label(frame, text="ViralMSA %s" % VERSION, font=('Arial',24))
        header.pack()

        # handle input FASTA selection
        button_seqs_prefix = "Input Sequences:\n"
        button_seqs_nofile = "<none selected>"
        def find_filename_seqs():
            fn = askopenfilename(title="Select Input Sequences (FASTA format)",filetypes=(("FASTA Files",("*.fasta","*.fas")),("All Files","*.*")))
            if len(fn) == 0:
                button_seqs.configure(text="%s%s" % (button_seqs_prefix,button_seqs_nofile))
            else:
                button_seqs.configure(text="%s%s" % (button_seqs_prefix,fn))
        button_seqs = Button(frame, text="%s%s" % (button_seqs_prefix,button_seqs_nofile), command=find_filename_seqs)
        button_seqs.pack(padx=3, pady=3)

        # handle reference
        dropdown_ref_default = "Select Reference"
        dropdown_ref_var = StringVar(frame)
        dropdown_ref_var.set(dropdown_ref_default)
        dropdown_ref = OptionMenu(frame, dropdown_ref_var, *sorted(REF_NAMES[k][v] for k in REF_NAMES for v in REF_NAMES[k]))
        dropdown_ref.pack()

        # handle user email
        entry_email_default = "Enter Email Address"
        entry_email = Entry(frame, width=30)
        entry_email.insert(END, entry_email_default)
        entry_email.pack()

        # handle output folder selection
        button_out_prefix = "Output Directory:\n"
        button_out_nofolder = "<none selected>"
        def find_directory_out():
            dn = askdirectory(title="Select Output Directory")
            if len(dn) == 0:
                button_out.configure(text="%s%s" % (button_out_prefix,button_out_nofolder))
            else:
                button_out.configure(text="%s%s" % (button_out_prefix,dn))
        button_out = Button(frame, text="%s%s" % (button_out_prefix,button_out_nofolder), command=find_directory_out)
        button_out.pack(padx=3, pady=3)

        # handle aligner selection
        dropdown_aligner_prefix = "Aligner: "
        dropdown_aligner_var = StringVar(frame)
        dropdown_aligner_var.set("%s%s" % (dropdown_aligner_prefix,DEFAULT_ALIGNER.lower()))
        dropdown_aligner = OptionMenu(frame, dropdown_aligner_var, *sorted(("%s%s" % (dropdown_aligner_prefix,a)) for a in ALIGNERS))
        dropdown_aligner.pack()

        # handle threads selection
        dropdown_threads_prefix = "Threads: "
        dropdown_threads_var = StringVar(frame)
        dropdown_threads_var.set("%s%s" % (dropdown_threads_prefix,DEFAULT_THREADS))
        dropdown_threads = OptionMenu(frame, dropdown_threads_var, *[("%s%d" % (dropdown_threads_prefix,i)) for i in range(1,DEFAULT_THREADS+1)])
        dropdown_threads.pack()

        # handle omit reference toggle
        check_omitref_var = IntVar(frame)
        check_omitref = Checkbutton(frame, text="Omit reference sequence from output MSA", variable=check_omitref_var, onvalue=1, offvalue=0)
        check_omitref.pack()

        # add run button
        def finish_applet():
            valid = True
            # check sequences
            try:
                if button_seqs['text'] == "%s%s" % (button_seqs_prefix,button_seqs_nofile):
                    gui_popup("ERROR: Input Sequences file not selected", title="ERROR"); valid = False
            except:
                pass
            # check reference
            try:
                if dropdown_ref_var.get() == dropdown_ref_default:
                   gui_popup("ERROR: Reference not selected", title="ERROR"); valid = False
            except:
                pass
            # check email
            try:
                if '@' not in entry_email.get():
                    gui_popup("ERROR: Email Address not entered", title="ERROR"); valid = False
            except:
                pass
            # check output directory
            try:
                if button_out['text'] == "%s%s" % (button_out_prefix,button_out_nofolder):
                    gui_popup("ERROR: Output Directory not selected", title="ERROR"); valid = False
                elif isdir(button_out['text'].lstrip(button_out_prefix).strip()):
                    gui_popup("ERROR: Output Directory already exists", title="ERROR"); valid = False
            except:
                pass
            # close applet to run ViralMSA
            if valid:
                sys.argv.append('-s'); sys.argv.append(button_seqs['text'].lstrip(button_seqs_prefix).strip())
                sys.argv.append('-r'); sys.argv.append([v for k in REF_NAMES for v in REF_NAMES[k] if REF_NAMES[k][v] == dropdown_ref_var.get()][0])
                sys.argv.append('-e'); sys.argv.append(entry_email.get())
                sys.argv.append('-o'); sys.argv.append(button_out['text'].lstrip(button_out_prefix).strip())
                sys.argv.append('-a'); sys.argv.append(dropdown_aligner_var.get().lstrip(dropdown_aligner_prefix).strip())
                sys.argv.append('-t'); sys.argv.append(dropdown_threads_var.get().lstrip(dropdown_threads_prefix).strip())
                if check_omitref_var.get() == 1:
                    sys.argv.append('--omit_ref')
                try:
                    root.destroy()
                except:
                    pass
        button_run = Button(frame, text="Run", command=finish_applet)
        button_run.pack(padx=3, pady=3)

        # add title and execute GUI
        root.title("ViralMSA %s" % VERSION)
        root.mainloop()
    except:
        print("ERROR: Unable to import Tkinter", file=sys.stderr); exit(1)
    if len(sys.argv) == 1:
        exit()

# parse user args
def parse_args():
    # check if user wants to run the GUI
    if len(sys.argv) == 1:
        run_gui()

    # check if user wants to update ViralMSA
    if '-u' in sys.argv or '--update' in sys.argv:
        update_viralmsa()

    # check if user just wants to list references
    if '-l' in sys.argv or '--list_references' in sys.argv:
        print("=== List of ViralMSA Reference Sequences ===")
        for v in sorted(REF_NAMES.keys()):
            print("* %s" % v)
            for r in sorted(REF_NAMES[v].keys()):
                print("  - %s: %s" % (r, REF_NAMES[v][r]))
        exit(0)

    # use argparse to parse user arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--sequences', required=True, type=str, help="Input Sequences (FASTA format, or SAM if already mapped)")
    parser.add_argument('-r', '--reference', required=True, type=str, help="Reference")
    parser.add_argument('-e', '--email', required=True, type=str, help="Email Address (for Entrez)")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
    parser.add_argument('-a', '--aligner', required=False, type=str, default=DEFAULT_ALIGNER, help="Aligner (options: %s)" % ', '.join(sorted(ALIGNERS.keys())))
    parser.add_argument('-t', '--threads', required=False, type=int, default=DEFAULT_THREADS, help="Number of Threads")
    parser.add_argument('-b', '--buffer_size', required=False, type=int, default=DEFAULT_BUFSIZE, help="File Stream Buffer Size (bytes)")
    parser.add_argument('-l', '--list_references', action="store_true", help="List all reference sequences")
    parser.add_argument('--omit_ref', action="store_true", help="Omit reference sequence from output alignment")
    parser.add_argument('--stdout', action="store_true", help="Write MSA to standard output instead of to file")
    parser.add_argument('-q', '--quiet', action="store_true", help="Suppress log output")
    parser.add_argument('--viralmsa_dir', required=False, type=str, default=abspath(expanduser("~/.viralmsa")), help="ViralMSA Cache Directory")
    parser.add_argument('-u', '--update', action="store_true", help="Update ViralMSA (current version: %s)" % VERSION)
    args = parser.parse_args()

    # check user args for validity
    if args.threads < 1:
        print("ERROR: Number of threads must be positive", file=sys.stderr); exit(1)
    if args.buffer_size < 1:
        print("ERROR: Output buffer size must be positive", file=sys.stderr); exit(1)
    args.aligner = args.aligner.lower()
    if args.aligner not in ALIGNERS:
        print("ERROR: Invalid aligner: %s (valid options: %s)" % (args.aligner, ', '.join(sorted(ALIGNERS.keys()))), file=sys.stderr); exit(1)
    args.sequences = abspath(expanduser(args.sequences))
    if not isfile(args.sequences):
        print("ERROR: Sequences file not found: %s" % args.sequences, file=sys.stderr); exit(1)
    if args.sequences.lower().endswith('.gz'):
        print("ERROR: Sequences cannot be compressed: %s" % args.sequences, file=sys.stderr); exit(1)
    args.output = abspath(expanduser(args.output))
    if isdir(args.output) or isfile(args.output):
        print("ERROR: Output directory exists: %s" % args.output, file=sys.stderr); exit(1)
    if isfile(args.reference):
        if count_IDs_fasta(args.reference, bufsize=args.buffer_size) != 1:
            print("ERROR: Reference file (%s) must have exactly 1 sequence in the FASTA format" % args.reference, file=sys.stderr); exit(1)
        ref_seq = ''.join(l.strip() for l in open(args.reference, 'r', args.buffer_size) if not l.startswith('>'))
        h = md5(ref_seq.encode()).hexdigest()
        fn = args.reference
        args.reference = '%s_HASH_%s' % (fn.split('/')[-1].strip(), h)
        args.ref_path = '%s/%s' % (args.viralmsa_dir, args.reference)
        args.ref_genome_path = '%s/%s' % (args.ref_path, fn.split('/')[-1].strip())
        if not isdir(args.ref_path):
            makedirs(args.ref_path)
            copy(fn, args.ref_genome_path)
    else:
        tmp = args.reference.lower().replace(' ','').replace('-','').replace('_','')
        if tmp in REFS:
            args.reference = REFS[tmp]
        args.reference = args.reference.upper()
        args.ref_path = '%s/%s' % (args.viralmsa_dir, args.reference)
        args.ref_genome_path = '%s/%s.fas' % (args.ref_path, args.reference)

    # user args are valid, so return
    return args

# download reference genome
def download_ref_genome(reference, ref_path, ref_genome_path, email, bufsize=DEFAULT_BUFSIZE):
    try:
        from Bio import Entrez
    except ModuleNotFoundError:
        print("ERROR: Unable to import Biopython, which is needed to download a reference genome.\nInstall with: pip install biopython", file=sys.stderr); exit(1)
    makedirs(ref_path, exist_ok=True); Entrez.email = email
    try:
        handle = Entrez.efetch(db='nucleotide', rettype='fasta', id=reference)
    except:
        raise RuntimeError("Encountered error when trying to download reference genome from NCBI. Perhaps the accession number is invalid?")
    seq = handle.read()
    if seq.count('>') != 1:
        print("ERROR: Reference genome must only have a single sequence", file=sys.stderr); exit(1)
    f = open(ref_genome_path, 'w', buffering=bufsize); f.write(seq.strip()); f.write('\n'); f.close()

# convert alignment (SAM/FASTA/PAF) to FASTA
def aln_to_fasta(out_aln_path, out_msa_path, ref_genome_path, omit_ref=False, bufsize=DEFAULT_BUFSIZE):
    if out_aln_path.lower().endswith('.sam'):
        aln_type = 'SAM'
    elif out_aln_path.lower().endswith('.fas'):
        aln_type = 'FASTA'
    elif out_aln_path.lower().endswith('.paf'):
        aln_type = 'PAF'
    else:
        print("ERROR: Invalid alignment extension: %s" % out_aln_path, file=sys.stderr); exit(1)
    if out_msa_path is None:
        msa = sys.stdout
    else:
        msa = open(out_msa_path, 'w', buffering=bufsize)
    ref_header, ref_seq = list(iter_fasta(ref_genome_path))[0]; ref_seq_len = len(ref_seq)
    if not omit_ref:
        msa.write('>%s\n%s\n' % (ref_header, ref_seq))
    num_output_IDs = 0
    if aln_type == 'SAM':
        for l in open(out_aln_path):
            if l == '\n' or l[0] == '@':
                continue
            parts = l.split('\t')
            flags = int(parts[1])
            if flags != 0 and flags != 16:
                continue
            ID = parts[0].strip(); num_output_IDs += 1
            ref_ind = int(parts[3])-1
            seq = parts[9].strip()
            cigar = parts[5].strip()
            edits = parse_cigar(cigar)
            msa.write(">%s\n" % ID)
            if ref_ind > 0:
                msa.write('-'*ref_ind) # write gaps before alignment
            ind = 0; seq_len = ref_ind
            for e, e_len in edits:
                if e == 'M' or e == '=' or e == 'X': # (mis)match)
                    msa.write(seq[ind:ind+e_len])
                    ind += e_len; seq_len += e_len
                elif e == 'D':                       # deletion (gap in query)
                    msa.write('-'*e_len)
                    seq_len += e_len
                elif e == 'I':                       # insertion (gap in reference; ignore)
                    ind += e_len
                elif e == 'S' or e == 'H':           # starting/ending segment of query not in reference (i.e., span of insertions; ignore)
                    ind += e_len
            if seq_len < ref_seq_len:
                msa.write('-'*(ref_seq_len-seq_len)) # write gaps after alignment
            msa.write('\n')
    elif aln_type == 'FASTA':
        for query in iter_fasta(out_aln_path):
            msa.write('>%s\n%s\n' % query); num_output_IDs += 1
    elif aln_type == 'PAF':
        raise RuntimeError("PAF alignments are not yet supported")
    else:
        raise ValueError("Unknown alignment type: %s" % aln_type)
    msa.close()
    return num_output_IDs

# main content
def main():
    # parse user args and prepare run
    INPUT_TYPE = None
    args = parse_args()
    makedirs(args.viralmsa_dir, exist_ok=True)
    makedirs(args.output)
    global QUIET; QUIET = args.quiet
    global LOGFILE; LOGFILE = open("%s/viralmsa.log" % args.output, 'w')
    if args.sequences.lower().endswith('.sam'):
        INPUT_TYPE = 'SAM'
        num_input_IDs = count_IDs_sam(args.sequences, bufsize=args.buffer_size)
    else: # assume FASTA input if not SAM
        INPUT_TYPE = 'FASTA'
        ALIGNERS[args.aligner]['check']()
        num_input_IDs = count_IDs_fasta(args.sequences, bufsize=args.buffer_size)

    # print run information
    print_log("===== RUN INFORMATION =====")
    print_log("ViralMSA Version: %s" % VERSION)
    print_log("Sequences: %s" % args.sequences)
    print_log("- %d sequences in input file" % num_input_IDs)
    print_log("Reference: %s" % args.reference)
    print_log("Email Address: %s" % args.email)
    print_log("Output Directory: %s" % args.output)
    if INPUT_TYPE == 'FASTA':
        print_log("Aligner: %s" % args.aligner)
    print_log("ViralMSA Cache Directory: %s" % args.viralmsa_dir)
    print_log()

    # download reference genome if not already downloaded
    print_log("===== REFERENCE GENOME =====")
    if isfile(args.ref_genome_path):
        print_log("Reference genome found: %s" % args.ref_genome_path)
    else:
        print_log("Downloading reference genome from NCBI...")
        download_ref_genome(args.reference, args.ref_path, args.ref_genome_path, args.email, bufsize=args.buffer_size)
        print_log("Reference genome downloaded: %s" % args.ref_genome_path)

    # align if FASTA input, otherwise use SAM input as-is
    if INPUT_TYPE == 'FASTA':
        # build aligner index (if needed)
        ALIGNERS[args.aligner]['build_index'](args.ref_genome_path, args.threads)
        print_log()

        # align viral genomes against referencea
        print_log("===== ALIGNMENT =====")
        if args.aligner in ALIGNERS_FAS:
            out_aln_path = '%s/%s.aln.fas' % (args.output, args.sequences.split('/')[-1])
        elif args.aligner in ALIGNERS_PAF:
            out_aln_path = '%s/%s.paf' % (args.output, args.sequences.split('/')[-1])
        else:
            out_aln_path = '%s/%s.sam' % (args.output, args.sequences.split('/')[-1])
        ALIGNERS[args.aligner]['align'](args.sequences, out_aln_path, args.ref_genome_path, args.threads, bufsize=args.buffer_size)
    elif INPUT_TYPE == 'SAM':
        out_aln_path = args.sequences

    # convert alignment (SAM/PAF) to MSA FASTA
    print_log("Converting alignment to FASTA...")
    if args.stdout:
        out_msa_path = None
    else:
        out_msa_path = '%s/%s.aln' % (args.output, args.sequences.split('/')[-1])
    num_output_IDs = aln_to_fasta(out_aln_path, out_msa_path, args.ref_genome_path, omit_ref=args.omit_ref, bufsize=args.buffer_size)
    print_log("Multiple sequence alignment complete: %s" % out_msa_path)
    if num_output_IDs < num_input_IDs:
        print_log("WARNING: Some sequences from the input are missing from the output. Perhaps try a different aligner or reference genome?")
        print_log("- Input: %d sequences" % num_input_IDs)
        print_log("- Output: %d sequences" % num_output_IDs)
    print_log()

    # print citations and finish
    print_log("===== CITATIONS =====")
    print_log(CITATION['viralmsa'])
    print_log(CITATION[args.aligner])
    LOGFILE.close()

# run tool
if __name__ == "__main__":
    main()
