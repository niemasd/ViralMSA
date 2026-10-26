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
VERSION = '1.1.51'
RELEASES_URL = 'https://api.github.com/repos/niemasd/ViralMSA/tags'
CIGAR_LETTERS = {'M','D','I','S','H','=','X'}
DEFAULT_BUFSIZE = 1048576 # 1 MB #8192 # 8 KB
DEFAULT_ALIGNER = 'rammap'
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
    'rammap':    'rammap: Wang JR, Li H (2026). "Memory-safe high-performance sequence mapping with rammap." bioRxiv. doi:10.64898/2026.05.26.726289',
    'seq-align': 'seq-align: Turner I (2015). "seq-align". https://github.com/noporpoise/seq-align',
    'star':      'STAR: Dobin A, Davis CA, Schlesinger F, Drehkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR (2013). "STAR: ultrafast universal RNA-seq aligner." Bioinformatics. 29(1):15-21. doi:10.1093/bioinformatics/bts635',
    'unimap':    'Unimap: Li H (2021). "Unimap: A fork of minimap2 optimized for assembly-to-reference alignment." https://github.com/lh3/unimap',
    'viralmsa':  'ViralMSA: Moshiri N (2021). "ViralMSA: Massively scalable reference-guided multiple sequence alignment of viral genomes." Bioinformatics. 37(5):714–716. doi:10.1093/bioinformatics/btaa743',
    'wfmash':    'wfmash: Jain C, Koren S, Dilthey A, Phillippy AM, Aluru S (2018). "A Fast Adaptive Algorithm for Computing Whole-Genome Homology Maps". Bioinformatics. 34(17):i748-i756. doi:10.1093/bioinformatics/bty597. Marco-Sola S, Moure JC, Moreto M, Espinosa A (2021). "Fast gap-affine pairwise alignment using the wavefront algorithm". Bioinformatics. 37(4):456-463. doi:10.1093/bioinformatics/btaa777',
    'winnowmap': 'Winnowmap2: Jain C, Rhie A, Hansen NF, Koren S, Phillippy AM (2022). "Long-read mapping to repetitive reference sequences using Winnowmap2". Nature Methods. 19:705-710. doi:10.1038/s41592-022-01457-8',
}

# reference genomes for common viruses
REFS_FAIL = False; REFS_JSON = None; REFS = None; REFS_JSON_URL = 'https://github.com/Niema-Lab/Reference-Genomes/releases/latest/download/REFS.json'
def get_refs_json():
    global REFS_FAIL
    if not REFS_FAIL:
        try:
            global REFS_JSON; REFS_JSON = jload(urlopen(REFS_JSON_URL))
            global REFS; REFS = {n:k for k in REFS_JSON for n in REFS_JSON[k]['shortname']}
        except:
            REFS_FAIL = True

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
    if hasattr(fa_path, 'read') or not isinstance(fa_path, str):
        fa_file = fa_path
    elif fa_path.lower().endswith('.gz'):
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

# check rammap
def check_rammap():
    try:
        o = subprocess.check_output(['rammap', '-h'])
    except:
        o = None
    if o is None or 'Usage: rammap' not in o.decode():
        print("ERROR: rammap is not runnable in your PATH", file=sys.stderr); exit(1)

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
    if o is None or 'wfmash' not in o.decode():
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
    index_path = '%s.minimap2.mmi' % ref_genome_path
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
    log = open('%s.NGMLR.log' % ref_genome_path, 'w'); subprocess.call(command, stdin=subprocess.DEVNULL, stderr=log, stdout=subprocess.DEVNULL); log.close()
    if verbose:
        print_log("NGMLR index built: %s and %s" % (enc_index_path, ht_index_path))

# build rammap index
def build_index_rammap(ref_genome_path, threads, verbose=True):
    index_path = '%s.rammap.mmi' % ref_genome_path
    if isfile(index_path):
        if verbose:
            print_log("rammap index found: %s" % index_path)
        return
    command = ['rammap', '-t', str(threads), '-d', index_path, ref_genome_path]
    if verbose:
        print_log("Building rammap index: %s" % ' '.join(command))
    log = open('%s.log' % index_path, 'w'); subprocess.call(command, stderr=log); log.close()
    if verbose:
        print_log("rammap index built: %s" % index_path)

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

# helper function to start a mapper and return its stdout stream + log handle
def start_alignment(command, out_aln_path, aligner_name, bufsize=DEFAULT_BUFSIZE, verbose=True, log_path=None):
    if verbose:
        print_log("Aligning using %s: %s" % (aligner_name, ' '.join(command)))
    if log_path is None:
        log_path = '%s.log' % out_aln_path
    log = open(log_path, 'w', buffering=bufsize)
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, universal_newlines=True, bufsize=bufsize)
    except Exception as e:
        log.close(); raise e
    return process, log

# helper function to finish mapping (wait for streamed mapper and report errors)
def finish_alignment(process, log, aligner_name, out_aln_path, verbose=True):
    returncode = process.wait()
    if (process.stdout is not None) and (not process.stdout.closed):
        process.stdout.close()
    if not log.closed:
        log.close()
    if getattr(process, 'flag_viralmsa_delete_log', False) and isfile('Log.out'):
        remove('Log.out')
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, process.args)
    if verbose:
        print_log("%s alignment complete: %s" % (aligner_name, out_aln_path))

# helper function to tee a stream (yield every line and also copy it to sink)
def tee_stream(stream, sink):
    try:
        for line in stream:
            sink.write(line); yield line
    finally:
        if hasattr(stream, 'close'):
            stream.close()
        if not sink.closed:
            sink.close()

# align genomes using bowtie2
def align_bowtie2(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    command = ['bowtie2', '--very-sensitive', '-p', str(threads), '-f', '-x', '%s.bowtie2' % ref_genome_path, '-U', seqs_path]
    return start_alignment(command, out_aln_path, 'bowtie2', bufsize=bufsize, verbose=verbose)

# align genomes using BWA MEM
def align_bwa(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    command = ['bwa', 'mem', '-t', str(threads), ref_genome_path, seqs_path]
    return start_alignment(command, out_aln_path, 'BWA', bufsize=bufsize, verbose=verbose)

# align genomes using DRAGMAP
def align_dragmap(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    tmp_fq_path = '%s.fastq.gz' % out_aln_path
    fasta2fastq(seqs_path, tmp_fq_path)
    command = ['dragen-os', '--num-threads', str(threads), '--Aligner.sw-all', '1', '-r', '%s.DRAGMAP' % ref_genome_path, '-1', tmp_fq_path]
    return start_alignment(command, out_aln_path, 'DRAGMAP', bufsize=bufsize, verbose=verbose)

# align genomes using HISAT2
def align_hisat2(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    command = ['hisat2', '--very-sensitive', '-p', str(threads), '-f', '-x', '%s.hisat2' % ref_genome_path, '-U', seqs_path]
    return start_alignment(command, out_aln_path, 'HISAT2', bufsize=bufsize, verbose=verbose)

# align genomes using LRA
def align_lra(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    lra_ref_genome_path = '%s.lra' % ref_genome_path
    command = ['lra', 'align', '-t', str(threads), '-CONTIG', '-p', 's', lra_ref_genome_path, seqs_path]
    return start_alignment(command, out_aln_path, 'LRA', bufsize=bufsize, verbose=verbose)

# align genomes using minigraph
def align_minigraph(seqs_path, out_paf_path, ref_genome_path, threads, verbose=True):
    command = ['minigraph', '-c', '-t', str(threads), '--secondary=no', '-l', '0', '-d', '0', '-L', '0', ref_genome_path, seqs_path]
    return start_alignment(command, out_paf_path, 'Minigraph', verbose=verbose)

# align genomes using minimap2
def align_minimap2(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    index_path = '%s.minimap2.mmi' % ref_genome_path
    command = ['minimap2', '-t', str(threads), '--score-N=0', '--secondary=no', '--sam-hit-only', '-a', index_path, seqs_path]
    return start_alignment(command, out_aln_path, 'Minimap2', bufsize=bufsize, verbose=verbose)

# align genomes using mm2-fast
def align_mm2fast(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    index_path = '%s.mm2fast.mmi' % ref_genome_path
    command = ['mm2-fast', '-t', str(threads), '--score-N=0', '--secondary=no', '--sam-hit-only', '-a', index_path, seqs_path]
    return start_alignment(command, out_aln_path, 'mm2-fast', bufsize=bufsize, verbose=verbose)

# align genomes using NGMLR
def align_ngmlr(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    command = ['ngmlr', '--skip-write', '-x', 'pacbio', '-i', '0', '--no-smallinv', '-t', str(threads), '-r', ref_genome_path, '-q', seqs_path, '-o', '/dev/stdout']
    return start_alignment(command, out_aln_path, 'NGMLR', bufsize=bufsize, verbose=verbose)

# align genomes using rammap
def align_rammap(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    index_path = '%s.rammap.mmi' % ref_genome_path
    command = ['rammap', '-t', str(threads), '--score-N=0', '--secondary=no', '--sam-hit-only', '-a', index_path, seqs_path]
    return start_alignment(command, out_aln_path, 'rammap', bufsize=bufsize, verbose=verbose)

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
    log = open('%s.log' % out_aln_path, 'w', buffering=bufsize)
    ref_header, ref_seq = list(iter_fasta(ref_genome_path))[0]

    # perform alignment
    pool = Pool(processes=threads)
    complete = False
    try:
        inputs = ((command, ref_seq, curr_header, curr_seq) for curr_header, curr_seq in iter_fasta(seqs_path))
        for curr in pool.imap_unordered(run_needleman_wunsch, inputs, chunksize=1):
            log.write("Finished %s\n" % curr[0]); yield '>%s\n' % curr[0]; yield '%s\n' % curr[1]
        complete = True
    finally:
        if complete:
            pool.close()
        else:
            pool.terminate()
        pool.join(); log.close()
    if verbose:
        print_log("seq-align complete: %s" % out_aln_path)

# align genomes using STAR
def align_star(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    delete_log = not isfile('Log.out')
    index_path = '%s.STAR' % ref_genome_path
    out_sam_dir, out_sam_fn = split(out_aln_path)
    out_file_prefix = '%s.' % '.'.join(out_aln_path.split('.')[:-1])
    command = ['STAR', '--runThreadN', str(threads), '--genomeDir', index_path, '--readFilesIn', seqs_path, '--outFileNamePrefix', out_file_prefix, '--outFilterMismatchNmax', '9999999999', '--outStd', 'SAM']
    process, log = start_alignment(command, out_aln_path, 'STAR', bufsize=bufsize, verbose=verbose, log_path='%s/STAR.log' % out_sam_dir)
    process.flag_viralmsa_delete_log = delete_log
    return process, log

# align genomes using unimap
def align_unimap(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    index_path = '%s.umi' % ref_genome_path
    command = ['unimap', '-t', str(threads), '--score-N=0', '--secondary=no', '--sam-hit-only', '-a', '--cs', index_path, seqs_path]
    return start_alignment(command, out_aln_path, 'Unimap', bufsize=bufsize, verbose=verbose)

# align genomes using wfmash
def align_wfmash(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    command = ['wfmash', '-N', '--sam-format', '--threads=%d' % threads, ref_genome_path, seqs_path]
    return start_alignment(command, out_aln_path, 'wfmash', bufsize=bufsize, verbose=verbose)

# align genomes using Winnowmap
def align_winnowmap(seqs_path, out_aln_path, ref_genome_path, threads, bufsize=DEFAULT_BUFSIZE, verbose=True):
    index_path = '%s.meryl.k%d.txt' % (ref_genome_path, WINNOWMAP_K)
    command = ['winnowmap', '-k', str(WINNOWMAP_K), '-W', index_path, '-t', str(threads), '--score-N=0', '--secondary=no', '--sam-hit-only', '-a', ref_genome_path, seqs_path]
    return start_alignment(command, out_aln_path, 'Winnowmap', bufsize=bufsize, verbose=verbose)

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

    'rammap': {
        'check':       check_rammap,
        'build_index': build_index_rammap,
        'align':       align_rammap,
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

# parse user args
def parse_args():
    # check if user wants to update ViralMSA
    if '-u' in sys.argv or '--update' in sys.argv:
        update_viralmsa()

    # check if user just wants to list references
    if '-l' in sys.argv or '--list_references' in sys.argv:
        get_refs_json()
        if REFS_JSON is None:
            print("ERROR: Listing references requires an internet connection", file=sys.stderr); exit(1)
        print("=== Reference Sequence Aliases (for use with -s) ===")
        for ID in sorted(REFS_JSON.keys(), key=lambda x: REFS_JSON[x]['name']):
            print('* %s' % REFS_JSON[ID]['name'])
            print("  - Aliases: %s" % ', '.join(REFS_JSON[ID]['shortname']))
            print()
        exit(0)

    # use argparse to parse user arguments
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--sequences', required=True, type=str, help="Input Sequences (FASTA format, or SAM if already mapped)")
    parser.add_argument('-r', '--reference', required=True, type=str, help="Reference")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
    parser.add_argument('-e', '--email', required=False, type=str, default=None, help="Email Address (for Entrez)")
    parser.add_argument('-a', '--aligner', required=False, type=str, default=DEFAULT_ALIGNER, help="Aligner (options: %s)" % ', '.join(sorted(ALIGNERS.keys())))
    parser.add_argument('-t', '--threads', required=False, type=int, default=DEFAULT_THREADS, help="Number of Threads")
    parser.add_argument('-b', '--buffer_size', required=False, type=int, default=DEFAULT_BUFSIZE, help="File Stream Buffer Size (bytes)")
    parser.add_argument('-l', '--list_references', action='store_true', help="List all reference sequences")
    parser.add_argument('-or', '--omit_ref', action='store_true', help="Omit reference sequence from output alignment")
    parser.add_argument('-km', '--keep_mapping', action='store_true', help="Keep the mapper output in the output directory")
    parser.add_argument('-q', '--quiet', action='store_true', help="Suppress log output")
    parser.add_argument('-vd', '--viralmsa_dir', required=False, type=str, default=abspath(expanduser("~/.viralmsa")), help="ViralMSA Cache Directory")
    parser.add_argument('-u', '--update', action='store_true', help="Update ViralMSA (current version: %s)" % VERSION)
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
        args.reference = args.reference.upper()
        args.ref_path = '%s/%s' % (args.viralmsa_dir, args.reference)
        args.ref_genome_path = '%s/%s.fas' % (args.ref_path, args.reference)
        if not isfile(args.ref_genome_path):
            if args.email is None:
                print("ERROR: Must provide email address to download reference genome", file=sys.stderr); exit(1)
            get_refs_json()
            if REFS is None:
                print("ERROR: Unable to download reference genome if offline", file=sys.stderr); exit(1)
            tmp = args.reference.lower().replace(' ','').replace('-','').replace('_','')
            if tmp in REFS:
                args.reference = REFS[tmp].upper()
                args.ref_path = '%s/%s' % (args.viralmsa_dir, args.reference)
                args.ref_genome_path = '%s/%s.fas' % (args.ref_path, args.reference)

    # user args are valid, so return
    return args

# download reference genome
def download_ref_genome(reference, ref_path, ref_genome_path, email, bufsize=DEFAULT_BUFSIZE):
    if email is None:
        print("ERROR: Must provide email address to download reference genome", file=sys.stderr); exit(1)
    get_refs_json()
    if REFS is None:
        print("ERROR: Unable to download reference genome if offline", file=sys.stderr); exit(1)
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
def aln_to_fasta(aln_stream, aln_type, out_msa_path, ref_genome_path, omit_ref=False, bufsize=DEFAULT_BUFSIZE):
    msa = open(out_msa_path, 'w', buffering=bufsize)
    ref_header, ref_seq = list(iter_fasta(ref_genome_path))[0]; ref_seq_len = len(ref_seq)
    if not omit_ref:
        msa.write('>%s\n%s\n' % (ref_header, ref_seq))
    num_output_IDs = 0
    if aln_type == 'SAM':
        for l in aln_stream:
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
        for query in iter_fasta(aln_stream, bufsize=bufsize):
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
    print_log("Output Directory: %s" % args.output)
    if args.email is not None:
        print_log("Email Address: %s" % args.email)
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

    # prepare for alignment
    mapper_process = None
    mapper_log = None
    samtools_process = None
    samtools_log = None
    mapping_sink = None
    input_aln_file = None
    aln_stream = None
    out_aln_path = None

    # perform alignment
    try:
        if INPUT_TYPE == 'FASTA':
            # build aligner index (if needed)
            ALIGNERS[args.aligner]['build_index'](args.ref_genome_path, args.threads)
            print_log()
            print_log("===== ALIGNMENT =====")
            if args.aligner in ALIGNERS_FAS:
                aln_type = 'FASTA'; out_aln_path = '%s/reference_mapping.fas' % args.output
            elif args.aligner in ALIGNERS_PAF:
                aln_type = 'PAF'; out_aln_path = '%s/reference_mapping.paf' % args.output
            else:
                aln_type = 'SAM'; out_aln_path = '%s/reference_mapping.bam' % args.output
            mapper_result = ALIGNERS[args.aligner]['align'](args.sequences, out_aln_path, args.ref_genome_path, args.threads, bufsize=args.buffer_size)
            if args.aligner == 'seq-align':
                aln_stream = mapper_result
            else:
                mapper_process, mapper_log = mapper_result; aln_stream = mapper_process.stdout
            if args.keep_mapping:
                if aln_type == 'SAM':
                    samtools_log = open('%s/samtools.log' % args.output, 'w', buffering=args.buffer_size)
                    try:
                        samtools_process = subprocess.Popen(['samtools', 'view', '-b', '-o', out_aln_path, '-'], stdin=subprocess.PIPE, stderr=samtools_log, universal_newlines=True, bufsize=args.buffer_size)
                    except FileNotFoundError as exc:
                        samtools_log.close(); samtools_log = None; raise RuntimeError("--keep_mapping with SAM output requires samtools in PATH") from exc
                    mapping_sink = samtools_process.stdin
                else:
                    mapping_sink = open(out_aln_path, 'w', buffering=args.buffer_size)
                aln_stream = tee_stream(aln_stream, mapping_sink)
        elif INPUT_TYPE == 'SAM':
            aln_type = 'SAM'; out_aln_path = args.sequences; input_aln_file = open(out_aln_path, 'r', buffering=args.buffer_size); aln_stream = input_aln_file
        print_log("Converting alignment to FASTA...")
        out_msa_path = '%s/multiple_sequence_alignment.aln' % args.output
        num_output_IDs = aln_to_fasta(aln_stream, aln_type, out_msa_path, args.ref_genome_path, omit_ref=args.omit_ref, bufsize=args.buffer_size)
        if mapper_process is not None:
            finish_alignment(mapper_process, mapper_log, args.aligner, out_aln_path, verbose=not args.quiet)
            mapper_process = None; mapper_log = None
        if samtools_process is not None:
            samtools_returncode = samtools_process.wait()
            if samtools_log is not None and not samtools_log.closed:
                samtools_log.close()
            samtools_log = None
            if samtools_returncode != 0:
                raise subprocess.CalledProcessError(samtools_returncode, samtools_process.args)
            samtools_process = None
    finally:
        if hasattr(aln_stream, 'close'):
            aln_stream.close()
        if mapping_sink is not None and not mapping_sink.closed:
            mapping_sink.close()
        if mapper_process is not None:
            try:
                if mapper_process.poll() is None:
                    mapper_process.terminate()
            except ProcessLookupError:
                pass
            mapper_process.wait()
            if mapper_process.stdout is not None and not mapper_process.stdout.closed:
                mapper_process.stdout.close()
            if mapper_log is not None and not mapper_log.closed:
                mapper_log.close()
        if samtools_process is not None:
            try:
                if samtools_process.poll() is None:
                    samtools_process.terminate()
            except ProcessLookupError:
                pass
            samtools_process.wait()
        if samtools_log is not None and not samtools_log.closed:
            samtools_log.close()
        if INPUT_TYPE == 'FASTA' and args.aligner == 'dragmap' and out_aln_path is not None:
            tmp_fq_path = '%s.fastq.gz' % out_aln_path
            if isfile(tmp_fq_path):
                remove(tmp_fq_path)
        if input_aln_file is not None:
            input_aln_file.close()
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
