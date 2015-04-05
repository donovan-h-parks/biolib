###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2014'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'

import os
import sys
import gzip
import traceback
import exceptions as exception


class InputFileException(exception.Exception):
    pass

    
protein_bases = {'a','r','n','d','c','q','e','g','h','i','l','k','m','f','p','s','t','w','y','v'}
nucleotide_bases {'a','c','g','t'}
insertion_bases = {'-', '.'}
    
def is_nucleotide(seq_file, req_perc=0.95, max_seqs_to_read=10):
    """Check if a file contains sequences in nucleotide space.
    
    The check is performed by looking for the characters in
    {a,c,g,t,n,.,-} and confirming that these comprise the
    majority of a sequences. A set number of sequences are
    read and the file assumed to be not be in nucleotide space 
    if none of these sequences are comprised primarily of the
    defined nucleotide set.

    Parameters
    ----------
    seq_file : str
        Name of fasta/q file to read.
    req_perc : float
        Percentage of bases in {a,c,g,t,n,.,-} before 
        declaring the sequences as being in nucleotide 
        space.
    max_seqs_to_read : int
        Maximum sequences to read before declaring 
        sequence file to not be in nucleotide space.

    Returns
    -------
    boolean
        True is sequences are in nucleotide space. 
    """
    
    for _seq_id, seq in read_seq(seq_file):
        seq = seq.lower()
        
        nt_bases = 0
        for c in (nucleotide_bases | {'n'} | insertion_bases):
            nt_bases = seq.count(c)
            
        if float(nt_bases) / len(seq) >= req_perc:
            return True
            
    return False

def is_protein(seq_file, req_perc=0.95, max_seqs_to_read=10):
    """Check if a file contains sequences in protein space.
    
    The check is performed by looking for the 20 amino acids,
    along with X, and the insertion characters '-' and '.', in
    order to confirm that these comprise the majority of a 
    sequences. A set number of sequences are read and the file 
    assumed to be not be in nucleotide space if none of these 
    sequences are comprised primarily of the defined nucleotide set.

    Parameters
    ----------
    seq_file : str
        Name of fasta/q file to read.
    req_perc : float
        Percentage of amino acid bases before 
        declaring the sequences as being in nucleotide 
        space.
    max_seqs_to_read : int
        Maximum sequences to read before declaring 
        sequence file to not be in nucleotide space.

    Returns
    -------
    boolean
        True is sequences are in protein space. 
    """
    
    for _seq_id, seq in read_seq(seq_file):
        seq = seq.lower()
        
        prot_bases = 0
        for c in (protein_bases | {'x'} | insertion_bases):
            prot_bases = seq.count(c)
            
        if float(prot_bases) / len(seq) >= req_perc:
            return True
            
    return False
    
def read(seq_file):
    """Read sequences from fasta/q file.

    Parameters
    ----------
    seq_file : str
        Name of fasta/q file to read.

    Returns
    -------
    dict : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    """
    if seq_file.endswith(('.fq.gz', '.fastq.gz', '.fq', '.fq.gz')):
        return read_fastq(seq_file)
    else:
        return read_fasta(seq_file)


def read_fasta(fasta_file):
    """Read sequences from fasta file.

    Parameters
    ----------
    fasta_file : str
        Name of fasta file to read.

    Returns
    -------
    dict : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    """

    if not os.path.exists(fasta_file):
        raise InputFileException

    if os.stat(fasta_file).st_size == 0:
        return

    try:
        open_file = open
        if fasta_file.endswith('.gz'):
            open_file = gzip.open

        seqs = {}
        for line in open_file(fasta_file):
            # skip blank lines
            if not line.strip():
                continue

            if line[0] == '>':
                seq_id = line[1:].split(None, 1)[0]
                seqs[seq_id] = []
            else:
                seqs[seq_id].append(line[0:-1])

        for seq_id, seq in seqs.iteritems():
            seqs[seq_id] = ''.join(seq)
    except:
        print traceback.format_exc()
        print ''
        print  "[Error] Failed to process sequence file: " + fasta_file
        sys.exit()

    return seqs


def read_fastq(fastq_file):
    """Read sequences from fastq file.

    Parameters
    ----------
    fastq_file : str
        Name of fastq file to read.

    Returns
    -------
    dict : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    """

    if not os.path.exists(fastq_file):
        raise InputFileException

    if os.stat(fastq_file).st_size == 0:
        return

    try:
        open_file = open
        if fastq_file.endswith('.gz'):
            open_file = gzip.open

        seqs = {}
        line_num = 0
        for line in open_file(fastq_file):
            line_num += 1

            if line_num == 1:
                seq_id = line[1:].split(None, 1)[0]
            elif line_num == 2:
                seqs[seq_id].seq = line[0:-1]
            elif line_num == 4:
                line_num = 0
    except:
        print traceback.format_exc()
        print ''
        print  "[Error] Failed to process sequence file: " + fastq_file
        sys.exit()

    return seqs


def read_seq(seq_file):
    """Generator function to read sequences from fasta/q file.

    This function is intended to be used as a generator
    in order to avoid having to have large sequence files
    in memory. Input file may be gzipped and in either
    fasta or fastq format. It is slightly more efficient
    to directly call read_fasta_seq() or read_fastq_seq()
    if the type of input file in known.

    Example:
    seq_io = SeqIO()
    for seq_id, seq in seq_io.read_seq(fasta_file):
        print seq_id
        print seq

    Parameters
    ----------
    seq_file : str
        Name of fasta/q file to read.

    Yields
    ------
    list : [seq_id, seq]
        Unique id of the sequence followed by the sequence itself.
    """

    if seq_file.endswith(('.fq.gz', '.fastq.gz', '.fq', '.fq.gz')):
        for rtn in read_fastq_seq(seq_file):
            yield rtn
    else:
        for rtn in read_fasta_seq(seq_file):
            yield rtn


def read_fasta_seq(fasta_file):
    """Generator function to read sequences from fasta file.

    This function is intended to be used as a generator
    in order to avoid having to have large sequence files
    in memory. Input file may be gzipped.

    Example:
    seq_io = SeqIO()
    for seq_id, seq in seq_io.read_fasta_seq(fasta_file):
        print seq_id
        print seq

    Parameters
    ----------
    fasta_file : str
        Name of fasta file to read.

    Yields
    ------
    list : [seq_id, seq]
        Unique id of the sequence followed by the sequence itself.
    """

    if not os.path.exists(fasta_file):
        raise InputFileException

    if os.stat(fasta_file).st_size == 0:
        return

    try:
        open_file = open
        if fasta_file.endswith('.gz'):
            open_file = gzip.open

        seq_id = None
        seq = None
        for line in open_file(fasta_file):
            # skip blank lines
            if not line.strip():
                continue

            if line[0] == '>':
                if seq_id != None:
                    yield seq_id, ''.join(seq)

                seq_id = line[1:].split(None, 1)[0]
                seq = []
            else:
                seq.append(line[0:-1])

        # report last sequence
        yield seq_id, ''.join(seq)
    except:
        print traceback.format_exc()
        print ''
        print  "[Error] Failed to process sequence file: " + fasta_file
        sys.exit()


def read_fastq_seq(fastq_file):
    """Generator function to read sequences from fastq file.

    This function is intended to be used as a generator
    in order to avoid having to have large sequence files
    in memory. Input file may be gzipped.

    Example:
    seq_io = SeqIO()
    for seq_id, seq in seq_io.read_fastq_seq(fastq_file):
        print seq_id
        print seq

    Parameters
    ----------
    fastq_file : str
        Name of fastq file to read.

    Yields
    ------
    list : [seq_id, seq]
        Unique id of the sequence followed by the sequence itself.
    """

    if not os.path.exists(fastq_file):
        raise InputFileException

    if os.stat(fastq_file).st_size == 0:
        return

    try:
        open_file = open
        if fastq_file.endswith('.gz'):
            open_file = gzip.open

        line_num = 0
        for line in open_file(fastq_file):
            line_num += 1

            if line_num == 1:
                seq_id = line[1:].split(None, 1)[0]
            elif line_num == 2:
                yield seq_id, line[0:-1]
            elif line_num == 4:
                line_num = 0
    except:
        print traceback.format_exc()
        print ''
        print  "[Error] Failed to process sequence file: " + fastq_file
        sys.exit()


def extract_seqs(fasta_file, seqs_to_extract):
    """Extract specific sequences from fasta file.

    Parameters
    ----------
    fasta_file : str
        Fasta file containing sequences.
    seqs_to_extract : set
        Ids of sequences to extract.

    Returns
    -------
    dict : dict[seq_id] -> seq
        Dictionary of sequences indexed by sequence id.
    """

    if not os.path.exists(fasta_file):
        raise InputFileException

    if os.stat(fasta_file).st_size == 0:
        return

    seqs = {}
    for line in open(fasta_file):
        if line[0] == '>':
            seq_id = line[1:].partition(' ')[0]

            seq_of_interest = False
            if seq_id in seqs_to_extract:
                seqs[seq_id] = []
                seq_of_interest = True
        elif seq_of_interest:
            seqs[seq_id].append(line[0:-1])

    for seq_id, seq in seqs.iteritems():
        seqs[seq_id] = ''.join(seq)

    return seqs


def seq_lengths(fasta_file):
    """Calculate length of each sequence.

    Parameters
    ----------
    fasta_file : str
        Fasta file containing sequences.

    Returns
    -------
    dict : d[seq_id] -> length
        Length of each sequence.
    """

    if not os.path.exists(fasta_file):
        raise InputFileException

    if os.stat(fasta_file).st_size == 0:
        return

    lens = {}
    for seq_id, seq in read_fasta_seq(fasta_file):
        lens[seq_id] = len(seq)

    return lens


def write_fasta(seqs, output_file):
    """Write sequences to fasta file.

    If the output file has the extension 'gz',
    it will be compressed using gzip.

    Parameters
    ----------
    seqs : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    output_file : str
        Name of fasta file to produce.
    """

    if output_file.endswith('.gz'):
        fout = gzip.open(output_file, 'wb')
    else:
        fout = open(output_file, 'w')

    for seq_id, seq in seqs.iteritems():
        fout.write('>' + seq_id + '\n')
        fout.write(seq + '\n')
    fout.close()
