# -*- coding: utf-8 -*-
"""
FINDS GENES

@author: Sam Myers

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###
def separate_codons(dna):
    """Splits a DNA sequence string at every third letter.

    >>> separate_codons('ATGCATGAATGTAG')
    ['ATG', 'CAT', 'GAA', 'TGT', 'AG']
    """
    codons = [dna[i:i+3] for i in range(0, len(dna), 3)]
    return codons


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return complement[nucleotide]


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    reverse_dna = dna[::-1] # uses extended slice to reverse the dna string
    new_nucleotides = [get_complement(nucleotide) for nucleotide in reverse_dna]
    reverse_complement = ''.join(new_nucleotides)
    return reverse_complement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stop_codons = ['TAA', 'TAG', 'TGA']
    codons = separate_codons(dna)
    ORF = []
    for codon in codons:
        if codon in stop_codons:
            break
        else:
            ORF.append(codon)
    return ''.join(ORF)


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("TCGTAGATGGAATGTATAG")
    ['ATGGAATGTATAG']
    """
    stop_codons = ['TAA', 'TAG', 'TGA']
    codons = separate_codons(dna)
    all_ORFs = []
    ORF = []
    for codon in codons:
        if codon in stop_codons and ORF != []:
            all_ORFs.append(''.join(ORF))
            ORF = []
        elif codon == 'ATG' or ORF != []:
            ORF.append(codon)
    if ORF != []:
        all_ORFs.append(''.join(ORF))
    return all_ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    frame1, frame2, frame3 = dna, dna[1:], dna[2:]
    all_ORFs = find_all_ORFs_oneframe(frame1) + find_all_ORFs_oneframe(frame2) + find_all_ORFs_oneframe(frame3)
    return all_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    strand1, strand2 = dna, get_reverse_complement(dna)
    all_ORFs = find_all_ORFs(strand1) + find_all_ORFs(strand2)
    return all_ORFs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_ORFs = find_all_ORFs_both_strands(dna)
    longest = max(all_ORFs, key=len)
    
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specified DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    ORFs = []
    
    for i in range(num_trials):
        ORFs.append(longest_ORF(shuffle_string(dna)))
    
    maxlength = len(max(ORFs, key=len))
    return maxlength


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    amino_acids = []
    codons = separate_codons(dna)
    for codon in codons:
        try:
            amino_acids.append(aa_table[codon])
        except KeyError:
            continue
    return ''.join(amino_acids)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    all_ORFs = find_all_ORFs_both_strands(dna)
    AA_sequences = []
    for ORF in all_ORFs:
        if len(ORF) >= threshold:
            AA_sequences.append(coding_strand_to_AA(ORF))
    return AA_sequences


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print gene_finder(dna)
