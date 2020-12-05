import sys


def reverse_complement(s):
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    output = []
    for c in s:
        if c in 'ACGT':
            output.append(complement_map[c])
    return ''.join(output[::-1])


def three_frame_translation(s):
    frame_0 = [s[i:i+3] for i in range(0, len(s), 3)]
    frame_1 = [s[1:][i:i+3] for i in range(0, len(s[1:]), 3)]
    frame_1.insert(0, s[:1])
    frame_2 = [s[2:][i:i+3] for i in range(0, len(s[2:]), 3)]
    frame_2.insert(0, s[:2])
    return [frame_0, frame_1, frame_2]


def six_frame_translation(s):
    output = []
    codon_map = {'UUU': 'F',
                 'UUC': 'F',
                 'UUA': 'L',
                 'UUG': 'L',
                 'UCU': 'S',
                 'UCC': 'S',
                 'UCA': 'S',
                 'UCG': 'S',
                 'UAU': 'Y',
                 'UAC': 'Y',
                 'UAA': '*',
                 'UAG': '*',
                 'UGU': 'C',
                 'UGC': 'C',
                 'UGA': '*',
                 'UGG': 'W',
                 'CUU': 'L',
                 'CUC': 'L',
                 'CUA': 'L',
                 'CUG': 'L',
                 'CCU': 'P',
                 'CCC': 'P',
                 'CCA': 'P',
                 'CCG': 'P',
                 'CAU': 'H',
                 'CAC': 'H',
                 'CAA': 'Q',
                 'CAG': 'Q',
                 'CGU': 'R',
                 'CGC': 'R',
                 'CGA': 'R',
                 'CGG': 'R',
                 'AUU': 'I',
                 'AUC': 'I',
                 'AUA': 'I',
                 'AUG': 'M',
                 'ACU': 'T',
                 'ACC': 'T',
                 'ACA': 'T',
                 'ACG': 'T',
                 'AAU': 'N',
                 'AAC': 'N',
                 'AAA': 'K',
                 'AAG': 'K',
                 'AGU': 'S',
                 'AGC': 'S',
                 'AGA': 'R',
                 'AGG': 'R',
                 'GUU': 'V',
                 'GUC': 'V',
                 'GUA': 'V',
                 'GUG': 'V',
                 'GCU': 'A',
                 'GCC': 'A',
                 'GCA': 'A',
                 'GCG': 'A',
                 'GAU': 'D',
                 'GAC': 'D',
                 'GAA': 'E',
                 'GAG': 'E',
                 'GGU': 'G',
                 'GGC': 'G',
                 'GGA': 'G',
                 'GGG': 'G', }
    translations = three_frame_translation(
        s) + three_frame_translation(reverse_complement(s))
    for translation in translations:
        current_amino_acid_seq = []
        for codon in translation:
            rna_codon = codon.replace("T", "U")
            if rna_codon in codon_map:
                current_amino_acid_seq.append(codon_map[rna_codon])
        output.append(''.join(current_amino_acid_seq))
    return output
