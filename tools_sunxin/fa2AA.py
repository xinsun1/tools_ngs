__author__ = 'sunxin'

import sys

def nu2AA(nu) :
    codon_dict = {
    "GCT" : "A",
    "GCC" : "A",
    "GCA" : "A",
    "GCG" : "A",
    "CGT" : "R",
    "CGC" : "R",
    "CGA" : "R",
    "CGG" : "R",
    "AGA" : "R",
    "AGG" : "R",
    "AAT" : "N",
    "AAC" : "N",
    "GAT" : "D",
    "GAC" : "D",
    "TGT" : "C",
    "TGC" : "C",
    "CAA" : "Q",
    "CAG" : "Q",
    "GAA" : "E",
    "GAG" : "E",
    "GGT" : "G",
    "GGC" : "G",
    "GGA" : "G",
    "GGG" : "G",
    "CAT" : "H",
    "CAC" : "H",
    "ATT" : "I",
    "ATC" : "I",
    "ATA" : "I",
    "TTA" : "L",
    "TTG" : "L",
    "CTT" : "L",
    "CTC" : "L",
    "CTA" : "L",
    "CTG" : "L",
    "AAA" : "K",
    "AAG" : "K",
    "ATG" : "M",
    "TTT" : "F",
    "TTC" : "F",
    "CCT" : "P",
    "CCC" : "P",
    "CCA" : "P",
    "CCG" : "P",
    "TCT" : "S",
    "TCC" : "S",
    "TCA" : "S",
    "TCG" : "S",
    "AGT" : "S",
    "AGC" : "S",
    "ACT" : "T",
    "ACC" : "T",
    "ACA" : "T",
    "ACG" : "T",
    "TGG" : "W",
    "TAT" : "Y",
    "TAC" : "Y",
    "GTT" : "V",
    "GTC" : "V",
    "GTA" : "V",
    "GTG" : "V",
    "TAA" : "Stop",
    "TGA" : "Stop",
    "TAG" : "Stop"
    }

    if str(nu).upper() in codon_dict.keys() :
        return  codon_dict[nu]
    else :
        return nu


def main() :
    print('take in three codon')

    if len(sys.argv) == 1 :
        exit(0)

    print (sys.argv[1] + '\t' + nu2AA(sys.argv[1]))



if __name__ == "__main__" :
    main()