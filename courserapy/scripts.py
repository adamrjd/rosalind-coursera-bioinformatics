from os import getcwd
# Super-useful scripts


def flatten(listoflists):
    return [sublist for lists in listoflists for sublist in lists]

# Manipulating graphs


def circularize(element):
    if len(element) == 1:
        return element
    else:
        return [element[i:] + element[:i] for i in range(0, len(element))]


def decircularize(peptides):
    peptides = sorted(peptides)
    for i in range(0, len(peptides) / len(peptides[0])):
        cyclopeptide = peptides[i]
        removable = circularize(cyclopeptide)
        [peptides.remove(pep) for pep in removable
         if pep != cyclopeptide and pep in peptides]
    return peptides


def assemble(paths):
    reads = []
    for path in paths:
        read = path[0] + "".join([kmer[len(kmer) - 1:] for kmer in path[1:]])
        reads.append(read)
    return reads

# DNA manipulation numerically


def skew(dna):
    s = [0]
    skew_dict = {'G': 1, 'C': -1, 'A': 0, 'T': 0}
    for i, let in enumerate(dna):
        s.append(s[i] + skew_dict[let])
    return s


def hammingdistance(p, q):
    if len(p) != len(q):
        return 0
    else:
        return sum([1 for i in range(len(p)) if p[i] != q[i]])


def numbertopattern(num, k):
    return "".join(
        [numnuc[int(num / 4**(k - i)) % 4] for i in range(1, k + 1)])


def patterntonumber(path):
    l = len(path)
    return sum([nucnum[path[i]] * 4**(l - i - 1) for i in range(0, l)])

# RNA scripts

base = 'ACGT'
pairs = {base[index]: base[len(base) - 1 - index]
         for index, nucleotide in enumerate(base)}
nucnum = {nucleotide: index for index, nucleotide in enumerate(base)}
numnuc = {nucnum[key]: key for key in nucnum}

dnarna = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'U', 'U': 'T'}

rbase = 'ACGU'
rpairs = {rbase[index]: rbase[len(rbase) - 1 - index]
          for index, rnucleotide in enumerate(rbase)}
with open(getcwd() + '/datfiles/' + 'codons.txt', 'r') as inputdata:
    data = [line.strip().split(" ") for line in inputdata.readlines()]
    codons = {item[0]: item[1] for item in data}


def reversecomplement(dna):
    return "".join([pairs[letter] for i, letter in enumerate(dna[::-1])])


def reversecomplementrna(rna):
    return "".join([rpairs[letter] for i, letter in enumerate(rna[::-1])])


def translate(rna):
    return "".join([codons[rna[i:i + 3]] for i in range(0, len(rna) - 2, 3)
                    if codons[rna[i:i + 3]] != '#'])


def transcribe(dna):
    return "".join([dnarna[let] for let in dna])


def reversetranscribe(rna):
    return "".join([dnarna[let] for let in rna])


def shared_kmers(k, dna1, dna2):
    kmers_dict = dict()
    shared_kmers_lst = list()
    enumerate_dna = lambda dna: enumerate(
        [dna[_:_ + k] for _ in range(len(dna) - k + 1)])
    for i, kmer in enumerate_dna(dna1):
        if kmer in kmers_dict:
            kmers_dict[kmer].append(i)
        else:
            kmers_dict[kmer] = list([i])
    for j, kmer in enumerate_dna(dna2):
        if kmer in kmers_dict:
            for x in kmers_dict[kmer]:
                shared_kmers_lst.append((x, j))
        elif reversecomplement(kmer) in kmers_dict:
            for x in kmers_dict[reversecomplement(kmer)]:
                shared_kmers_lst.append((x, j))
    return list(sorted(shared_kmers_lst, key=lambda l: l[1]))

if __name__ == '__main__':
    with open('text.txt', 'r') as f:
        length, str1, str2 = [l.strip() for l in f.readlines()]
    with open('output.txt', 'w') as f:
        for item in shared_kmers(int(length), str1, str2):
            f.write(str(item) + '\n')
