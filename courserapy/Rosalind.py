base_pairs = 'ACGT'
rbase_pairs = 'ACGU'


def count_bp(s):
    base_dict = {b: 0 for b in base_pairs}
    for bp in s:
        base_dict[bp] += 1
    return [base_dict[bp] for bp in base_pairs]


def transcribe(dna):
    rna_dict = {a: b for a, b in zip(base_pairs, rbase_pairs)}
    return ''.join([rna_dict[bp] for bp in dna])


def reverse_complement(dna):
    dna_dict = {a: b for a, b in zip(base_pairs, list(reversed(base_pairs)))}
    return ''.join([dna_dict[bp] for bp in dna[::-1]])

with open('input.txt', mode='r') as f:
    print(reverse_complement(f.readline().strip()))
