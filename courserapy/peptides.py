from dataimport import masstable


def fragment(peptide, cyclic=False):
    fragments = []
    if cyclic == True:
        for string in [peptide[i:] + peptide[:i] for i in range(len(peptide))]:
            [fragments.append(string[0:j]) for j in range(1, len(peptide))]
        fragments.append(peptide)
    else:
        for i in range(len(peptide)):
            [fragments.append(peptide[i:i + j])
             for j in range(1, len(peptide) - i + 1)]
    return fragments


def theoreticalspectrum(peptide, cyclic=False):
    theory = [0]
    fragments = fragment(peptide, cyclic)
    mass = lambda frag: sum([acidsmasses[ch] for ch in frag])
    [theory.append(mass(frag)) for frag in fragments]
    return sorted(theory)


def encoding_peptides(dna, peptide):
    rna = transcribe(dna)
    k = len(peptide) * 3
    rng = range(0, (len(dna) - k + (1 * len(peptide))) * 2)
    return [reversetranscribe(rstrand) for rstrand in ["".join(rna)[i:i + k] for i in rng] if (translate(rstrand) == peptide or translate(reversecomplementrna(rstrand)) == peptide)]


def peptide_sequencing(spectrum, cyclic=False):
    aminoacids = [massesacids[m] for m in spectrum]
    length = len(aminoacids)
    peptides = itertools.permutations(aminoacids)
    return [[acidsmasses[pep] for pep in peptide] for peptide in peptides if theoreticalspectrum(peptide, cyclic) == spectrum]


def spectrum_score(peptide, spectrum, cyclic=False):
    theoretical = theoreticalspectrum(peptide, cyclic)
    if spectrum[-1] < theoretical[-1]:
        return -1
    else:
        return sum([min(theoretical.count(protein), spectrum.count(protein)) for protein in set(spectrum)])


'''def leaderboard_sequencing(spectrum, N, cyclic=False): need to translate this to just using masses
    lookup = masstable()
    acids_masses = lookup.masses
    peptides = dict()
    best_seqs = list([0, set()])
    for a in acids_masses:
        score = spectrum_score(a, spectrum, cyclic)
        if score in peptides:
            peptides[score].add(a)
        else:
            peptides[score] = {a}
    while len(peptides) > 0:
        growths = dict()
        for vals in peptides.values():
            for seq in vals:
                for a in aa:
                    growth = seq + a
                    if mass(growth) in spectrum:
                        score = spectrum_score(growth, spectrum, cyclic)
                        if score != -1:
                            if score in growths:
                                growths[score].add(growth)
                            else:
                                growths[score] = {growth}

        if len(growths.keys()) > 0:
            peptides = growths
            max_score = max(peptides.keys())
        else:
            break

        if max_score >= best_seqs[0]:
            best_seqs[0] = max_score
            best_seqs[1] = peptides[max_score]
        total = sum([len(l) for l in peptides.values()])
        if sum([len(l) for l in peptides.values()]) > N:
            min_score = min(peptides.keys())
            while total > N:
                if len(peptides[min_score]) != 0:
                    peptides[min_score].pop()
                    total -= 1
                else:
                    del peptides[min_score]
                    min_score = min(peptides.keys())
    return '-'.join([str(acidsmasses[ch]) for ch in best_seqs[1].pop()])


def spectral_convolution_with_multiplicities(spectrum):
    convolution = dict()
    for i, peak in enumerate(spectrum):
        for mass in spectrum[:i]:
            diff = peak - mass
            if diff >= 57 and diff <= 200:
                if diff in convolution:
                    convolution[diff] += 1
                else:
                    convolution[diff] = 1
    return convolution


def spectral_convolution(spectrum):
    return list(
        filter(lambda x: x != None,
               [peak - mass if peak - mass > 0 else None
                for mass in spectrum[:i]
                for i, peak in enumerate(spectrum)])
    )


def convolution_cyclopeptide_sequencing(M, N, spectrum):
    convolution = spectral_convolution_with_multiplicities(spectrum)
    multiplicities = list()
    total = 0
    temp = sorted(convolution.values(),reversed=True)
    for i, multiplicity in enumerate(temp):
        total += multiplicity
        if total > M:
            multiplicities = set(temp[:i])
            break
    leaderboard_spectrum = leaderboard_sequencing(
        [k for k in convolution if convolution[k] in multiplicities], N, True)'''


if __name__ == "__main__":
    #print(convolution_cyclopeptide_sequencing(20, 60,[int(x) for x in '57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493'.split(' ')]))
