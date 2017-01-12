from os import getcwd
CWD = getcwd()
DATDIR = '/datfiles/'


class masstable(object):
    '''Canned lists and methods for amino acid masses in peptide sequencing problems'''

    masses = list()
    acids = list()
    __mass_dict = dict()

    def __init__(self):
        self.parse(CWD + DATDIR + 'masstable.txt')

    def parse(self, file):

        with open(file, 'r') as f:
            data = [line.strip().split(' ') for line in f.readlines()]
        for item in data:
            self.acids.append(item[0])
            self.masses.append(int(item[1]))
            self.__mass_dict[item[0]] = int(item[1])

    def get_mass(self, acid):
        return self.__mass_dict[acid]


class Score(object):
    '''Canned scoring matrices for amino acid sequence alignment problems'''
    NEGATIVE_IDENTITY = dict()
    INVERSE_IDENTITY = dict()
    IDENTITY = dict()
    ALL_ONES = dict()
    BLOSUM62 = dict()
    PAM250 = dict()

    def __init__(self):
        # identity matrices for all amino acids
        self.NEGATIVE_IDENTITY = Score.dictify(
            CWD + DATDIR + 'NEGATIVE_IDENTITY.txt')
        self.INVERSE_IDENTITY = Score.dictify(
            CWD + DATDIR + 'INVERSE_IDENTITY.txt')
        self.IDENTITY = Score.dictify(CWD + DATDIR + 'IDENTITY.txt')
        self.ALL_ONES = Score.dictify(CWD + DATDIR + 'ALL_ONES.txt')
        # biologically significant matrices
        self.BLOSUM62 = Score.dictify(CWD + DATDIR + 'BLOSUM62.txt')
        self.PAM250 = Score.dictify(CWD + DATDIR + 'PAM250.txt')

    @staticmethod
    def isnumeric(s):
        ''' builtin str.isnumeric() does not account for '-' '''
        l = "-0123456789"
        return all([x in list(l) for x in list(s)])

    @staticmethod
    def dictify(file):
        '''Canned method to create dict of dicts...definitely better ways to do this'''
        with open(file, 'r') as input_data:
            data = [[int(x) if Score.isnumeric(x) else x
                     for x in line.strip().split(' ') if x != '']
                    for line in input_data.readlines()]
        data_dict = dict()
        matrix = [[[let, [data[j][0], data[i + 1][j]]]
                   for j in range(1, len(data))]
                  for i, let in enumerate(data[0])]
        for item in matrix:
            for k, kv in item:
                if k in data_dict:
                    l = [[key, data_dict[k][key]] for key in data_dict[k]]
                    l.append([kv[0], kv[1]])
                    d = {}
                    d[k] = {ky: vl for ky, vl in l}
                    data_dict.update(d)
                else:
                    data_dict[k] = {kv[0]: kv[1]}
        return data_dict

if __name__ == "__main__":
    test = Score().INVERSE_IDENTITY
    print(test)
