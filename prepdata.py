import re
import collections
import math
import os

AA = 'ACDEFGHIKLMNPQRSTVWY'
diPeptides = [aa1 + aa2 for aa1 in AA for aa2 in AA]
AADict = {}
for i in range(len(AA)):
    AADict[AA[i]] = i


def AAC(sequence):
    count = collections.Counter(sequence)
    code = [count[aa] / len(sequence) for aa in AA]
    return code


group1 = {
    'hydrophobicity_PRAM900101': 'RKEDQN',
    'hydrophobicity_ARGP820101': 'QSTNGDE',
    'hydrophobicity_ZIMJ680101': 'QNGSWTDERA',
    'hydrophobicity_PONP930101': 'KPDESNQT',
    'hydrophobicity_CASG920101': 'KDEQPSRNTG',
    'hydrophobicity_ENGD860101': 'RDKENQHYP',
    'hydrophobicity_FASG890101': 'KERSQD',
    'normwaalsvolume': 'GASTPDC',
    'polarity': 'LIFWCMVY',
    'polarizability': 'GASDT',
    'charge': 'KR',
    'secondarystruct': 'EALMQKRH',
    'solventaccess': 'ALFCGIVW'
}
group2 = {
    'hydrophobicity_PRAM900101': 'GASTPHY',
    'hydrophobicity_ARGP820101': 'RAHCKMV',
    'hydrophobicity_ZIMJ680101': 'HMCKV',
    'hydrophobicity_PONP930101': 'GRHA',
    'hydrophobicity_CASG920101': 'AHYMLV',
    'hydrophobicity_ENGD860101': 'SGTAW',
    'hydrophobicity_FASG890101': 'NTPG',
    'normwaalsvolume': 'NVEQIL',
    'polarity': 'PATGS',
    'polarizability': 'CPNVEQIL',
    'charge': 'ANCQGHILMFPSTWYV',
    'secondarystruct': 'VIYCWFT',
    'solventaccess': 'RKQEND'
}
group3 = {
    'hydrophobicity_PRAM900101': 'CLVIMFW',
    'hydrophobicity_ARGP820101': 'LYPFIW',
    'hydrophobicity_ZIMJ680101': 'LPFYI',
    'hydrophobicity_PONP930101': 'YMFWLCVI',
    'hydrophobicity_CASG920101': 'FIWC',
    'hydrophobicity_ENGD860101': 'CVLIMF',
    'hydrophobicity_FASG890101': 'AYHWVMFLIC',
    'normwaalsvolume': 'MHKFRYW',
    'polarity': 'HQRKNED',
    'polarizability': 'KMHFRYW',
    'charge': 'DE',
    'secondarystruct': 'GNPSD',
    'solventaccess': 'MSPTHY'
}
groups = [group1, group2, group3]
properties = (
    'hydrophobicity_PRAM900101', 'hydrophobicity_ARGP820101', 'hydrophobicity_ZIMJ680101', 'hydrophobicity_PONP930101',
    'hydrophobicity_CASG920101', 'hydrophobicity_ENGD860101', 'hydrophobicity_FASG890101', 'normwaalsvolume',
    'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')


def CTDC(sequence):
    def Count(seq1, seq2):
        s = 0
        for aa in seq1:
            s = s + seq2.count(aa)
        return s

    code = []
    for p in properties:
        c1 = Count(group1[p], sequence) / len(sequence)
        c2 = Count(group2[p], sequence) / len(sequence)
        c3 = 1 - c1 - c2
        code = code + [c1, c2, c3]
    return code


def CTDD(sequence):
    def Count(aaSet, sequence):
        number = 0
        for aa in sequence:
            if aa in aaSet:
                number = number + 1
        cutoffNums = [1, math.floor(0.25 * number), math.floor(0.50 * number), math.floor(0.75 * number), number]
        cutoffNums = [i if i >= 1 else 1 for i in cutoffNums]

        code = []
        for cutoff in cutoffNums:
            myCount = 0
            for i in range(len(sequence)):
                if sequence[i] in aaSet:
                    myCount += 1
                    if myCount == cutoff:
                        code.append((i + 1) / len(sequence) * 100)
                        break
            if myCount == 0:
                code.append(0)
        return code

    code = []
    for p in properties:
        code = code + Count(group1[p], sequence) + Count(group2[p], sequence) + Count(group3[p], sequence)
    return code


def CTDT(sequence):
    code = []
    aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
    for p in properties:
        c1221, c1331, c2332 = 0, 0, 0
        for pair in aaPair:
            if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                c1221 = c1221 + 1
                continue
            if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                c1331 = c1331 + 1
                continue
            if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                c2332 = c2332 + 1
        code = code + [c1221 / len(aaPair), c1331 / len(aaPair), c2332 / len(aaPair)]
    return code


def DDE(sequence):
    myCodons = {'A': 4, 'C': 2, 'D': 2, 'E': 2, 'F': 2, 'G': 4, 'H': 2, 'I': 3, 'K': 2, 'L': 6, 'M': 1, 'N': 2, 'P': 4,
                'Q': 2, 'R': 6, 'S': 6, 'T': 4, 'V': 4, 'W': 1, 'Y': 2}
    myTM = []
    for pair in diPeptides:
        myTM.append((myCodons[pair[0]] / 61) * (myCodons[pair[1]] / 61))
    code = []
    tmpCode = [0] * 400
    for j in range(len(sequence) - 2 + 1):
        tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[
            sequence[j + 1]]] + 1
    if sum(tmpCode) != 0:
        tmpCode = [i / sum(tmpCode) for i in tmpCode]

    myTV = []
    for j in range(len(myTM)):
        myTV.append(myTM[j] * (1 - myTM[j]) / (len(sequence) - 1))

    for j in range(len(tmpCode)):
        tmpCode[j] = (tmpCode[j] - myTM[j]) / math.sqrt(myTV[j])

    code = code + tmpCode
    return code


def DPC(sequence):
    code = []
    tmpCode = [0] * 400
    for j in range(len(sequence) - 2 + 1):
        tmpCode[AADict[sequence[j]] * 20 + AADict[sequence[j + 1]]] = tmpCode[AADict[sequence[j]] * 20 + AADict[
            sequence[j + 1]]] + 1
    if sum(tmpCode) != 0:
        tmpCode = [i / sum(tmpCode) for i in tmpCode]
    code = code + tmpCode
    return code

def test():
    with open("threshold.txt") as f:
        sequences = [line.split()[0] for line in f.readlines()]
        with open('um/DPC.txt', 'w') as g:
            for i in sequences:
                g.write(i + ' ')
                code = DPC(i)
                for j in code: g.write(str(j) + ' ')
                g.write('\n')
        with open('um/AAC.txt', 'w') as g:
            for i in sequences:
                g.write(i + ' ')
                code = AAC(i)
                for j in code: g.write(str(j) + ' ')
                g.write('\n')
        with open('um/CTDC.txt', 'w') as g:
            for i in sequences:
                g.write(i + ' ')
                code = CTDC(i)
                for j in code: g.write(str(j) + ' ')
                g.write('\n')
        with open('um/CTDT.txt', 'w') as g:
            for i in sequences:
                g.write(i + ' ')
                code = CTDT(i)
                for j in code: g.write(str(j) + ' ')
                g.write('\n')
        with open('um/CTDD.txt', 'w') as g:
            for i in sequences:
                g.write(i + ' ')
                code = CTDD(i)
                for j in code: g.write(str(j) + ' ')
                g.write('\n')
        with open('um/DDE.txt', 'w') as g:
            for i in sequences:
                g.write(i + ' ')
                code = DDE(i)
                for j in code: g.write(str(j) + ' ')
                g.write('\n')
def work(sequence):
    code = AAC(sequence)
    code.extend(CTDC(sequence))
    code.extend(CTDD(sequence))
    code.extend(CTDT(sequence))
    code.extend(DDE(sequence))
    code.extend(DPC(sequence))
    return code


