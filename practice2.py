import functools
import time


# cache = {}

# @functools.lru_cache
# def EditDistance(s1, s2, n):
#     # def new_EditDistance(s1, s2, n):
#     n += 1
#     print("call number: ", n)
#     if not s1 or not s2:
#         return len(s1) + len(s2)
#     else:
#         if s1[0] == s2[0]:
#             return EditDistance(s1[1:], s2[1:], n)
#         else:
#             return min(1 + EditDistance(s1, s2[1:], n),
#                        1 + EditDistance(s1[1:], s2, n),
#                        1 + EditDistance(s1[1:], s2[1:], n))


# key = s1 + "#" + s2
# if key in cache:
#     return cache[key]
# else:
#     result = new_EditDistance(s1, s2, n)
#     cache[key] = result
#     return result
import numpy as np


def simpleMatch(a, b):
    return 1 if a == b else -1

def distanceMatch(a, b):
    return 0 if a == b else -1

def linearGap(n):
    return -1 * n

def alignmentScore(s1, s2, gapPenalty, match):
    if not s1 or not s2:
        return gapPenalty(len(s1)) + gapPenalty(len(s2))
    else:
        return max(gapPenalty(1) + alignmentScore(s1, s2[1:], gapPenalty, match),
                   gapPenalty(1) + alignmentScore(s1[1:], s2, gapPenalty, match),
                   match(s1[0], s2[0]) + alignmentScore(s1[1:], s2[1:], gapPenalty, match))


# print(alignmentScore("GAAACACCCGGAGCATATGCTG", "GCAAACATCCGAGCATATGCTG", linearGap, simpleMatch))



def alignmentScoreDP(s1, s2, gapPenalty, match):
    m = np.zeros((len(s1) + 1, len(s2) + 1)) # creates an n by m matrix of zeros, using numpy
    m[0, 0] = 0
    for i in range(1, len(s1) + 1):
        m[i, 0] = gapPenalty(i)
    for j in range(1, len(s2) + 1):
        m[0, j] = gapPenalty(j)
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            m[i, j] = max(gapPenalty(1) + m[i, j - 1],
                          gapPenalty(1) + m[i - 1, j],
                          match(s1[i - 1], s2[j - 1]) + m[i - 1, j - 1])
    return m
print(alignmentScoreDP("ATGCCC", "GACTGGG", linearGap, simpleMatch))

def readAlignmentG(s1, s2, m, gapPenalty, match):
    i = len(s1)
    j = len(s2)
    s1a = ""
    s2a = ""
    score = 0
    while i > 0 or j > 0:
        if i > 0 and j > 0 and m[i, j] == m[i - 1, j - 1] + match(s1[i - 1], s2[j - 1]):
            i = i - 1
            j = j - 1
            s1a = s1[i] + s1a
            s2a = (s2[j] if s1[i] == s2[j] else s2[j].lower()) + s2a
            score += match(s1[i], s2[j])
        else:
            foundit = False
            for g in range(1, i + 1):
                if m[i, j] == m[i - g, j] + gapPenalty(g):
                    s1a = s1[i - g:i] + s1a
                    s2a = ('-' * g) + s2a
                    i = i - g
                    score += gapPenalty(g)
                    foundit = True
                    break
            if not foundit:
                for g in range(1, j + 1):
                    if m[i, j] == m[i, j - g] + gapPenalty(g):
                        s1a = ('-' * g) + s1a
                        s2a = s2[j - g:j] + s2a
                        j = j - g
                        score += gapPenalty(g)
                        foundit = True
                        break
            assert foundit
    return s1a, s2a, score