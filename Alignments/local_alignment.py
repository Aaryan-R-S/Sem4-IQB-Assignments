# --------------------------- [Module Imports] -------------------------
import numpy as np

# --------------------------- [Scoring Scheme] -------------------------
reward_match = 2
penalty_mismatch = -1
penalty_gap = -1

# --------------------------- [Input Sequences] -------------------------
# seq1 = input("Enter Sequence No.1: ")
# seq2 = input("Enter Sequence No.2: ")

aligned_seq1 = ""
aligned_seq2 = ""

# -> SAMPLE INPUT-1
seq1 = "TTGTATC"
seq2 = "ACCGGTAT"

# -> SAMPLE OUTPUT-1
# 8
# GTAT
# GTAT

# --------------------------- [Matrix Initialization] -------------------------
match_matrix = np.zeros((len(seq1), len(seq2)))
score_matrix = np.zeros((len(seq1)+1, len(seq2)+1))

# For match_matrix
for i in range(len(seq1)):
    for j in range(len(seq2)):
        if(seq1[i]==seq2[j]):
            match_matrix[i][j] = reward_match
        else:
            match_matrix[i][j] = penalty_mismatch
            
# For score_matrix
for i in range(len(seq1)+1):
    score_matrix[i][0] = max(0, i*penalty_gap)

for j in range(len(seq2)+1):
    score_matrix[0][j] = max(0, j*penalty_gap)
    
# --------------------------- [Fill Score Matrix] -------------------------
for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):
        score_matrix[i][j] = max(
            0,
            score_matrix[i-1][j-1] + match_matrix[i-1][j-1],
            score_matrix[i-1][j] + penalty_gap,
            score_matrix[i][j-1] + penalty_gap,
        )
        
# --------------------------- [Find Max ii and jj] -------------------------
ii = 0
jj = 0

score = score_matrix[ii][jj]

for i in range(0, len(seq1)+1):
    for j in range(0, len(seq2)+1):
        if (score_matrix[i][j]>score):
           score = score_matrix[i][j]
           ii = i
           jj = j

# --------------------------- [Print Alignment and Score] -------------------------
while (score_matrix[ii][jj]>0):
    aligned_seq1 = seq1[ii-1] + aligned_seq1
    aligned_seq2 = seq2[jj-1] + aligned_seq2
    ii -= 1
    jj -= 1

print(score)
print(aligned_seq1)
print(aligned_seq2)