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
seq1 = "CTGTATC"
seq2 = "GTTACTGT"

# -> SAMPLE OUTPUT-1
# 8
# ----CTGTTATC
# GTTACTGT----

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
            score_matrix[i-1][j-1] + match_matrix[i-1][j-1],
            score_matrix[i-1][j] + penalty_gap,
            score_matrix[i][j-1] + penalty_gap,
        )
        
# --------------------------- [Find Max ii and jj] -------------------------
ii = 0
jj = len(seq2)

score = score_matrix[ii][jj]

for i in range(0, len(seq1)+1):
    if (score_matrix[i][len(seq2)]>score):
        score = score_matrix[i][len(seq2)]
        ii = i
        jj = len(seq2)

for j in range(0, len(seq2)+1):
    if (score_matrix[len(seq1)][j]>score):
        score = score_matrix[len(seq1)][j]
        ii = len(seq1)
        jj = j

# --------------------------- [Print Alignment and Score] -------------------------
temp_i = ii
temp_j = jj

while not (ii==0 or jj==0):
    if (ii>0 and jj>0 and score_matrix[ii][jj] == score_matrix[ii-1][jj-1] + match_matrix[ii-1][jj-1]):
        aligned_seq1 = seq1[ii-1] + aligned_seq1
        aligned_seq2 = seq2[jj-1] + aligned_seq2
        ii -= 1
        jj -= 1
    
    elif (ii>0 and score_matrix[ii][jj] == score_matrix[ii-1][jj] + penalty_gap):
        aligned_seq1 = seq1[ii-1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        ii -= 1
    
    else:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[jj-1] + aligned_seq2
        jj -= 1

while(ii!=0):
    aligned_seq1 =  seq1[ii-1] + aligned_seq1
    aligned_seq2 = "-" + aligned_seq2 
    ii -= 1
    
while(jj!=0):
    aligned_seq1 = "-" + aligned_seq1 
    aligned_seq2 =  seq2[jj-1] + aligned_seq2    
    jj -= 1

ii = temp_i
jj = temp_j
if (jj == len(seq2)):
    while (ii!=len(seq1)+1):
        aligned_seq1 = aligned_seq1 + seq1[ii-1]
        aligned_seq2 = aligned_seq2 + "-"
        ii += 1
else:
    while (jj!=len(seq2)+1):
        aligned_seq1 = aligned_seq1 + "-"
        aligned_seq2 = aligned_seq2 + seq2[jj-1]
        jj += 1
        
print(score)
print(aligned_seq1)
print(aligned_seq2)