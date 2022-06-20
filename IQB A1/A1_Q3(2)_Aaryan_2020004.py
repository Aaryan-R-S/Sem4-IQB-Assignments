# --------------------------- [Module Imports] -------------------------
import numpy as np

# --------------------------- [Scoring Scheme] -------------------------
reward_match = 2
penalty_mismatch = -1
penalty_gap = -2

# --------------------------- [Input Sequences] -------------------------
seq1 = "ATCAGAGTA"
seq2 = "TTCAGTA"

# --------------------------- [Matrix Initialization] -------------------------
match_matrix = np.zeros((len(seq1), len(seq2)))
substitution_matrix = np.zeros((len(seq1)+1, len(seq2)+1))

# Initializing match_matrix - A 2D matrix to store scores in case of match/mismatch of different pairs of characters each from one of the sequence provided
for i in range(len(seq1)):
    for j in range(len(seq2)):
        if(seq1[i]==seq2[j]):
            match_matrix[i][j] = reward_match
        else:
            match_matrix[i][j] = penalty_mismatch
            
# Initializing substitution_matrix - A bidimensional array to store the scores of the best local alignment where substitution_matrix[i][j] = score for optimal local alignment considering only first i characters of sequence 1 and first j characters of sequence 2
for i in range(len(seq1)+1):
    substitution_matrix[i][0] = max(0, i*penalty_gap)

for j in range(len(seq2)+1):
    substitution_matrix[0][j] = max(0, j*penalty_gap)
    
# --------------------------- [Fill Substitution Matrix] -------------------------
for i in range(1, len(seq1)+1):
    for j in range(1, len(seq2)+1):
        # Main Algortithm for filling Substitution Matrix in case of Local Alignment
        substitution_matrix[i][j] = max(
            0,
            substitution_matrix[i-1][j-1] + match_matrix[i-1][j-1],
            substitution_matrix[i-1][j] + penalty_gap,
            substitution_matrix[i][j-1] + penalty_gap,
        )

# --------------------------- [Print Substitution Matrix] -------------------------
print()
print("--------------------------------------------------------------------------------------------------------------------")
print('{0:>60}'.format("Substitution Matrix"))
print("--------------------------------------------------------------------------------------------------------------------")
for i in range(0, len(seq1)+1):
    if(i==0):
        print(end="                ")
    else:
        print('{0: >5}'.format(seq1[i-1]), end="   ")
print()

for j in range(0, len(seq2)+1):
    if(j==0):
        print(end="        ")
    else:
        print('{0: >5}'.format(seq2[j-1]), end="   ")
    for i in range(0, len(seq1)+1):
        print('{0: >5}'.format(substitution_matrix[i][j]), end="   ")
    print()
        
# --------------------------- [Find Max Score] -------------------------
score = 0
for i in range(0, len(seq1)+1):
    for j in range(0, len(seq2)+1):
        if (substitution_matrix[i][j]>score):
           score = substitution_matrix[i][j]
           
# --------------------------- [Construct Optimal Alignments by Retracing] -------------------------
def get_optimal_alignments(ii, jj, aligned_seq1, aligned_seq2, score):
    # Base Case
    if(substitution_matrix[ii][jj]==0):
        return [[aligned_seq1, aligned_seq2, score]]
    
    match_mismatch_alignments = []
    char_gap_alignments = []
    gap_char_alignments = []
    
    # Case 1: Aligned (ii-1)th character of seq1 with (jj-1)th character of seq2
    if (ii>0 and jj>0 and substitution_matrix[ii][jj] == substitution_matrix[ii-1][jj-1] + match_matrix[ii-1][jj-1]):
        match_mismatch_alignments = get_optimal_alignments(ii-1, jj-1, seq1[ii-1]+aligned_seq1, seq2[jj-1]+aligned_seq2, score+match_matrix[ii-1][jj-1])

    # Case 2: Aligned (ii-1)th character of seq1 with a Gap
    if (ii>0 and substitution_matrix[ii][jj] == substitution_matrix[ii-1][jj] + penalty_gap):
        char_gap_alignments = get_optimal_alignments(ii-1, jj, seq1[ii-1]+aligned_seq1, "-"+aligned_seq2, score+penalty_gap)
  
    # Case 3: Aligned a Gap with (jj-1)th character of seq2
    if (jj>0 and substitution_matrix[ii][jj] == substitution_matrix[ii][jj-1] + penalty_gap):
        gap_char_alignments = get_optimal_alignments(ii, jj-1, "-"+aligned_seq1, seq1[jj-1]+aligned_seq2, score+penalty_gap)
    
    return match_mismatch_alignments+char_gap_alignments+gap_char_alignments

# --------------------------- [Print Optimal Alignments and Scores] -------------------------
print()
print("--------------------------------------------------------------------------------------------------------------------")
print('{0:>67}'.format("Optimal Alignments with Scores"))
print("--------------------------------------------------------------------------------------------------------------------")

optimal_alignments_list = []

for i in range(0, len(seq1)+1):
    for j in range(0, len(seq2)+1):
        if (substitution_matrix[i][j]==score):
            optimal_alignments_list = optimal_alignments_list + get_optimal_alignments(i, j, "", "", 0)

for i in range(len(optimal_alignments_list)):
    print()
    print(optimal_alignments_list[i][0])
    for k in range(len(optimal_alignments_list[i][0])):
        if(optimal_alignments_list[i][0][k]==optimal_alignments_list[i][1][k]):
            print("|", end="")
        else:
            print(" ", end="")
    print()
    print(optimal_alignments_list[i][1])
    print("-----------")
    print("Score: " + str(optimal_alignments_list[i][2]))
    print("-----------")

print()
print("------------------------------------------------ END OF PROGRAM ------------------------------------------------------")
print()