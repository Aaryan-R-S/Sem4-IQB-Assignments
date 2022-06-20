# Mapping of Amino acids to their respective alpha-helix propensity values ----------------------------------------
p_helix = {
  "E": 1.53,  
  "A": 1.45,  
  "L": 1.34,  
  "H": 1.24,  
  "M": 1.20,  
  "Q": 1.17,  
  "W": 1.14,  
  "V": 1.14,  
  "F": 1.12,  
  "K": 1.07,  
  "I": 1.00,  
  "D": 0.98,  
  "T": 0.82,  
  "S": 0.79,  
  "R": 0.79,  
  "C": 0.77,  
  "N": 0.73,  
  "Y": 0.61,  
  "P": 0.59,  
  "G": 0.53  
}

# Mapping of Amino acids to their respective beta-strand propensity values ----------------------------------------
p_strand = {
  "M": 1.67,  
  "V": 1.65,  
  "I": 1.60,  
  "C": 1.30,  
  "Y": 1.29,  
  "F": 1.28,  
  "Q": 1.23,  
  "L": 1.22,  
  "T": 1.20,  
  "W": 1.19,  
  "A": 0.97,  
  "R": 0.90,  
  "G": 0.81,  
  "D": 0.80,  
  "K": 0.74,  
  "S": 0.72,  
  "H": 0.71,  
  "N": 0.65,  
  "P": 0.62,  
  "E": 0.26  
}

# Input Sequence --------------------------------------------------------------------------------------------------
# seq = "WHGCITVYWMTV"
seq = "SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDMLNPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF"
# print(seq)

# prefix_sum_p_helix(or strand)[i] stores the sum of first i amino acids' helix (or strand) propensity values appearing in the sequence
prefix_sum_p_helix = [0]
prefix_sum_p_strand = [0]

# num_former_p_helix(or strand)[i] stores the number of amino acids with helix (or strand) propensity values >= 1 from the first i amino acids appearing in the sequence
num_former_p_helix = [0]
num_former_p_strand = [0]

# Set the prefix_sum_p_helix(or strand) and num_former_p_helix(or strand) arrays ----------------------------------
for i in range(len(seq)):
    prefix_sum_p_helix.append(prefix_sum_p_helix[-1] + p_helix[seq[i]])
    prefix_sum_p_strand.append(prefix_sum_p_strand[-1] + p_strand[seq[i]])
    num_former_p_helix.append(num_former_p_helix[-1] + (1 if p_helix[seq[i]]>=1 else 0))
    num_former_p_strand.append(num_former_p_strand[-1] + (1 if p_strand[seq[i]]>=1 else 0))
    
# print(prefix_sum_p_helix)
# print(prefix_sum_p_strand)
# print(num_former_p_helix)
# print(num_former_p_strand)

# Predicting helices ----------------------------------------------------------------------------------------------
possibilities_helix = []        # To store all possible nucleation sites
nucleation_sites_helix = []     # To store all nucleation sites

# For assignment of the secondary structures (helices)
is_helix = []
for i in range(len(seq)):
    is_helix.append(" ")

# Look for each window of size 6 and check if it has atleast 4 residues with P(H)>= 1 then its a nucleation site
# Then expand window to left and then to right until average P(H)<1 in each case
for i in range(len(seq)-6+1):
    num_former_helix = num_former_p_helix[i+6]-num_former_p_helix[i]
    pos = ""
    for j in range(i, i+6):
        pos += seq[j]
    possibilities_helix.append(pos)
    if num_former_helix >= 4:
        # Assign current window residues the secondary structures
        for j in range(i, i+6):
            is_helix[j] = "H"
        nucleation_sites_helix.append(pos) 
        # Extend to left and assign the secondary structures
        k = i-1
        while k>=0:
            sum_former_helix = prefix_sum_p_helix[k+4] - prefix_sum_p_helix[k]
            if sum_former_helix<4:
                break
            else:
                is_helix[k] = "H"
            k-=1
        # Extend to right and assign the secondary structures
        k = i+6
        while k<len(seq):
            sum_former_helix = prefix_sum_p_helix[k+1] - prefix_sum_p_helix[k-3]
            if sum_former_helix<4:
                break
            else:
                is_helix[k] = "H"
            k+=1
        
# print(possibilities_helix)
# print(nucleation_sites_helix)
# print(is_helix)

# Predicting strands ---------------------------------------------------------------------------------------------
possibilities_strand = []        # To store all possible nucleation sites
nucleation_sites_strand = []     # To store all nucleation sites

# For assignment of the secondary structure (strands)
is_strand = []
for i in range(len(seq)):
    is_strand.append(" ")

# Look for each window of size 5 and check if it has atleast 3 residues with P(S)>= 1 then its a nucleation site
# Then expand window to left and then to right until average P(S)<1 in each case
for i in range(len(seq)-5+1):
    num_former_strand = num_former_p_strand[i+5]-num_former_p_strand[i]
    pos = ""
    for j in range(i, i+5):
        pos += seq[j]
    possibilities_strand.append(pos)
    if num_former_strand >= 3:
        # Assign current window residues the secondary structures
        for j in range(i, i+5):
            is_strand[j] = "S"
        nucleation_sites_strand.append(pos) 
        # Extend to left and assign the secondary structures
        k = i-1
        while k>=0:
            sum_former_strand = prefix_sum_p_strand[k+4] - prefix_sum_p_strand[k]
            if sum_former_strand<4:
                break
            else:
                is_strand[k] = "S"
            k-=1
        # Extend to right and assign the secondary structures
        k = i+5
        while k<len(seq):
            sum_former_strand = prefix_sum_p_strand[k+1] - prefix_sum_p_strand[k-3]
            if sum_former_strand<4:
                break
            else:
                is_strand[k] = "S"
            k+=1
        
# print(possibilities_strand)
# print(nucleation_sites_strand)
# print(is_strand)

# Resolve conflicts and assign final secondary structures --------------------------------------------------------
final_ss = []       # For assignment of the final secondary structure
i = 0

while i < len(seq):
    # Keep blank if not assigned any structure
    if is_helix[i] == " " and is_strand[i] == " ":
        final_ss.append(" ")
        
    # Assign strand if not assigned helix but assigned a strand
    elif is_helix[i] == " ":
        final_ss.append("S")
        
    # Assign helix if not assigned strand but assigned a helix
    elif is_strand[i] == " ":
        final_ss.append("H")
        
    # If assigned both helix and strand then resolve conflict by checking which is greater P(S) or P(H) for the conflicting part of the sequence
    # If P(S)<P(H) assign every amino acid of that part with helices
    # Otherwise assign every amino acid of that part a strand secondary structure
    else:
        k = i+1
        while k<len(seq) and is_helix[k] == "H" and is_strand[k] == "S":
            k+=1
        p_sum_helix = prefix_sum_p_helix[k] - prefix_sum_p_helix[i]
        p_sum_strand = prefix_sum_p_strand[k] - prefix_sum_p_strand[i]
        to_set = " "
        if p_sum_helix > p_sum_strand :
            to_set = "H"
        else:
            to_set = "S"
        for j in range(k-i):
            final_ss.append(to_set)
        i = k-1
    i+=1

# for i in range(len(seq)):
#     print(seq[i]+" "+final_ss[i])

print()
print("Assignment of Secondary Structure")
print(seq)
print("".join(is_helix))
print("".join(is_strand))
print()
print("Final Answer")
print(seq)
print("".join(final_ss))
print()