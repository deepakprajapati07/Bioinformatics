# Custom Multiple Sequences Alignment Program
# Progressive alignment which uses pairwise alignment for all possible combinations

# Requirements: Set of Nucleotides sequences

# Execute the program using the following command: python custom_msa.py

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    # Initialization
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    for i in range(m + 1):
        dp[i][0] = i * gap
    for j in range(n + 1):
        dp[0][j] = j * gap
    
    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_score = match if seq1[i-1] == seq2[j-1] else mismatch
            dp[i][j] = max(
                dp[i-1][j-1] + match_score,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )
    
    # Traceback
    aligned_seq1, aligned_seq2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + gap:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:  # j > 0 and dp[i][j] == dp[i][j-1] + gap
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    return dp[m][n], aligned_seq1, aligned_seq2

def propagate_gaps_to_msa(msa, aligned_seq, seq_index):
    for i, char in enumerate(aligned_seq):
        if char == "-":  # Gap detected in aligned sequence
            for j in range(len(msa)):
                if j != seq_index:
                    msa[j] = msa[j][:i] + "-" + msa[j][i:]
    msa[seq_index] = aligned_seq  # Update aligned sequence in MSA
    return msa

def msa(sequences, match=1, mismatch=-1, gap=-2):
    msa_result = [sequences[0]]
    for seq in sequences[1:]:
        _, aligned_seq1, aligned_seq2 = needleman_wunsch(msa_result[0], seq, match, mismatch, gap)
        msa_result = propagate_gaps_to_msa(msa_result, aligned_seq1, 0)
        msa_result.append(aligned_seq2)

    # Ensure all sequences are the same length by padding with gaps
    max_length = max(len(seq) for seq in msa_result)
    msa_result = [seq.ljust(max_length, "-") for seq in msa_result]

    # Compute final score
    final_score = 0
    for i in range(len(msa_result)):
        for j in range(i + 1, len(msa_result)):
            score, _, _ = needleman_wunsch(msa_result[i], msa_result[j], match, mismatch, gap)
            final_score += score

    return final_score, msa_result

def get_user_sequences():
    """
    Prompt the user to input sequences via the terminal.
    """
    print("\n-----Input Instructions-----")
    print("Enter your sequences one by one. Press Enter after each sequence.")
    print("When done, press Enter on an empty line to finish.\n")

    sequences = []
    while True:
        seq = input(f"Enter sequence {len(sequences) + 1} (or press Enter to finish): ").strip()
        if not seq:
            break
        sequences.append(seq)

    if len(sequences) < 2:
        print("At least two sequences are required for alignment.")
        return []

    print(f"\nYou entered {len(sequences)} sequences.")
    return sequences

# Main Program
if __name__ == "__main__":
    print("---------------------------------------------------------")
    print("\n*** Welcome to Custom Sequence Alignment Tool ***\n")
    print("---------------------------------------------------------")
    sequences = get_user_sequences()

    if sequences:
        final_score, aligned_sequences = msa(sequences)

        print("\nPairwise Alignments and Scores:")
        for i in range(len(aligned_sequences)):
            for j in range(i + 1, len(aligned_sequences)):
                score, aligned1, aligned2 = needleman_wunsch(aligned_sequences[i], aligned_sequences[j])
                print(f"{aligned_sequences[i]} vs {aligned_sequences[j]} => Score: {score}")

        print(f"\nFinal MSA Score: {final_score}")
        print("Aligned Sequences:")
        for seq in aligned_sequences:
            print(seq)
