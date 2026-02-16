from Bio import SeqIO
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from collections import Counter
import math

input_file = "cleaned_sequences.fasta"   # Use preprocessed file
output_file = "nr_kmer_cosine.fasta"

k = 6
threshold = 0.995

# ---------------------------
# 1. Load Sequences
# ---------------------------
records = list(SeqIO.parse(input_file, "fasta"))

filtered = []
for record in records:
    seq = str(record.seq).upper()
    if len(seq) > 0:
        filtered.append(record)

print("Original sequences:", len(records))
print("Sequences used for clustering:", len(filtered))

# ---------------------------
# 2. Build k-mer Vocabulary
# ---------------------------
def get_kmers(seq, k):
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

vocab = set()
kmer_counts = []

for record in filtered:
    seq = str(record.seq)
    kmers = get_kmers(seq, k)
    counts = Counter(kmers)
    kmer_counts.append(counts)
    vocab.update(counts.keys())

vocab = list(vocab)

# ---------------------------
# 3. Convert to Frequency Matrix
# ---------------------------
matrix = []

for counts in kmer_counts:
    vec = [counts.get(kmer, 0) for kmer in vocab]
    matrix.append(vec)

matrix = np.array(matrix, dtype=float)

# Safe normalization
norms = np.linalg.norm(matrix, axis=1, keepdims=True)
norms[norms == 0] = 1
norm_matrix = matrix / norms

# ---------------------------
# 4. Cosine Similarity Matrix
# ---------------------------
sim_matrix = cosine_similarity(norm_matrix)

# ---------------------------
# 5. Greedy Clustering
# ---------------------------
n = len(filtered)
visited = set()
clusters = []

for i in range(n):
    if i in visited:
        continue
    
    cluster = [i]
    visited.add(i)
    
    for j in range(i + 1, n):
        if j not in visited and sim_matrix[i][j] >= threshold:
            cluster.append(j)
            visited.add(j)
    
    clusters.append(cluster)

print("Number of clusters:", len(clusters))

# ---------------------------
# 6. Save Representatives
# ---------------------------
representatives = [filtered[cluster[0]] for cluster in clusters]
SeqIO.write(representatives, output_file, "fasta")

# ---------------------------
# 7. Evaluation Metrics
# ---------------------------

# 1️⃣ Redundancy Reduction Rate
RRR = 1 - (len(representatives) / len(filtered))

# 2️⃣ Intra-cluster Similarity
intra_sims = []
for cluster in clusters:
    if len(cluster) > 1:
        sims = []
        for i in cluster:
            for j in cluster:
                if i != j:
                    sims.append(sim_matrix[i][j])
        intra_sims.append(np.mean(sims))

avg_intra = np.mean(intra_sims) if intra_sims else 0

# 3️⃣ Inter-cluster Similarity
rep_indices = [cluster[0] for cluster in clusters]
inter_sims = []

for i in range(len(rep_indices)):
    for j in range(i+1, len(rep_indices)):
        inter_sims.append(sim_matrix[rep_indices[i]][rep_indices[j]])

avg_inter = np.mean(inter_sims) if inter_sims else 0

# 4️⃣ Shannon Diversity Index
cluster_sizes = [len(cluster) for cluster in clusters]
total = sum(cluster_sizes)

H = 0
for size in cluster_sizes:
    p = size / total
    H -= p * math.log(p)

# ---------------------------
# Print Results
# ---------------------------
print("\n--- Evaluation Metrics ---")
print("Redundancy Reduction Rate:", round(RRR, 4))
print("Average Intra-cluster Similarity:", round(avg_intra, 4))
print("Average Inter-cluster Similarity:", round(avg_inter, 4))
print("Shannon Diversity Index:", round(H, 4))