# 🧬 RSV-A Mutation Prediction & Protein Stability Analysis

> Deep learning-based nucleotide mutation forecasting for Respiratory Syncytial Virus A (RSV-A) using LSTM, DNABERT-6 transformer embeddings, rule-based ΔΔG protein stability estimation, and full B-cell + T-cell epitope prediction across all 11 RSV genes.

<br>

## 📌 Overview

Respiratory Syncytial Virus A (RSV-A) mutates rapidly due to the lack of proofreading in RNA polymerases, making vaccine development and antiviral drug design continuously challenging. This project presents a unified computational pipeline that:

- Constructs a **non-redundant RSV-A genome dataset** from NCBI using alignment-free 6-mer cosine similarity clustering
- Trains a **two-layer LSTM deep learning model** to predict future nucleotide mutations from genomic sequences
- Extracts **contextual sequence embeddings** using the pre-trained DNABERT-6 transformer model
- Performs **per-gene protein analysis** across all 11 RSV open reading frames (NS1, NS2, N, P, M, SH, G, F, M2-1, M2-2, L)
- Estimates **protein stability change (ΔΔG)** using an expanded 15-rule physico-chemical classification framework
- Predicts **B-cell and T-cell epitopes** (MHC-I and MHC-II) for all 11 RSV proteins
- Flags **mutations overlapping epitope regions** as high-impact immune escape candidates
- Delivers all results through an **interactive Gradio GUI** with CSV export

<br>

## 🏆 Key Results

| Metric | Value |
|---|---|
| Validation Accuracy | **97.49%** |
| Precision | **98.92%** |
| Recall | **98.92%** |
| F1-Score | **98.91%** |
| Sequences (raw → non-redundant) | 5,466 → 377 |
| Redundancy Reduction Rate | 0.931 |
| Shannon Diversity Index | 3.71 |
| Total ΔΔG (proxy) | 16.2 kcal/mol |
| Mutations in epitope regions | 19 / 41 (46.3%) |

<br>

## 🗂️ Project Structure

```
RSV_Project/
│
├── data/
│   ├── sequence.fasta          ← training sequences (download from NCBI)
│   ├── test_sequence.fasta     ← test RSV-A genome
│   └── wild_type.fasta         ← reference wild-type genome
│
├── data_preprocessing.py       ← length filter, N-filter, duplicate removal
├── redundancy_removal.py       ← 6-mer vectorisation + cosine clustering
├── feature_extraction.py       ← one-hot + position + GC encoding, sliding window
├── bert_feature_extraction.py  ← DNABERT-6 transformer embeddings (768-dim)
├── model.py                    ← LSTM architecture + training + plots
├── evaluation.py               ← confusion matrix, precision, recall, F1
├── mutation_prediction.py      ← sliding window inference + confidence scoring
├── protein_analysis.py         ← per-gene CDS extraction + AA comparison
├── ddg_prediction.py           ← 15-rule ΔΔG estimation + stability suggestions
├── epitope_prediction.py       ← B-cell (Parker, Kolaskar) + MHC-I/II epitopes
├── visualization.py            ← 5 publication-quality plots
├── main.py                     ← training pipeline (runs steps 1–6)
├── run_test.py                 ← full analysis pipeline (runs steps 7–14)
└── gradio_app.py               ← interactive 7-tab web GUI
```

<br>

## ⚙️ Pipeline — Step by Step

```
NCBI GenBank (5,466 sequences)
        ↓
  Data Preprocessing
  (length filter · N-filter · deduplication)
        ↓
  Redundancy Removal
  (6-mer vectors · cosine similarity · threshold 0.995)
  → 377 representative genomes
        ↓
  Feature Extraction
  (one-hot · relative position · GC content · window=50)
        ↓
  ┌─────────────────────┬──────────────────────┐
  │   LSTM Model        │   DNABERT-6           │
  │   (train / infer)   │   (gene embeddings)   │
  └─────────┬───────────┴──────────────────────┘
            ↓
  Mutation Prediction
  (confidence threshold 0.60 · gene-level annotation)
            ↓
  Per-Gene Protein Analysis
  (NS1 · NS2 · N · P · M · SH · G · F · M2-1 · M2-2 · L)
            ↓
  ΔΔG Stability Analysis
  (15 physico-chemical rules · stability suggestions)
            ↓
  Epitope Prediction
  (B-cell Parker + Kolaskar · MHC-I NetMHCpan · MHC-II NetMHCIIpan)
  → Mutation-epitope overlap detection
            ↓
  Visualisations + CSV Export + Gradio GUI
```

<br>

## 🚀 How to Run (Google Colab)

### Step 1 — Mount Google Drive and install dependencies

```python
from google.colab import drive
drive.mount('/content/drive', force_remount=True)

import os
os.chdir('/content/drive/MyDrive/RSV_Project')

!pip install biopython tensorflow scikit-learn seaborn pandas matplotlib
!pip install transformers torch sentencepiece requests gradio
```

### Step 2 — Write all module files

Run each `%%writefile` cell in order (CELL_02 through CELL_15) to save all `.py` files to your Drive folder.

### Step 3 — Add your FASTA files

Upload your sequences to `RSV_Project/data/`:
- `sequence.fasta` — training sequences from NCBI
- `test_sequence.fasta` — the RSV-A sequence you want to analyse
- `wild_type.fasta` — your reference wild-type genome

### Step 4 — Train the LSTM model

```bash
!python main.py
```

This trains the model and saves `model.keras` to your Drive. Outputs: `accuracy_plot.png`, `loss_plot.png`, `confusion_matrix.png`.

### Step 5 — Run the full analysis

```bash
!python run_test.py
```

This runs: mutation prediction → DNABERT embeddings → per-gene protein analysis → ΔΔG → epitope prediction → all visualisations → all CSV exports.

### Step 6 — Launch the Gradio web interface

```bash
!python gradio_app.py
```

Opens a public shareable link. Upload any RSV FASTA file and get all results in real time across 7 tabs.

<br>

## 🧠 Model Architecture

| Layer | Units | Dropout |
|---|---|---|
| Input | (50, 6) | — |
| LSTM Layer 1 | 256 | 0.4 |
| LSTM Layer 2 | 128 | 0.3 |
| Dense | 128 (ReLU) | — |
| Dense | 64 (ReLU) | — |
| Output (Softmax) | 4 | — |

- **Optimizer:** Adam
- **Loss:** Categorical cross-entropy
- **Batch size:** 128
- **Max epochs:** 25 (early stopping, patience=3)
- **Validation split:** 10%

<br>

## 🔬 New Features Added (Research-Level Extensions)

### 1. DNABERT-6 Transformer Embeddings
Pre-trained bidirectional transformer from HuggingFace (`zhihan1996/DNA_bert_6`). Each RSV gene sequence is tokenised into overlapping 6-mers and encoded into a 768-dimensional contextual [CLS] embedding. Complements the LSTM's local window with global bidirectional sequence understanding.

### 2. Per-Gene Protein Analysis (All 11 RSV Genes)
CDS coordinates from the reference genome (GenBank KX765698) are used to extract and translate each gene individually. Wild-type and mutated proteins are compared gene by gene, reporting amino acid changes and stop codon mutations per gene.

### 3. Expanded ΔΔG Rule-Based Framework (15 Rules)
Covers: hydrophobic core disruption, charge polarity inversions (positive ↔ negative), disulfide bond cysteine loss, proline helix-breaking, glycine-induced flexibility, aromatic stacking loss, and conservative substitutions. Each mutation receives a kcal/mol proxy value and an interpretation label.

### 4. Epitope Prediction Module
- **B-cell (Parker Hydrophilicity Scale)** — sliding window scoring, top 75th percentile
- **B-cell (Kolaskar-Tongaonkar Antigenicity Scale)** — propensity-based, top 75th percentile
- **MHC-I** — NetMHCpan via IEDB REST API (HLA-A\*02:01, 9-mer); rule-based anchor fallback
- **MHC-II** — NetMHCIIpan via IEDB REST API (DRB1\*01:01, 15-mer); rule-based amphipathicity fallback
- **Mutation-epitope overlap detection** — flags any predicted mutation inside an epitope region as HIGH IMPACT

### 5. Five New Visualisations
- Mutation frequency plot with RSV gene annotation strip
- Mutation confidence heatmap
- Per-gene mutation count bar chart
- Per-gene ΔΔG heatmap across all 11 genes
- Full epitope map across all 11 RSV proteins

<br>

## 📊 Output Files (saved to Google Drive)

| File | Description |
|---|---|
| `model.keras` | Trained LSTM model weights |
| `mutations.csv` | All predicted mutations with position, confidence, gene, severity |
| `ddg_results.csv` | ΔΔG values and interpretation for each amino acid change |
| `epitopes.csv` | All predicted B-cell and T-cell epitopes per gene |
| `stability_suggestions.csv` | Residue substitution recommendations for destabilising mutations |
| `epitope_disruptions.csv` | Mutations that overlap with predicted epitope regions |
| `accuracy_plot.png` | Training and validation accuracy curves |
| `loss_plot.png` | Training and validation loss curves |
| `confusion_matrix.png` | LSTM classification confusion matrix |
| `mutation_frequency.png` | Mutation frequency across genome positions |
| `mutation_heatmap.png` | Mutation confidence heatmap |
| `gene_mutation_chart.png` | Mutations per RSV gene |
| `ddg_gene_heatmap.png` | Cumulative ΔΔG per gene |
| `epitope_map.png` | Epitope regions across all 11 RSV proteins |

<br>

## 🖥️ Gradio GUI — 7 Tabs

| Tab | Contents |
|---|---|
| Overview | Model metrics, genome stats, variant classification, full summary report |
| Mutation Analysis | All predicted mutations table, top-5 high-confidence, CSV download |
| Per-Gene Proteins | AA changes and stop codon mutations for each of the 11 RSV genes |
| Protein Stability (ΔΔG) | Full ΔΔG table, stop codon table, stability suggestions, CSV download |
| BERT Embeddings | DNABERT-6 embedding summary for NS1 and NS2 |
| Epitope Analysis | All epitopes per gene, mutation-epitope overlap table, CSV download |
| Visualisations | All 5 analysis plots + training accuracy, loss, and confusion matrix |

<br>

## 📦 Dependencies

```
biopython
tensorflow
scikit-learn
pandas
numpy
matplotlib
seaborn
transformers
torch
sentencepiece
requests
gradio
```

Install all at once:
```bash
pip install biopython tensorflow scikit-learn pandas numpy matplotlib seaborn transformers torch sentencepiece requests gradio
```

<br>

## 📚 References

- Collins & Graham (2008) — RSV pathogenesis. *Journal of Virology*
- Hochreiter & Schmidhuber (1997) — LSTM. *Neural Computation*
- Ji et al. (2021) — DNABERT. *Bioinformatics*
- Sims et al. (2009) — k-mer alignment-free genome comparison. *PNAS*
- Kolaskar & Tongaonkar (1990) — B-cell epitope prediction. *FEBS Letters*
- Larsen et al. (2007) — NetMHCpan MHC-I prediction. *BMC Bioinformatics*
- Dehouck et al. (2009) — ΔΔG prediction. *Bioinformatics*

<br>

## 👩‍💻 Authors

**G. S. Akshaya · S. Bhumika · C. Sireesha**
Department of Artificial Intelligence and Engineering
Amrita School of Engineering, Bengaluru
Amrita Vishwa Vidyapeetam, India

<br>

## 📄 License

This project is developed for academic research purposes.

---

> *"From 9,211 raw nucleotide differences, the LSTM model identifies 41 biologically meaningful mutations — demonstrating context-aware deep learning prioritisation over naive sequence comparison."*
