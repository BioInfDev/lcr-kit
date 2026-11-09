# lcr-kit
Tools for detecting low-complexity regions (LCR) in protein sequences. Implements two independent methods: **Shannon Entropy** and **Dotplot-Convolve-FDR**.  

The project consists of two parts:
* core — a library with the algorithm implementations;
* cli — command-line tools `lcr-shannon` and `lcr-dotplot` that output results as TSV.

# Methods
## Shannon Entropy

The sequence is scanned with a sliding window of fixed size. For each window, Shannon entropy is computed from the symbol frequency composition. Windows whose entropy is below a given threshold are marked as low-complexity and merged into final regions.

The method is fast and predictable, well suited as a baseline and for screening large sequence sets.

## Dotplot-Convolve-FDR

A self-dotplot of the sequence is constructed, a convolution with a kernel is applied to the values along the diagonal, and the significance of the resulting peaks is assessed with FDR control.

The method is more robust to tandem repeats and regions with heterogeneous composition than a simple entropy threshold, at the cost of longer computation time.

# Installing
1. Cloning repository
```bash 
git clone https://github.com/BioInfDev/lcr-kit.git
```
2. Change directory to lcr-kit
```bash
cd lcr-kit
```
3. Create and activate venv
```bash
python3 -m venv venv
source venv/bin/activate
```  
4. Building
```bash
pip install .
```

# CLI Usage
Both tools read sequences in FASTA format and write TSV.

## lcr-shannon
Detects low-complexity regions using a sliding-window Shannon entropy threshold.

```bash
lcr-shannon <fasta> <protein> [--window N] [--min_length N]
```

**Arguments**

| Argument | Type | Description |
|----------|------|-------------|
| `fasta` | path | Path to a FASTA file with a protein sequence. |
| `protein` | string | Gene or protein name. |

**Options**

| Flag | Default | Description |
|------|---------|-------------|
| `-w, --window` | `3` | Sliding window size. |
| `-l, --min_length` | `4` | Minimum length of a reported LCR. |

**Example**

```bash
lcr-shannon P12345.fasta P12345 --window 3 --min_length 4
```

## lcr-dotplot
Detects low-complexity regions by convolving a self-dotplot along the diagonal and applying FDR control.

```bash
lcr-dotplot <fasta> <protein> [--kernel-size N] [--fdr-rate F] [--min_length N]
```

**Arguments**

| Argument | Type | Description |
|----------|------|-------------|
| `fasta` | path | Path to a FASTA file with a protein sequence. |
| `protein` | string | Gene or protein name. |

**Options**

| Flag | Default | Description |
|------|---------|-------------|
| `-k, --kernel-size` | `3` | Kernel size for the dotplot convolution. |
| `-t, --fdr-rate` | `0.01` | Target FDR rate. |
| `-l, --min_length` | `4` | Minimum length of a reported LCR. |

**Example**

```bash
lcr-dotplot P12345.fasta P12345 --kernel-size 3 --fdr-rate 0.01 --min_length 4
```

## Output TSV
File structure:
1. Title cell:
   ```Low-complexity regions of <protein> | Method: 'Dotplot-Convolve-FDR' (kernel = <kernel_size>x<kernel_size>, fdr-rate = <fdr_rate>) | Minimum LCR length = <min_length>```
2. Header row:
   ```Sequence	Start	End```
3. LCR per row.

### Example file (`lcr_shannon_P12345.tsv`):

```tsv
Low-complexity regions of P12345 | Method: 'Shannon Entropy' (window = 3) | Minimum LCR length = 4
Sequence	Start	End
AAAAAAAAAA  45  78
GGGGGGGGGGG  112  140
```

### Example file (`lcr_dotplot_P12345.tsv`):

```tsv
Low-complexity regions of P12345 | Method: 'Dotplot-Convolve-FDR' (kernel = 3x3, fdr-rate = 0.01) | Minimum LCR length = 4
Sequence	Start	End
AAAAAAAAAA	45	78
GGGGGGGGGGG	112	140
```

# Library Usage
```Python
from lcr_kit.core.detector import ShannonDetector, DotplotDetector


seq = "MSETTAPPPQQQPPPP"
min_length = 4

window = 3
shannon_detector = ShannonDetector(sequence=seq, window=window, min_length=min_length)
shannon_detector.run()

kernel_size, fdr_rate = 3, 0.01
dotplot_detector = DotplotDetector(sequence=seq, kernel_size=kernel_size, fdr_rate=fdr_rate, min_length=min_length)
dotplot_detector.run()
```
