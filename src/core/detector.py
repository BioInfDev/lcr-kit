import random
from itertools import groupby
from typing import NamedTuple, List, Optional

import numpy as np
from scipy.stats import entropy
from scipy.signal import convolve2d

import src.core.seq2smth as seq2smth

class LCR(NamedTuple):
    seq: str
    start: int
    end: int

class ShannonDetector:
    '''
    Identifies regions of low complexity in a protein sequence. 
    Algorithm:
        0. Slicing protein sequence to subsequences using the sliding window method
        1. Calculates the Shannon entropy of subsequences. 
        2. Finds subsequences with entropy lower than random — low complexity points (LCPs).
        3. Low-complexity points are converted into a sequence of amino acids and trimmed.
        4. Low-complexity regions (LCRs) are formed by merging.
    '''

    def __init__(self, sequence: str, window: int, min_length: int):
        '''
        Initializes the Shannon entropy-based low complexity region detector.
        '''
        self.sequence = sequence
        self.window = window
        self.min_length = min_length
        self.rand_entropy = entropy([1/self.window]*self.window)
        self.lcr: List[LCR] = []
                
    def determine_lcr(self, lcp):
        '''
        Converts contiguous Low Complexity Points into Low Complexity Blobs (LCBs).
        
        Groups consecutive LCPs and extends them by window size to form continuous
        amino acid sequences. Each LCB represents a candidate low-complexity region.
        '''
        for _, group in groupby(enumerate(lcp), lambda x: x[0] - x[1]):
            lcp_group = [pos for _, pos in group]
            start = lcp_group[0] 
            end = lcp_group[-1] + self.window - 1
            seq = self.sequence[start: end + 1]
            self.lcr.append(self._trimming_lcr(LCR(seq, start, end)))
    
    def _trimming_lcr(self, lcr: LCR) -> LCR:
        '''
        Trims Low Complexity Regions by removing flanking unique residues.
        '''
        seq, start, end = lcr.seq, lcr.start, lcr.end
        if seq.find(seq[0], 1) == -1:
            seq, start = seq[1:], start + 1
        if seq.find(seq[-1], 0, -1) == -1:
            seq, end = seq[:-1], end - 1
        return LCR(seq, start, end)
    
    def merge_lcr(self):
        '''
        Merges overlapping low complexity regions.
        '''
        i = 0
        while i < len(self.lcr) - 1:
            cur_lcr, nxt_lcr = self.lcr[i], self.lcr[i + 1]
            if cur_lcr.end - nxt_lcr.start == 0:
                self.lcr[i] = LCR(self.sequence[cur_lcr.start: nxt_lcr.end + 1], cur_lcr.start, nxt_lcr.end)
                del self.lcr[i + 1]
            else:
                i += 1

    def filter_by_len(self):
        '''
        Removes LCRs shorter than the configured minimal length.
        '''
        self.lcr = [lcr for lcr in self.lcr if lcr.end - lcr.start + 1 >= self.min_length]

    def run(self):
        '''
        Executes the complete low complexity region detection pipeline
        '''
        entropies = seq2smth.seq2entropy(self.sequence, self.window)
        lcp = [i for i, aa_entropy in enumerate(entropies) if aa_entropy < self.rand_entropy]
        self.determine_lcr(lcp)
        self.merge_lcr()
        self.filter_by_len()

class DotplotDetector:
    '''
    Identifies regions of low complexity in a protein sequence using a dotplot
    self-similarity matrix with 2D kernel convolution and an FDR-based threshold.
    Algorithm:
        0. Builds a binary dotplot: 1 where amino acids match, 0 otherwise.
        1. Convolves the dotplot with a uniform kernel (kernel_size x kernel_size),
           so each cell of the result holds the local match density. Only the main
           diagonal is used, because it reflects self-similarity of the sequence
           with itself.
        2. Estimates the null distribution: shuffles the sequence, recomputes the
           convolved diagonal, and pools the values across iterations. The
           (1 - fdr_rate) percentile becomes the detection threshold.
        3. Marks diagonal positions with density at or above the threshold as low
           complexity points (LCPs) and groups contiguous ones into low complexity
           regions (LCRs).
        4. Filters LCRs shorter than min_length.
    '''

    def __init__(self, sequence: str, kernel_size: int, fdr_rate: float = 0.01, min_length: int = 4):
        '''
        Initializes the dotplot-based low complexity region detector.
        '''
        self.sequence = sequence
        self.kernel = np.ones((kernel_size, kernel_size), dtype=np.float32)
        self.fdr_rate = fdr_rate
        self.min_length = min_length

        self.array = np.array(list(sequence))
        self.dotplot = (self.array[:, None] == self.array[None, :]).astype(np.float32)
        self.convolved = convolve2d(self.dotplot, self.kernel, mode='same', boundary='fill')

        self.threshold: Optional[float] = None
        self.lcr: List[LCR] = []

    def calculate_fdr_threshold(self, iterations: int = 50):
        '''
        Estimates a false discovery rate (FDR) threshold from shuffled sequences.

        Each iteration shuffles the input sequence, rebuilds the convolved dotplot,
        and collects its main diagonal. All collected diagonals are pooled, and the
        (1 - fdr_rate)-th percentile is stored as the detection threshold.
        '''
        sequence_list = list(self.sequence)
        diag = []
        for _ in range(iterations):
            random.shuffle(sequence_list)
            shuffled_sequence = ''.join(sequence_list)
            shuffled_array = np.array(list(shuffled_sequence))
            shuffled_dotplot = (shuffled_array[:, None] == shuffled_array[None, :]).astype(np.float32)
            shuffled_convolved = convolve2d(shuffled_dotplot, self.kernel)
            diag.append(np.diag(shuffled_convolved))
        diag = np.concatenate(diag)

        self.threshold = float(np.percentile(diag, 100 * (1 - self.fdr_rate)))

    def determine_lcr(self):
        '''
        Converts contiguous Low Complexity Points on the main diagonal into
        Low Complexity Regions (LCRs).

        Diagonal cells whose density is at or above the threshold are marked as
        LCPs. Contiguous LCPs are grouped and converted into LCR objects.
        '''
        lcp = np.diag(self.convolved)
        mask = lcp >= self.threshold

        start = None
        for i, is_lcr in enumerate(mask):
            if is_lcr and start is None:
                start = i
            elif not is_lcr and start is not None:
                self.lcr.append(LCR(self.sequence[start: i], start, i - 1))
                start = None

        if start is not None:
            self.lcr.append(LCR(self.sequence[start: len(mask)], start, len(mask) - 1))

    def filter_by_len(self):
        '''
        Removes LCRs shorter than the configured minimal length.
        '''
        self.lcr = [lcr for lcr in self.lcr if lcr.end - lcr.start + 1 >= self.min_length]

    def run(self):
        '''
        Executes the complete low complexity region detection pipeline.
        '''
        self.calculate_fdr_threshold(fdr_rate=self.fdr_rate)
        self.determine_lcr()
        self.filter_by_len()