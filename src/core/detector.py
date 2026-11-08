from itertools import groupby
from typing import NamedTuple, List

from scipy.stats import entropy

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