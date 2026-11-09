from typing import Dict, List, Literal

from scipy.stats import entropy


Charge = Literal[-1, 0, 1]
AAType = Literal[0, 1, 2, 3, 4]

def seq2freq(sequence: str) -> Dict[str, float]:
    '''
    Calculate amino acid frequencies in a protein sequence.
    '''
    letters = ('D', 'K', 'I', 'Y', 'G', 'R', 'M', 'E', 'L', 'W', 'P', 'F', 'H', 'T', 'N', 'A', 'C', 'V', 'S', 'Q')
    seq_len = sum(sequence.count(aa) for aa in letters)
    return {aa: sequence.count(aa) / seq_len for aa in letters}

def seq2charge(sequence: str) -> List[Charge]:
    '''
    Convert a protein sequence to a sequence of charge sequence.
    '''
    aa2charge = {
        'S': 0, 'A': 0, 'V': 0, 'L': 0, 'I': 0, 'M': 0, 'P': 0,'G': 0, 'T': 0, 'C': 0, 'Q': 0, 'N': 0, 'F': 0, 'Y': 0, 'W': 0,
        'K': 1, 'R': 1, 'H': 1, 
        'E': -1, 'D': -1
    }
    return [aa2charge[aa] for aa in sequence if aa in aa2charge]

def seq2type(sequence: str) -> List[AAType]:
    '''
    Convert a protein sequence to a sequence of amino acid types.
    '''
    aa2type = {
        'A': 0, 'V': 0, 'L': 0, 'I': 0, 'M': 0, 'P': 0,
        'G': 1, 'T': 1, 'C': 1, 'Q': 1, 'N': 1, 'S': 1, 
        'E': 2, 'D': 2, 
        'K': 3, 'R': 3, 'H': 3,
        'F': 4, 'Y': 4, 'W': 4
    }
    return [aa2type[aa] for aa in sequence if aa in aa2type]

def seq2slices(sequence: str, window: int) -> List[str]:
    '''
    Generate sliding window slices of a protein sequence.
    '''
    upper_bound = len(sequence) - window + 1
    return [sequence[l_:l_ + window] for l_ in range(upper_bound)]

def seq2entropy(sequence: str, window: int) -> List[float]:
    '''
    Calculate Shannon entropy for sliding window slices of a protein sequence.
    '''
    slices = seq2slices(sequence, window)
    return [entropy(list(freqs.values())) for freqs in (seq2freq(seq) for seq in slices)]