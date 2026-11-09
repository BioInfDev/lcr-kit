import csv
from pathlib import Path
from typing import Annotated

import typer
from Bio import SeqIO

from lcr_kit.core.detector import ShannonDetector

app = typer.Typer()


@app.command()
def detector(fasta: Annotated[Path, typer.Argument(exists = True, dir_okay = False, help = 'Path to the FASTA with protein sequence')],
            protein: Annotated[Path, typer.Argument(help = 'Gene name or Protein name')],
            window: Annotated[int, typer.Option('--window', '-w', help = 'Sliding window size')] = 3,
            min_length: Annotated[int, typer.Option('--min_length', '-l', help = 'Minimal LCR length')] = 4):
    record = next(SeqIO.parse(str(fasta), 'fasta'))

    shannon_detector = ShannonDetector(sequence=record.seq, window=window, min_length=min_length)
    shannon_detector.run()

    header = [f"Low-complexity regions of {protein} | Method: 'Shannon Entropy' (window = {window}) | Minimum LCR length = {min_length}"]
    columns = ['Sequence', 'Start', 'End']
    lcr = [[lcr.seq, lcr.start, lcr.end] for lcr in shannon_detector.lcr]

    with open(f'lcr_shannon_{protein}.tsv', 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f,  delimiter='\t')
        writer.writerows([header, columns, *lcr])    


if __name__ == '__main__':
    app()