import csv
from pathlib import Path
from typing import Annotated

import typer
from Bio import SeqIO

from lcr_kit.core.detector import DotplotDetector

app = typer.Typer()


@app.command()
def detector(fasta: Annotated[Path, typer.Argument(exists = True, dir_okay = False, help = 'Path to the FASTA with protein sequence')],
            protein: Annotated[Path, typer.Argument(help = 'Gene name or Protein name')],
            kernel_size: Annotated[int, typer.Option('--kernel-size', '-k', help = 'Kernel size for convolution')] = 3,
            fdr_rate: Annotated[float, typer.Option('--fdr-rate', '-t', help = 'Fdr rate')] = 0.01,
            min_length: Annotated[int, typer.Option('--min_length', '-l', help = 'Minimal LCR length')] = 4):
    record = next(SeqIO.parse(str(fasta), 'fasta'))

    dotplot_detector = DotplotDetector(sequence=record.seq, kernel_size=kernel_size, fdr_rate=fdr_rate, min_length=min_length)
    dotplot_detector.run()

    header = [f"Low-complexity regions of {protein} | Method: 'Dotplot-Convolve-FDR' (kernel = {kernel_size}x{kernel_size}, fdr-rate = {fdr_rate}) | Minimum LCR length = {min_length}"]
    columns = ['Sequence', 'Start', 'End']
    lcr = [[lcr.seq, lcr.start, lcr.end] for lcr in dotplot_detector.lcr]

    with open(f'lcr_dotplot_{protein}.tsv', 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f,  delimiter='\t')
        writer.writerows([header, columns, *lcr]) 


if __name__ == '__main__':
    app()