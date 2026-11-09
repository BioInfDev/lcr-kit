from pathlib import Path
from typing import Annotated
import json

from Bio import SeqIO
import typer
import pandas as pd

from src.core.detector import DotplotDetector, LCR


app = typer.Typer()


@app.command()
def detector(fasta: Annotated[Path, typer.Argument(exists = True, dir_okay = False, help = 'Path to the FASTA with protein sequence')],
            output: Annotated[Path, typer.Option('--output', '-o', exists=False, dir_okay=True, help='Path to output')],
            kernel_size: Annotated[int, typer.Option('--kernel-size', '-k', help = 'Kernel size for convolution')] = 3,
            fdr_rate: Annotated[float, typer.Option('--fdr-rate', '-t', help = 'Fdr rate')] = 0.01,
            min_length: Annotated[int, typer.Option('--min_length', '-l', help = 'Minimal LCR length')] = 4):
    if not output.exists():
        output.mkdir()
    
    record = next(SeqIO.parse(str(fasta), 'fasta'))
    record_output = output / '_'.join(record.id.split('|')[1:])
    params_json = output / 'params.json'

    dotplot_detector = DotplotDetector(sequence=record.seq, kernel_size=kernel_size, fdr_rate=fdr_rate, min_length=min_length)
    dotplot_detector.run()

    lcr = pd.DataFrame(dotplot_detector.lcr, columns=LCR._fields)
    lcr.to_csv(record_output, sep="\t", index=False)

    params = {'method': 'Shannon Detector', 'kernel': f'{kernel_size}x{kernel_size}', 'fdr_rate': fdr_rate,'min-length': min_length}
    with open(params_json, 'w') as f:
        json.dump(params, f)
    

if __name__ == '__main__':
    app()