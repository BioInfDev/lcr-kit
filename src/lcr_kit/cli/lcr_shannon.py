from pathlib import Path
from typing import Annotated
import json

from Bio import SeqIO
import typer
import pandas as pd

from lcr_kit.core.detector import ShannonDetector, LCR


app = typer.Typer()


@app.command()
def detector(fasta: Annotated[Path, typer.Argument(exists = True, dir_okay = False, help = 'Path to the FASTA with protein sequence')],
            output: Annotated[Path, typer.Option('--output', '-o', exists=False, dir_okay=True, help='Path to output')],
            window: Annotated[int, typer.Option('--window', '-w', help = 'Sliding window size')] = 3,
            min_length: Annotated[int, typer.Option('--min_length', '-l', help = 'Minimal LCR length')] = 4):
    if not output.exists():
        output.mkdir()
    
    record = next(SeqIO.parse(str(fasta), 'fasta'))
    record_output = output / '_'.join(record.id.split('|')[1:])
    params_json = output / 'params.json'

    shannon_detector = ShannonDetector(sequence=record.seq, window=window, min_length=min_length)
    shannon_detector.run()

    lcr = pd.DataFrame(shannon_detector.lcr, columns=LCR._fields)
    lcr.to_csv(record_output, sep="\t", index=False)

    params = {'method': 'Shannon Detector', 'window': window, 'min-length': min_length}
    with open(params_json, 'w') as f:
        json.dump(params, f)
    

if __name__ == '__main__':
    app()