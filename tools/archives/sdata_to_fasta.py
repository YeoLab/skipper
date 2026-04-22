from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import seqdata as sd

def write_fasta(test_sdata, output_filepath, mask = 32):
    """
    Write a FASTA file from an iterable of objects with .name and .seq attributes.
    
    Parameters
    ----------
    test_sdata
        seqdata object with attribute seq and name
    output_filepath : str
        Path to write the FASTA file to.
    """
    records = []
    for i in range(len(test_sdata['name'])):
        record = SeqRecord(Seq(''.join([c.decode() for c in test_sdata['seq'][i].values.tolist()][mask:-mask])), 
                           id=str(test_sdata['name'][i].values.tolist()), description="")
        records.append(record)
    # Write all records in one go
    SeqIO.write(records, output_filepath, "fasta")

if __name__ == '__main__':
    input_sdata = sys.argv[1]
    test_sdata = sd.open_zarr(input_sdata).load()
    outf = sys.argv[2]

    write_fasta(test_sdata, outf)