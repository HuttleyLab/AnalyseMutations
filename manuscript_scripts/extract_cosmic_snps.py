import os
import sys
import gzip

from cytoolz.itertoolz import pluck, get
import click
from cogent import LoadTable
from scitrack import CachingLogger

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2015, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

LOGGER = CachingLogger(create_dir=False)

# for Click documentation see http://click.pocoo.org/
@click.command()
@click.option('--cosmic_raw', required=True, help='Path to COSMIC mutant data.')
@click.option('--output_path', required=True, default=None, help='Path to write COSMIC_SNPs.txt')
def main(cosmic_raw, output_path):
    if not os.path.isfile(cosmic_raw):
        raise ValueError("COSMIC file does not exist at %s" % cosmic_raw)
    
    if not os.path.isdir(output_path):
        raise ValueError("Create the output dir first")
    
    runlog_path = os.path.join(output_path, 'cosmic_export.log')
    LOGGER.log_file_path = runlog_path
    
    LOGGER.input_file(cosmic_raw)
    
    infile = gzip.open(cosmic_raw)
    header = infile.readline().rstrip('\n').split('\t')
    col_indices = map(header.index, ['Mutation ID', 'Mutation genome position',
            'Mutation strand', 'Mutation Description', 'Mutation CDS',
            'Primary histology'])
    
    records = []
    type_index = header.index('Mutation Description')
    is_snp = lambda x: x[type_index].startswith('Substitution')

    for i, line in enumerate(infile):
        line = line.rstrip('\n').split('\t')
        if is_snp(line):
            records.append(get(col_indices, line))
    
    infile.close()
    
    # rename columns for succinctness and consistency
    header = ['Symbol', 'Location', 'Strand', 'Type', 'Change', 'Histology']
    table = LoadTable(header=header, rows=records)

    print repr(table)
    outfile_path = os.path.join(output_path, 'COSMIC_SNPs.txt')
    table.writeToFile(outfile_path, sep='\t')
    
    LOGGER.output_file(outfile_path)
    print "Done!"

if __name__ == "__main__":
    main()
