import os
import sys

import click
from cogent.db.ensembl import HostAccount, Genome
from scitrack import CachingLogger

__author__ = "Bob Buckley, Gavin Huttley"
__copyright__ = "Copyright 2015, Bob Buckley"
__credits__ = ["Bob Buckley", "Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

LOGGER = CachingLogger(create_dir=False)

account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split()) if 'ENSEMBL_ACCOUNT' in os.environ else None

chnames = tuple(map(str, range(1, 23)))+('X', 'Y', 'MT')

# for Click documentation see http://click.pocoo.org/
@click.command()
@click.option('--release', required=True, type=int, default=79, help='Ensembl release.')
@click.option('--output_path', required=True, default=None, help='Path to write output to file.')
def dumpchroms(release, output_path):
    """dump the chromosome sequences to separate files.

    output to files Ch??_n.fa - ?? is COSMIC 2 digit name/number used as ID
    and n is the name in the database.
    """
    global account, chnames
    if not os.path.isdir(output_path):
        raise ValueError("Create the output dir first")
    
    runlog_path = os.path.join(output_path, 'dump_genome.log')
    LOGGER.log_file_path = runlog_path
    
    human = Genome(Species='human', Release=release, account=account)
    for i,n in enumerate(chnames, start=1):
        nx = "Chr%02d_%s" % (i, n) + ".fa"
        print "dumping chromosome file", nx
        print "   fetching ..."
        sys.stdout.flush()
        ch = human.getRegion(CoordName=n)
        ch.Seq.Name = "Chr" + n
        print "   dumping ..."
        sys.stdout.flush()

        # output to temporary file name ... then rename once output completes
        # .. safer when working with large files.
        xoutpath = os.path.join(output_path, "z-"+nx)
        outpath  = os.path.join(output_path, nx)
        with open(xoutpath, "w") as dst:
            dst.write(ch.Seq.toFasta()+"\n")
        os.rename(xoutpath, outpath)
        LOGGER.output_file(outpath)
        
    return

if __name__=="__main__":
    dumpchroms()
    sys.stdout.flush()
    print "Done."
