import logging, sys, platform
from traceback import format_exc
from socket import gethostname
import hashlib

RUN_VARS_TOKEN = "vars:"
DELIM = ' : '


def set_logger(logfile_name, level=logging.DEBUG):
    """setup logging"""
    logging.basicConfig(filename=logfile_name, filemode='w', level=level,
                format='%(asctime)s\t%(levelname)s\t%(message)s',
                datefmt="%Y-%m-%d %H:%M:%S")
    logging.info('system_details: system=%s:python=%s' % (platform.version(),
                                                  platform.python_version()))

def get_file_hexdigest(filename):
    '''returns the md5 hexadecimal checksum of the file'''
    # from http://stackoverflow.com/questions/1131220/get-md5-hash-of-big-files-in-python
    with open(filename) as infile:
        md5 = hashlib.md5()
        while True:
            data = infile.read(128)
            if not data:
                break
            
            md5.update(data)
    return md5.hexdigest()

if __name__ == "__main__":
    infile = "/Users/gavin/DevRepos/AnalyseMutations/data/aligns/A/intergenic/freq_All-chrom_A-GC_All-AtoT.fasta.gz"
    md5 = get_file_hexdigest(infile)
    print md5