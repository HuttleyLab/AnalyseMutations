import logging, sys, platform
from traceback import format_exc
from socket import gethostname
import hashlib

RUN_VARS_TOKEN = "vars:"
DELIM = ' : '

class CachingLogger(object):
    """stores log messages until a log filename is provided"""
    def __init__(self, log_file_path=None):
        super(CachingLogger, self).__init__()
        self.log_file_path = log_file_path
        self._started = False
        self._messages = []
    
    def write(self, msg):
        """writes a log message"""
        
        if not self._started:
            self._messages.append(msg)
        else:
            logging.info(msg)
        
        if not self._started and self.log_file_path:
            # start the logger and flush the message cache
            set_logger(self.log_file_path)
            for m in self._messages:
                logging.info(m)
            
            self._messages = []
            self._started = True
            
    

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