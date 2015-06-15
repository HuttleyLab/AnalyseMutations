import collections
import csv

__author__ = "Bob Buckley, Gavin Huttley"
__copyright__ = "Copyright 2015, Bob Buckley"
__credits__ = ["Bob Buckley", "Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

def reader(fn, name='NamedTuple', rf=None, ignoreblanklines=True, trimextrafields=True, dialect=None):
    """
    return a generator that returns 'namedtuples' from filename fn.
    The fields named in the rf list are required (rf parameter is required!)
    The file is read using the csv module.
    """
    
    with open(fn, 'rb') as srcfd:
        if not dialect:
            dialect = csv.Sniffer().sniff(srcfd.read(4096))
            srcfd.seek(0)
        src = csv.reader(srcfd, dialect)
        hdr = src.next()
        assert all(x in hdr for x in rf)
        mktuple = collections.namedtuple(name, hdr)
        lenhdr = len(hdr)	
        for rec in src:
            if ignoreblanklines and not len(rec):
                continue
            data = rec[:lenhdr] if trimextrafields and len(rec)>lenhdr else rec
            yield mktuple(*data)
    return

if __name__=="__main__":
    src = reader("test1/manifest.tab", rf='Run download_path TaxID ScientificName SampleName LibraryLayout'.split())
    for rec in src:
        print rec.TaxID, rec.ScientificName
    print 'Done'