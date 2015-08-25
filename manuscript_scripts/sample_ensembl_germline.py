import os, cPickle

import sqlalchemy as sql
from cogent.db.ensembl import Genome, HostAccount

# GDUSERV account details, makes execution that much faster
account = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
human = Genome('human', Release=79, account=account, pool_recycle=10000)

# sample all non-somatic SNPs whose reported flanks match the reference
# genome
variation_feature = human.VarDb.getTable("variation_feature")

condition = sql.and_(variation_feature.c.somatic==0,
            variation_feature.c.alignment_quality==1)

query = sql.select([variation_feature], condition)

# for debugging
outfile_name = "../data/ensembl_snps_79/GERMLINE_flanks_match_ref.txt"
outfile = open(outfile_name, 'w')

# use limit() and offset() to specify a range of records
# query = query.limit(3)
for record in query.execute():
    out = cPickle.dump(record, outfile)
outfile.close()

print
print "DONE!"
print "Results written to", outfile_name
