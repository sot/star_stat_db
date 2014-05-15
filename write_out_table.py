#!/usr/bin/env python
import os
import sys
import Ska.DBI

def write_table(dbfile, tablename, outfile):
    db = Ska.DBI.DBI(dbi='sqlite', server=dbfile)
    data = db.fetchall("select * from {}".format(tablename))
    cols = data.dtype.names
    out = open(outfile, 'w')
    for row in data:
        for col in cols:
            out.write("{} {} {} {}\n".format(
                    row['obsid'], row['slot'], col, row[col]))
    out.close()

if __name__ == '__main__':
    write_table(sys.argv[1], sys.argv[2], sys.argv[3])
