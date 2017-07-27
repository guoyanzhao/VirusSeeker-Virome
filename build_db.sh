#!/bin/bash -x

echo "CREATE TABLE gi_taxid_nucl (gi integer PRIMARY KEY, tax_id integer);" | sqlite3 vhunter.db
echo -e '.separator "\t"\n.import /scratch/dwlab/databases/taxdump_20141021/gi_taxid_nucl.dmp gi_taxid_nucl\n'| sqlite3 vhunter.db
