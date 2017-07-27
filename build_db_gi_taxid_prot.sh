#!/bin/bash -x

echo "CREATE TABLE gi_taxid_prot (gi integer PRIMARY KEY, tax_id integer);" | sqlite3 vhunter.db
echo -e '.separator "\t"\n.import /scratch/dwlab/databases/taxdump_20160802/gi_taxid_prot.dmp gi_taxid_prot\n'| sqlite3 vhunter.db
