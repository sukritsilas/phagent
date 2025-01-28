rm -r temp/
mkdir mmseqs_db
mkdir mmseqs_out
mv SPFMphages_CDS_mmdb* mmseqs_db
mv SPFMphages_CDS_clustered.* mmseqs_out
cd mmseqs_out
mv *.tsv ../
cd ..
