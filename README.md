# phagent
 Pipeline to extract accessory regions from phage genomes

This code will find accessory regions from the Salmonella phages in SPFMphages.gb

You must have biopython (v1.79), MMseqs2 (see References), and pandas (v1.1.5) installed before running this script 

Run the pipeline by entering the following code on the Unix command prompt sh RunMe.sh

This will execute all the steps for finding accessory regions in SPFMphages.gb in order

Non-redundant accessory regions will be in file: SPFMphages_accRegions_nr_all.txt

The last program in RunMe.sh (subsetAccRegionMatrix.sh) subsets an individual accessory
region and extracts it from the various genome genbank files

This program can be modified and used to extract any accessory region in
SPFMphages_accRegions_nr_all.txt, like so --

sh subsetAccRegionMatrix.sh "(503,5737)" 1.691phamily17
Note, no space between anchors 503,5737

By default, RunMe.sh automatically extracts the first accessory region in SPFMphages_accRegions_nr_all.txt
This is 1.691phamily17 between anchors 503,5737. 

You only need to run RunMe.sh once. To extract other accessory regions, run only the
subsetAccRegionMatrix.sh program again with the phamily name and the anchors substituted
appropriately with the desired accessory region from SPFMphages_accRegions_nr_all.txt

To visualize the accessory regions, upload the extracted genbank files into Clinker
https://github.com/gamcil/clinker

To use a different genome dataset (.gb file), replace SPFMphages with the name of your dataset in all .sh files-
RunMe.sh
mmseqsCluster.sh
mmSeqsCleanup.sh
makeFastANI.sh
subsetAccRegionMatrix.sh