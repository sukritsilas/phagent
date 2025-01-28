mkdir fastANImatrix
mv fastANI fastANImatrix/
mv splitGenbanksIntoFastaFiles.py fastANImatrix/
mv parseFastANIresults.py fastANImatrix/
cd fastANImatrix/
python splitGenbanksIntoFastaFiles.py ../SPFMphages.gb 
ls *.fa | sort > SPFMphagesANIlist.txt
./fastANI --ql SPFMphagesANIlist.txt --rl SPFMphagesANIlist.txt -o SPFMphagesANIout.txt
python parseFastANIresults.py SPFMphagesANIout.txt
mv SPFMphagesSimilarityMatrix.p ../
mkdir phageFasta/
mv *.fa phageFasta/
cd ..
