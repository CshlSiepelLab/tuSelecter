# tuSelecter

Warning: This application is somewhat memory intensive. Reccomend on platform with no less than 16GB of RAM. 

Required R Packages:
argparse
data.table
depmixS4
rtracklayer
doParallel

This application requires a preprocessed form of transcripts obtained from a gtf file. The build_txtable.sh script will do this for you.

./build_txtable.sh gencode.v19.annotation.gtf > txtable.out

Then you run the tuSelecter.R program.
For more information run: 
tuSelecter.R --help 



 