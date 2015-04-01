#!/bin/bash
##################################################
## Kyle Hernandez
## Copyright 2014, Center for Research Informatics
## University of Chicago
##################################################

function usage(){
    cat<<EOF
--------------------------------------------------------------------------------
$(basename $0): A tool to format the COSMIC VCF file for ExScaliburSMD
Copyright 2014, Center for Research Informatics <cri.uchicago.edu>
University of Chicago
--------------------------------------------------------------------------------
usage: -i <CosmicCodingMuts.vcf> -o <FormattedCosmicMuts.vcf> -d <ucsc.hg19.dict> 

-i	Input CosmicCodingMuts.vcf file as downloaded from COSMIC
-o	Output formatted COSMIC VCF file
-d	Path to the GATK hg19 bundle's sequence dictionary file
EOF
}

# Globals
INVCF=
OVCF=
REF=

# Define options
while getopts "hi:o:d:" OPTION
do
    case $OPTION in
        h)
            usage
            exit 0
            ;;
        i)
            INVCF=$OPTARG
            ;;
        o)
            OVCF=$OPTARG
            ;;
        d)
            REF=$OPTARG
            ;;
        ?)
            usage
            exit 0
            ;;
    esac
done

# Test for required args
if [[ -z $INVCF ]]||[[ -z $OVCF ]]||[[ -z $REF ]]; then
    usage
    exit 1
fi

if [ ! -e $INVCF ]; then echo "Missing input file $INVCF"; exit 1; fi
if [ ! -e $REF ]; then echo "Missing sequence diciontary file $REF"; exit 1; fi

# TMP file 
TMP="${INVCF}.tmp"

# Add the chr to the beginning of the line
echo "Formatting chromosome names..." 
trap 'echo "ERROR!!"; exit 1' ERR
set -o pipefail
cat $INVCF | awk '{if($1 ~ /^##reference/){print "##reference=hg19"}else if($1 ~ /^#/){print $0}else if($1~/^MT/){sub(/^MT/,"chrM"); print $0}else{print "chr" $0}}' > $TMP

# Sort and don't create index because it isn't compatible with MuTect
echo "Sorting the new VCF file..."
java -Xmx4G -jar ${PICARD}/SortVcf.jar CREATE_INDEX=false I=$TMP O=$OVCF SD=$REF

# Cleanup
echo "Removing temporary files..."
rm $TMP
