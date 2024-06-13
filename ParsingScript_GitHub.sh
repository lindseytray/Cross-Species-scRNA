#Script to convert the OrthoFinder output with multiple EMBL IDs per Orthogroup into one per line
#Usage:
#bash ParsingScript.sh HsTFAA__v__lamprey_protein_coding_merged.linear.tsv

orthofinderInfile=$1

sed 's/ //g' $orthofinderInfile | awk '{print $1"        "$2}' | sed 's/,/\n.\t/g' | awk '{print $2}' > EMBL_Terms

echo "Orthogroup	EMBL_ID	GeneSymbol	LampreyTranscript" > ${orthofinderInfile%.linear.tsv}.parsed.tsv
cat EMBL_Terms | while read embl; do

	OG=`grep $embl $orthofinderInfile | awk '{print $1}'`
	geneName=`grep $embl location/of/HsTFAA.fasta | head -n 1 | sed 's/ /\n/g' | grep "gene_symbol" | sed 's/:/\t/g' | awk '{print $2}'`
	lampreyInfo=`grep $embl $orthofinderInfile | sed 's/ //g' | awk '{print $3}'`
        echo "$OG       $embl	$geneName	$lampreyInfo" >> ${orthofinderInfile%.linear.tsv}.parsed.tsv
	sed -i 's/		/	.	/g' ${orthofinderInfile%.linear.tsv}.parsed.tsv

done
rm EMBL_Terms
