#These notes contain the pathfinding commands for parsing (specifically) the Lamprey gtf for the relevant orthogroup crossreference IDs. To get mapping files for the other model taxa, replace "HsTFAA__v__lamprey_protein_coding_merged.parsed.tsv" with the tsv for that model, as well as the "lamprey_protein_coding_merged.fasta.transdecoder.genome_gene_name_MT.gtf" with the correct GTF.

#WARNING: Some of the sed commands (like sed 's/-RA\.p[0-9]//g'), were written specifically to trim extensions within the lamprey ortho IDs, and may not perfectly apply to the names in the files for the other model organisms.

#Get sorted, unique list of Lamprey Accessions (including "PMZ" and "MSTG" IDs)
awk '{print $4}' HsTFAA__v__pmz.proteins.parsed.tsv | sed 's/,/\n/g' | sed 's/-P//g' | sed 's/\./\t/g' | awk '{print $1"       "$2}' | sed 's/\t/\./g' | sed 's/\.$//g' | grep -v "amprey" | sort | uniq > LampreyAccessionsPmz

#Take list of accessions and find them in the "lamprey" GTF. Parse results for just the value in the "gene_name" position
>CrossRefLinesInGFFPmz; cat LampreyAccessionsPmz | while read accession; do crossID=`grep $accession pmz.genes.chr.names.gff | grep "Name" | head -n 1 | awk '{print $9}' | sed 's/;/\t/g' | awk '{print $2}' | sed 's/Name=//g'`; echo "$accession	$crossID" >> CrossRefLinesInGFFPmz; done


#Create an almost-finished cross-reference file that contains the old and new values from lamprey and blast results
>DummyTempFilePmz; cat CrossRefLinesInGFFPmz | while read old new; do thirdandfourthcolumn=`grep $old HsTFAA__v__pmz.proteins.parsed.tsv | head -n 1 | awk '{print $3"	"$4}'`; echo "$old	$new	$thirdandfourthcolumn" >> DummyTempFilePmz; done

#Create a lise where hits with multiple PmzIDs are expaned to one row
awk '{print $4}' DummyTempFilePmz | sed 's/,/\n/g' | sort | uniq > ExpandedCSVListPmz

>OneLinePerIDPmz; cat ExpandedCSVListPmz | while read line; do cols123=`grep "$line" DummyTempFilePmz | awk -v var="$line" '{print $1"	"$2"	"$3"	"var}'`; echo "$cols123" >> OneLinePerIDPmz; done


#Put each PMZ and MSTG term on its own line that contains the short-form name of the ortholog
>LampreyMappingFilePmz.tsv; cat ExpandedCSVList | while read line; do cols123=`grep "$line" DummyTempFilePmz | awk -v var="$line" '{print $1"	"$2"	"$3"	"var}'`; echo "$cols123" >> LampreyMappingFilePmz.tsv; done
