import pandas as pd
import math

Input = "Variants_CLC_table.txt"
GFF3_FILE ="sequence.gff3"
List_of_Genes = ['Gene_01_accession.1','Gene_02_accession.1','Gene_03_accession.1','Gene_04_accession.1']
List_of_Genes_end = [item + "_end" for item in List_of_Genes]


## Format feature table to pandas DF.................................................
column_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
FT = pd.read_csv(GFF3_FILE, sep="\t", quoting=3, skiprows = 5, names= column_names)
## Clean up DF by removing unnecessary columns and expanding the attributes column
FT = FT.drop(columns=['source', 'score'])
FT[['ID', 'B']] = FT['attributes'].str.split(';', n=1, expand=True)
FT = FT.drop(columns = ['B', 'attributes'])
FT['ID'] = FT['ID'].str.replace("ID=", "")
FT['ID'] = FT['ID'].str.replace("cds-", "")
FT[["start", "end"]] = FT[["start", "end"]].apply(pd.to_numeric)

## Format variant Table...............................................................
DF = pd.read_csv(Input, sep = "\t")
DF[['A', 'Position', 'C']] = DF['PosMutation'].str.split('_', n=2, expand=True)
DF = DF.drop(columns = ["A", "C"])
## Replace .. and ^ with a temp variable and then split column into two
DF['Position'] = DF['Position'].str.replace("^", "#")
DF['Position'] = DF['Position'].str.replace("..", "#")
DF[['StartPosition', 'EndPosition']] = DF['Position'].str.split("#", n=1, expand=True)
DF[["StartPosition", "EndPosition"]] = DF[["StartPosition", "EndPosition"]].apply(pd.to_numeric)
## Replace NA's with 0, needed for later calculations
DF = DF.fillna(0)

## Create CDS table...................................................................
CDS_DF = FT[FT["type"] == 'CDS'].reset_index()

## check if starting and ending postion are within range of each gene........................
## and add a new column with gene name and TRUE or FALSE
for row in CDS_DF.itertuples():
    DF[row.ID] = DF['StartPosition'].between(int(row.start), int(row.end), inclusive="both")
    DF[row.ID + "_end"] = DF['EndPosition'].between(int(row.start), int(row.end), inclusive="both")

## Set tables and variables
DF['Codon'] = 0
DF['Codon_end'] = 0
DF['CDS'] = ""

for Gene in List_of_Genes:
    ## Subtract gene start from variant and subtract 3
    Gene_Start = CDS_DF[CDS_DF['ID'] == Gene]['start'].max() - 3
    ## Perform this calcuation only if the variant is within the coding region of the gene
    DF['Codon'] = DF.apply(lambda row: math.floor((row['StartPosition'] - Gene_Start) / 3) if row[Gene] else row['Codon'], axis=1)
    DF['Codon_end'] = DF.apply(lambda row: math.floor((row['EndPosition'] - Gene_Start) / 3) if row[Gene + "_end"] else row['Codon_end'], axis=1)
    ## create another column that prints the gene name
    DF['CDS'] = DF.apply(lambda row: str(Gene) if row[Gene] else row['CDS'], axis = 1)

## Clean up columns and empty datapoints that are no longer needed
DF= DF.drop(columns = List_of_Genes_end)
DF= DF.drop(columns = List_of_Genes)
DF['Codon']= DF['Codon'].replace(0, pd.NA)
DF['Codon_end']= DF['Codon_end'].replace(0, pd.NA)


## Create mRNA table...................................................................
mRNA_DF = FT[FT["type"] == 'mRNA'].reset_index()

for row in mRNA_DF.itertuples():
    DF[row.ID] = DF['StartPosition'].between(int(row.start), int(row.end), inclusive="both")
    DF[row.ID + "_end"] = DF['EndPosition'].between(int(row.start), int(row.end), inclusive="both")
mRNA_list = mRNA_DF['ID'].values.tolist()
DF['mRNA'] = ""
for mRNA in mRNA_list:
    DF['mRNA'] = DF.apply(lambda row: str(mRNA) if row[mRNA] else row['mRNA'], axis = 1)

## Clean up
DF= DF.drop(columns = mRNA_list)
mRNA_list_end = [item + "_end" for item in mRNA_list]
DF= DF.drop(columns = mRNA_list_end)
DF = DF.drop(columns = ['Position'])
DF['EndPosition']= DF['EndPosition'].replace(0, pd.NA)
DF = DF.rename({'EndPosition':'second nucleotide position'}, axis=1)
DF = DF.rename({'StartPosition':'nucleotide position'}, axis=1)
DF.to_csv('firstoutput.txt', sep = "\t", index= False)