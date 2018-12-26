
# get all case names from uvm dataset
uvmCases=read.csv("./data/cases_all.txt", sep="\t", stringsAsFactors = F)
uvmCases=uvmCases[-c(1:4),]
uvmCases[1] ="TCGA-RZ-AB0B-01"
str(uvmCases)
head(uvmCases)

# read mutation data
mutation = read.csv("./data/data_mutations_uniprot.txt", sep="\t", stringsAsFactors = F)
head(mutation)
colnames(mutation)
bap1Index=which(mutation$Hugo_Symbol=="BAP1")
bap1mut=mutation[bap1Index,]
head(bap1mut)
bap1mut2=bap1mut[,c("Hugo_Symbol", "Entrez_Gene_Id", "Tumor_Sample_Barcode","Consequence", "Variant_Classification", "Variant_Type", "Reference_Allele","Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "HGVSp", "HGVSp_Short","Hotspot", "SwissProt_entry_Id","cDNA_Change","Codon_Change", "Protein_Change","HGVS_genomic_change" )]
head(bap1mut2)
row.names(bap1mut2)
bap1mutsamples=unique(bap1mut2$Tumor_Sample_Barcode)


#read RNAseq data
rna = read.csv("./data/data_RNA_Seq_v2_expression_median.txt", sep="\t", stringsAsFactors = F)
head(rna)
head(rna[1:5,1:5])
colnames(rna)
rna2=rna[,-c(1:2)]  #keep all the numeric, otherwise it won't work for transpose(keeping numeric)

rnat=t(rna2)

# rnat=t(rna)
head(rnat[1:5,1:5])
colnames(rnat)=rna[,1]
head(rnat[1:5,1:5])
"BAP1" %in% colnames(rnat)

rm(rnat2)



























































