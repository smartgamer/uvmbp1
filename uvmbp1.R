
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
bap1mutsamples=unique(bap1mut2$Tumor_Sample_Barcode)  #samples that have bap1 mutations


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
rnat=as.data.frame(rnat)
head(rnat[1:5,1:5])

#add one column:bap1 status. "wt" or "mut"
rnat$bap1Status=c("wt") 
head(rnat[1:5,1:5])
rnat$bap1Status
which(bap1mutsamples[1] %in% row.names(rnat))
bap1mutsamples[1] %in% row.names(rnat)
library(stringr)
bap1mutsamples2=str_replace_all(bap1mutsamples, "-", ".") #change to the same format
bap1mutsamples2[1] %in% row.names(rnat)
which(bap1mutsamples2[1] %in% row.names(rnat))
mutInd=which(row.names(rnat) %in% bap1mutsamples2)
rnat[mutInd,]$bap1Status="mut"
rnat$bap1Status
library(tidyverse)
rnat2=rnat 
rnat2$bap1Status=ifelse(rnat2$bap1Status=="wt", 0, 1)
rnat2$bap1Status
head(rnat2[1:5,1:5])

#train rain forest model
library(caret)
library(randomForest)
set.seed(1001)
# fit_rf = randomForest(bap1Status~., data=rnat2)
#Error in terms.formula(formula, data = data) : duplicated name 'UBE2Q2P3' in data frame using '.'
rnat3=rnat2[,!duplicated(colnames(rnat2))]  #https://stackoverflow.com/questions/24142942/how-to-remove-duplicated-column-names-in-r

nzv=nearZeroVar(rnat3)
rnat4=rnat3[,-nzv]
nearZeroVar(rnat4)
save(rnat, rnat2, rnat3, rnat4, file = "rnat.RData") 

# correlations = cor(rnat4[,1:18354])
correlations = cor(rnat4[,1:5000])
dim(correlations) 
highCorr = findCorrelation(correlations, cutoff=0.75)
rnat5 = rnat4[,-highCorr]
#repeat
correlations = cor(rnat5[,2200:8000])
highCorr = findCorrelation(correlations, cutoff=0.75)
rnat5 = rnat5[,-highCorr]

correlations = cor(rnat5[,6000:12000])
highCorr = findCorrelation(correlations, cutoff=0.75)
rnat5 = rnat5[,-highCorr]

correlations = cor(rnat5[,6000:10936])
highCorr = findCorrelation(correlations, cutoff=0.75)
rnat5 = rnat5[,-highCorr]

correlations = cor(rnat5[,1:8556])
highCorr = findCorrelation(correlations, cutoff=0.75)
rnat6 = rnat5[,-highCorr]

"bap1Status" %in% colnames(rnat6)
head(rnat6[, 4359:4366])
save(rnat6, file="rnat6.RData")


fit_rf = randomForest(bap1Status~., data=rnat3) #Error: protect(): protection stack overflow
# options(expressions = 5e5)
# ?options
# fit_rf = randomForest(bap1Status~., data=rnat3[,20300:20501])
#https://stackoverflow.com/questions/12767432/how-can-i-tell-when-my-dataset-in-r-is-going-to-be-too-large


rnat6$"ANKHD1-EIF4EBP3"   #Error: object 'EIF4EBP3' not found
rnat6$"ANKHD1-EIF4EBP3"
rnat6$bap1Status

fit_rf = randomForest(bap1Status~., data=rnat6) #Error in eval(predvars, data, env) : object 'ANKHD1-EIF4EBP3' not found
head(colnames(rnat6), 220)
# colnames(rnat6)=str_replace_all(colnames(rnat6), "-", "_") #remove -
# colnames(rnat6)=str_replace_all(colnames(rnat6), "/", "_")
# colnames(rnat6)=str_replace_all(colnames(rnat6), "~", "_")
# colnames(rnat6)=str_replace_all(colnames(rnat6), "(", "_")
# colnames(rnat6)=str_replace_all(colnames(rnat6), ")", "_")

nam=paste("gene", as.character(seq_along(1:4365)), sep = "", collapse = NULL) #colnames have many symbols that are not allowed by random forest
nam=c(nam, "bap1Status")
colnames(rnat6)=nam
tail(colnames(rnat6))

fit_rf = randomForest(bap1Status~., data=rnat6)
# fit_rf = randomForest(bap1Status~., data=rnat6[, 4300:4366]) 
# rnat6$bap1Status=as.factor(rnat6$bap1Status) #change numeric to factor
# fit_rf = randomForest(bap1Status~., data=rnat6[, 4300:4366]) 
# fit_rf = randomForest(bap1Status~., data=rnat6[, 3500:4366]) 


fit_rf
# Create an importance based on mean decreasing gini
importance(fit_rf)
varImp(fit_rf)

# Create a plot of importance scores by random forest
varImpPlot(fit_rf)


colnames(rnat)
load("~/mygit/uvmbp1/rnat6.RData")
head(rnat6[, 4359:4366])
colnames(rnat6[c(2432, 3430,3310,3938,2368,2176,2776,3109,3976,2003,1213,1835)]) 
# [1] "SNF1LK" "TMEM25" "TIL"    "VDAC1"  "SFXN3"  "CCO"    "SNX22"  "TAS1R1" "VSTM2L" "RIBC1"  "NR6A1" 
# [12] "RADIL" 
save(fit_rf, file="fit_rf.RData")


genenames=colnames(rnat6) #the last column is bap1Status
colnames(rnat6)=nam
str(rnat6$bap1Status)
rnat6$bap1Status=as.factor(rnat6$bap1Status)
str(rnat6$bap1Status)
#use caret to improve model performance
library(caret)
ctrl = trainControl(method="repeatedcv", number=10, repeats=10)
grid_rf=expand.grid(.mtry=c(2,4,8,16,32))
set.seed(1001)
m_rf=train(bap1Status ~., data = rnat6, method="rf", metric="Kappa", trControl=ctrl, tuneGrid=grid_rf )

#C5.0
grid_c50=expand.grid(.model = "tree",
  .trials=c(2,4,8,16,32),
  .winnow = "FALSE")
set.seed(1001)
m_c50=train(bap1Status ~., data = rnat6, method="C5.0", metric="Kappa", trControl=ctrl, tuneGrid=grid_c50 )
warnings()

m_rf
save(m_rf, file="m_rf.RData")
varImp(m_rf)
varImp(m_rf, scale=F)
plot(varImp(m_rf,scale=F), top = 20)
genenames[c(3719,1998,2263,1926,372,4349,2520,3025,3976,872,3938,1312,1093,3283,2384,3479,875,1057,2003,3638)]
# [1] "TTC21A"    "RHOF"      "SDHA"      "RCOR2"     "HOGA1"     "ZSCAN1"    "SLC25A30"  "SUMF2"    
# [9] "VSTM2L"    "LOH3CR2A"  "VDAC1"     "PAIP2B"    "LOC158856" "TIMELESS"  "SGSM2"     "TMPPE"    
# [17] "LPAR1"     "MYEOV"     "RIBC1"     "TRIM65"   

m_c50
save(m_c50, file="m_c50.RData")
varImp(m_c50)
varImp(m_c50,scale=F)
plot(varImp(m_c50,scale=F), top = 20)

#https://topepo.github.io/caret/variable-importance.html


