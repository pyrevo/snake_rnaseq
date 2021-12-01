library(data.table)

fc <- './Results/featureCounts_genes.txt'

# Rearrange the table
dt <- fread(fc)
dt <- dt[, c(1,12,7,10,9,8,11) ]
colnames(dt) <- c("gene_id","WT1.1","WT2.1","WT3.1","miR.218.2_KO1.1","miR.218.2_KO2.1","miR.218.2_KO3.1")
dt
write.table(dt, file='./Results/featureCounts_genes_mod.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Print out the sampleInfo table
sampleInfo <- data.frame("name"=c("WT1.1", "WT2.1", "WT3.1", "miR.218.2_KO1.1", "miR.218.2_KO2.1", "miR.218.2_KO3.1"), "condition"=c('WT','WT','WT','KO','KO','KO'))
sampleInfo
write.table(sampleInfo, file='./Results/sampleInfo.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Print out the comparison table
comps <- data.frame("comp"=c("KO_vs_WT", "WT_vs_KO"), "cond1"=c("KO", "WT"), "cond2"=c('WT', 'KO'))
comps
write.table(comps, file='./Results/compsTab.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

