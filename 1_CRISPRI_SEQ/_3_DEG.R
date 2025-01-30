setwd('directory')
nameWT = 'nameWT'
nameMUT = 'nameMUT'
scoresWT = read.csv(paste0('2.5_',nameWT,'_Gene.csv'), header=TRUE)
scoresMUT = read.csv(paste0('2.5_',nameMUT,'_Gene.csv'), header=TRUE)

required_packages <- c("ggrepel", "ggplot2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)}
  library(pkg, character.only = TRUE)}
#-------------------------------------------------------------------------------
scoresWT = scoresWT[!duplicated(scoresWT$gene),]
scoresMUT = scoresMUT[!duplicated(scoresMUT$gene),]
lm_data = merge(scoresWT[,c('target','gene','median','padj')],
                scoresMUT[,c('target','gene','median','padj','operon')],
                by='gene')
colnames(lm_data) = c('gene','WT_target','WT_median','WT_padj',
                      'MUT_target','MUT_median','MUT_padj','Operon')
lm=summary(lm(lm_data$MUT_median~lm_data$WT_median))
remove(pkg,required_packages,scoresWT,scoresMUT)
gene_est <- { lm_data[order(lm_data$gene),]; lm_data$pvalue <- 1; lm_data$estimate <- 0; lm_data }
for (i in 1:dim(gene_est)[1]){
  lm_data$hostfactor=0
  lm_data[lm_data$gene==gene_est$gene[i] & !is.na(lm_data$gene),]$hostfactor=1
  l=summary(lm(lm_data$MUT_median~lm_data$WT_median+lm_data$hostfactor))
  gene_est$pvalue[i]=l$coefficients[3,4]
  gene_est$estimate[i]=l$coefficients[3,1]
}
gene_est$FDR=p.adjust(gene_est$pvalue,method='fdr')
#-------------------------------------------------------------------------------
highlight=c()
suppressor=subset(gene_est, FDR<0.05 & estimate>0)
enhancer=subset(gene_est, FDR<0.05 & estimate<0)
p=ggplot(lm_data,aes(WT_median, MUT_median, label=gene))+
  geom_point(size=0.125,colour='#666666')+
  geom_point(data=suppressor, aes(WT_median,MUT_median), shape=21, stroke=0.15,size=1, fill='#7f7fff', color='black')+
  geom_point(data=enhancer, aes(WT_median,MUT_median), shape=21, stroke=0.15,size=1, fill='#ff7f7f', color='black')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_abline(intercept = lm$coefficients[1,1],
              slope=lm$coefficients[2,1],
              linetype='dashed',
              colour='#414042',
              size=0.5)+
  theme(legend.position="none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),axis.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        panel.background = element_blank())+
  geom_text_repel(data=subset(gene_est, FDR<0.05 | gene %in% highlight ),size=1)
print(p)
pdf(paste0("3.1_",nameWT,"_vs_",nameMUT,"_LinReg.pdf"),
    width=3,
    height=3,
    colormodel='rgb',
    paper = 'A4')
print(p)
dev.off()

#-------------------------------------------------------------------------------
comparison=gene_est
comparison=comparison[order(-comparison$MUT_median),]
finalall=subset(comparison, 
                select = c('gene','Operon','estimate','FDR',
                           'WT_median','MUT_median','WT_target','MUT_target'))
finaltrim=subset(comparison, FDR<0.05, 
                 select = c('gene','Operon','estimate','FDR',
                            'WT_median','MUT_median','WT_target','MUT_target'))
print(nrow(finaltrim))
write.csv(finalall, file=paste0("3.2_",nameWT,"_vs_",nameMUT,"_Linreg.csv"))
write.csv(finaltrim, file=paste0("3.3_",nameWT,"_vs_",nameMUT,"_Linreg_FDR.csv"))
rm(highlight,i,nameMUT,nameWT,suppressor,p,lm_data,lm,l,gene_est,finaltrim,finalall,enhancer,comparison)