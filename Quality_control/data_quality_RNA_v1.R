#RNAseq数据质量控制----
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(ggpubr)
library(patchwork)
library(GenomicRanges)
library(ggh4x)
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Quality_control')
dir()
#FastQC---------
fastqc <- read.table('FastQC_multiqc_data/multiqc_general_stats.txt',sep = '\t',header = T)
fastqc$Sample <- gsub('MGI_PolyA_HCC1395_0313','MGI_PolyA-T_HCC1395',fastqc$Sample) %>% 
  gsub('MGI_PolyA_HCC1395BL_0313','MGI_PolyA-T_HCC139BL',.)
fastqc_2 <- read.table('FastQC_multiqc_data/fastqc_per_sequence_quality_scores_plot.tsv',sep = '\t',header = T,fill = T) %>% t()
colnames(fastqc_2) <- fastqc_2[1,]
fastqc_2 <- fastqc_2[-1,] %>% as.data.frame()
fastqc_2$q30 <- rowSums(fastqc_2[,11:ncol(fastqc_2)],na.rm = TRUE)/rowSums(fastqc_2[,1:ncol(fastqc_2)],na.rm = TRUE)*100

fastqc_2$Sample <- gsub('MGI_PolyA_HCC1395_0313','MGI_PolyA-T_HCC1395',rownames(fastqc_2)) %>% 
  gsub('MGI_PolyA_HCC1395BL_0313','MGI_PolyA-T_HCC139BL',.)

fastqc.all <- merge(fastqc[,c('Sample','FastQC_mqc.generalstats.fastqc.percent_gc','FastQC_mqc.generalstats.fastqc.total_sequences','FastQC_mqc.generalstats.fastqc.percent_duplicates')],
      fastqc_2[,c('q30','Sample')],by='Sample')

fastqc.all$Batch <- paste0(str_split(fastqc.all$Sample,'\\_',simplify = T)[,1],'_',
                           str_split(fastqc.all$Sample,'\\_',simplify = T)[,2],'_',
                           str_split(fastqc.all$Sample,'\\_',simplify = T)[,3])

fastqc.all$Sample_uniq <- fastqc.all$Sample %>% gsub('_R1|_R2|BL','',.)

fastqc.mean.sd <- fastqc.all %>%
  dplyr::group_by(Batch) %>%
  dplyr::summarise(
    mean_bq = mean(q30, na.rm = TRUE),
    sd_bq = sd(q30, na.rm = TRUE),
    mean_gc = mean(FastQC_mqc.generalstats.fastqc.percent_gc, na.rm = TRUE),
    sd_gc = sd(FastQC_mqc.generalstats.fastqc.percent_gc, na.rm = TRUE),
    mean_total = mean(FastQC_mqc.generalstats.fastqc.total_sequences,na.rm = TRUE),
    sd_total =sd(FastQC_mqc.generalstats.fastqc.total_sequences,na.rm = TRUE)
  )

fastqc.uniq <- fastqc.all %>%
  dplyr::group_by(Sample_uniq) %>%
  dplyr::summarise(
    mean_bq = mean(q30, na.rm = TRUE),
    mean_gc = mean(FastQC_mqc.generalstats.fastqc.percent_gc, na.rm = TRUE),
    mean_total = mean(FastQC_mqc.generalstats.fastqc.total_sequences,na.rm = TRUE),
    mean_dup = mean(FastQC_mqc.generalstats.fastqc.percent_duplicates,na.rm = TRUE)
  )


fastqc.mean.sd$mean_total <- fastqc.mean.sd$mean_total/1000000
fastqc.mean.sd$sd_total <- fastqc.mean.sd$sd_total/1000000
fastqc.mean.sd[,-1] <- fastqc.mean.sd[,-1] %>% round(.,2)
write.csv(fastqc.mean.sd,file = 'fastqc.info.csv',row.names = F)

#FastScreen---------
fs <- read.table('Fastqscreen_multiqc_data/multiqc_fastq_screen.txt',sep = '\t',header = T)
colnames(fs)
fs['pollution'] <- fs$EColi.percentage+fs$Vector.percentage+fs$Virus.percentage+fs$Yeast.percentage+fs$Phix.percentage
fs['remnant'] <- fs$rRNA.percentage+fs$Mitoch.percentage
fs$Sample <- gsub('MGI_PolyA_HCC1395_0313','MGI_PolyA-T_HCC1395',fs$Sample) %>% 
  gsub('MGI_PolyA_HCC1395BL_0313','MGI_PolyA-T_HCC1395BL',.)
fs$Batch <- paste0(str_split(fs$Sample,'\\_',simplify = T)[,1],'_',
                   str_split(fs$Sample,'\\_',simplify = T)[,2],'_',
                   str_split(fs$Sample,'\\_',simplify = T)[,3])

fs$Sample_uniq <- fs$Sample %>% gsub('_R1_trimmed_screen|_R2_trimmed_screen|BL','',.)
fs.uniq <- fs %>%
  dplyr::group_by(Sample_uniq) %>%
  dplyr::summarise(
    mean_pl = mean(pollution, na.rm = TRUE),
    mean_rn = mean(remnant, na.rm = TRUE)
  )

fs.mean.sd <- fs %>%
  dplyr::group_by(Batch) %>%
  dplyr::summarise(
    mean_pl = mean(pollution, na.rm = TRUE),
    sd_pl = sd(pollution, na.rm = TRUE),
    mean_rn = mean(remnant, na.rm = TRUE),
    sd_rn = sd(remnant, na.rm = TRUE)
  )
fs.mean.sd[,-1] <- fs.mean.sd[,-1] %>% round(.,2)
write.csv(fs.mean.sd,file = 'fastqscreen.info.csv',row.names = F)

#BamQC & RNAseqQC---------------
bamqc <- read.table('BamQC_multiqc_data/multiqc_general_stats.txt',sep = '\t',header = T)
all_smp <- c(paste0('MGI_PolyA-T_HCC1395',c('_1','_2','_3')),paste0('MGI_PolyA-T_HCC1395BL',c('_1','_2','_3')),
             paste0('MGI_PolyA-T_HCC1395',c('_1','_2','_3'),'_bamqc'),paste0('MGI_PolyA-T_HCC1395BL',c('_1','_2','_3'),'_bamqc'))
names(all_smp) <- c(paste0('MGI_PolyA_HCC1395',c('_1','_2','_3')),paste0('MGI_PolyA_HCC1395BL',c('_1','_2','_3')),
                    paste0('MGI_PolyA_HCC1395',c('_1','_2','_3'),'_bamqc'),paste0('MGI_PolyA_HCC1395BL',c('_1','_2','_3'),'_bamqc'))

bamqc$Sample <- sapply(bamqc$Sample, function(x) {
  if(x %in% names(all_smp)) {
    return(all_smp[x])
  } else {
    return(x)
  }
})

bamqc.sub1 <- subset(bamqc,grepl('_bamqc',bamqc$Sample))


bamqc.sub1$Batch <-  paste0(str_split(bamqc.sub1$Sample,'\\_',simplify = T)[,1],'_',
                            str_split(bamqc.sub1$Sample,'\\_',simplify = T)[,2],'_',
                            str_split(bamqc.sub1$Sample,'\\_',simplify = T)[,3])

bamqc.sub2 <- subset(bamqc,!grepl('_bamqc',bamqc$Sample))
bamqc.sub2$Batch <- paste0(str_split(bamqc.sub2$Sample,'\\_',simplify = T)[,1],'_',
                           str_split(bamqc.sub2$Sample,'\\_',simplify = T)[,2],'_',
                           str_split(bamqc.sub2$Sample,'\\_',simplify = T)[,3])
bamqc.sub2$map_rate <- bamqc.sub2$QualiMap_mqc.generalstats.qualimap.mapped_reads/
  bamqc.sub2$QualiMap_mqc.generalstats.qualimap.total_reads*100

rnaseq <- read.table('RNAseqQC_multiqc_data/qualimap_rnaseq_genome_results.txt',sep='\t',header = T) 
all_smp <- c(paste0('MGI_PolyA-T_HCC1395',c('_1','_2','_3')),paste0('MGI_PolyA-T_HCC1395BL',c('_1','_2','_3')),
             paste0('MGI_PolyA-T_HCC1395',c('_1','_2','_3'),'_bamqc'),paste0('MGI_PolyA-T_HCC1395BL',c('_1','_2','_3'),'_bamqc'))
names(all_smp) <- c(paste0('MGI_PolyA_HCC1395',c('_1','_2','_3')),paste0('MGI_PolyA_HCC1395BL',c('_1','_2','_3')),
                    paste0('MGI_PolyA_HCC1395',c('_1','_2','_3'),'_bamqc'),paste0('MGI_PolyA_HCC1395BL',c('_1','_2','_3'),'_bamqc'))

rnaseq$Sample <- sapply(rnaseq$Sample, function(x) {
  if(x %in% names(all_smp)) {
    return(all_smp[x])
  } else {
    return(x)
  }
})


rnaseq$Batch <- paste0(str_split(rnaseq$Sample,'\\_',simplify = T)[,1],'_',
                       str_split(rnaseq$Sample,'\\_',simplify = T)[,2],'_',
                       str_split(rnaseq$Sample,'\\_',simplify = T)[,3])
rnaseq$intronic_rate <- rnaseq$reads_aligned_intergenic/rnaseq$total_alignments*100

bamqc.sub1$Sample <- gsub('_bamqc','',bamqc.sub1$Sample)
bamqc.all <- merge(bamqc.sub1[,c('Sample','Batch','QualiMap_mqc.generalstats.qualimap.median_insert_size')],
                   bamqc.sub2[,c('Sample','Batch','map_rate')],by='Sample') %>% 
  merge(.,rnaseq[,c('Sample','Batch','intronic_rate','X5_3_bias')],by='Sample')

bamqc.all$Sample_uniq <- gsub('BL','',bamqc.all$Sample)
bamqc.uniq <- bamqc.all %>%
  dplyr::group_by(Sample_uniq) %>%
  dplyr::summarise(
    mean_ins = mean(QualiMap_mqc.generalstats.qualimap.median_insert_size, na.rm = TRUE),
    mean_mr = mean(map_rate, na.rm = TRUE),
    mean_ir = mean(intronic_rate, na.rm = TRUE),
    mean_b53 = mean(X5_3_bias, na.rm = TRUE)
  )


bamqc.mean.sd <- bamqc.all %>%
  dplyr::group_by(Batch) %>%
  dplyr::summarise(
    mean_ins = mean(QualiMap_mqc.generalstats.qualimap.median_insert_size, na.rm = TRUE),
    sd_ins = sd(QualiMap_mqc.generalstats.qualimap.median_insert_size, na.rm = TRUE),
    mean_mr = mean(map_rate, na.rm = TRUE),
    sd_mr = sd(map_rate, na.rm = TRUE),
    mean_ir = mean(intronic_rate, na.rm = TRUE),
    sd_ir = sd(intronic_rate, na.rm = TRUE),
    mean_b53 = mean(X5_3_bias, na.rm = TRUE),
    sd_b53 = sd(X5_3_bias, na.rm = TRUE)
  )
bamqc.mean.sd[,-1] <- bamqc.mean.sd[,-1] %>% round(.,2)
write.csv(bamqc.mean.sd,file = 'BamQC.info.csv',row.names = F)

qc_all_sample <- merge(fastqc.uniq,fs.uniq) %>% merge(.,bamqc.uniq)

write.csv(qc_all_sample,file = 'QC.info.all_sample.RNAseq.csv',row.names = F)

#Duplicates------------
colors_1=c('#D6DEFF' ,'#FFEDC2')
names(colors_1) <- c('PolyA','RiboZero')
  
fastqc <- read.csv("FastQC_multiqc_data/multiqc_fastqc.txt",sep = '\t')
fastqc.stats <- read.csv("FastQC_multiqc_data/multiqc_general_stats.txt",sep='\t')
dup.df = fastqc.stats[,c('Sample','FastQC_mqc.generalstats.fastqc.percent_duplicates')]

fastqc$Platform <- str_split(fastqc$Sample,'\\_',simplify = T)[,1]
fastqc$lib_type <- str_split(fastqc$Sample,'\\_',simplify = T)[,2]
fastqc$Sample_uniq <- gsub("\\_\\d+$", "", fastqc$Sample) %>% gsub('_R1|_R2','',.)

fastqc$Total.Bases <- fastqc$Total.Bases %>% gsub(' Gbp','',.) %>% as.numeric()
fastqc_sum <- fastqc %>% 
  dplyr::group_by(Sample_uniq,Platform,lib_type) %>% 
  dplyr::summarize(Total.Bases=sum(Total.Bases),
                   Total.Sequences=sum(Total.Sequences),
                   median_sequence_length=mean(median_sequence_length),
                   GC=mean(X.GC)
  )

fastqc_sum$Total.Sequences <- fastqc_sum$Total.Sequences/1000000
get_boxplot <- function(df,d_type,y_title,y_lim,yintercept=F,h_line.p=10,h_line.r=20){
  xlim=c('NVG','ARG','SEQ','BMK','MGI')
  df[,'Base_quality'] = df[,d_type]
  df.p <- ggplot(data=df,aes(x=Platform,y=Base_quality,fill=lib_type)) +
    geom_boxplot(cex=0.5,alpha=0.7) +
    geom_signif(comparisons = list(c('Q30','Q40')),test.args =c(exact=FALSE),
                map_signif_level = T,textsize = 5,y = 43)+
    scale_fill_manual(values = colors_1)+
    theme_bw()+
    theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
          axis.text.x =element_text(color = "black",size = 10,face = "bold",angle = 0),
          axis.text.y = element_text(color = "black",size = 10),
          axis.title.x = element_text(color = "black",face = "bold",size = 14),
          axis.title.y = element_text(color = "black",face = "bold",size = 14),
          #axis.ticks.x = element_text(color = "black",face = "bold",size = 14),
          legend.title = element_text(face = "bold",size = 14),
          legend.text = element_text(face = "bold",size = 10,colour = "black"),
          legend.position = "bottom",legend.background = element_blank(),
          panel.grid.minor  = element_line(colour = NA),
          panel.grid.major.x   = element_line(colour = NA),
          panel.background = element_rect(fill="transparent",colour = NA)
    )+
    theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
    theme(strip.text.x = element_text(colour = "black",face = "bold",size = 14)) + 
    theme(strip.background.y = element_rect(fill = "white", colour = "white")) +
    theme(strip.text.y = element_text(colour = "black",face = "bold",size = 14)) + 
    theme(strip.placement = "inside") +
    theme(strip.switch.pad.grid = unit(1, "inch"))+
    labs(x="",y=y_title)+
    expand_limits(y=y_lim,x=xlim)+
    guides(fill=guide_legend(title = "",nrow = 1,byrow = FALSE))
  #+guides(fill=FALSE)
  
  if (yintercept){
    df.p=df.p+ geom_hline(yintercept=h_line.p,color = "red", size = 0.5)+geom_hline(yintercept=h_line.r,color = "blue", size = 0.5)
  }
  df.p;return(df.p)
}
#数据量G------
fastqc_sum <- subset(fastqc_sum,Total.Bases > 7)
total_base.p <- get_boxplot(fastqc_sum,'Total.Bases','Total Bases (Gpb)',c(9,22),yintercept=T);total_base.p

#Reads量------
total_reads.p <- get_boxplot(fastqc_sum,'Total.Sequences','Total Reads (Million)',c(50,150),yintercept=T,
                             h_line.p=60,h_line.r=120);total_reads.p

#GC含量-------
gc.p <- get_boxplot(fastqc_sum,'GC','GC (%)',c(42,52),yintercept=F);gc.p


#碱基质量-------
bq.df <- read.csv('FastQC_multiqc_data/fastqc_per_base_sequence_quality_plot.tsv',sep = '\t')

bq.df.long <- tidyr::pivot_longer(bq.df, 
             cols = colnames(bq.df)[-1], # 列1和列2是标识符列，不需要被转换
             names_to = "Sample", # 新的列名
             values_to = "Base_quality") 

bq.df.long$Platform <- str_split(bq.df.long$Sample,'\\_',simplify = T)[,1]
bq.df.long$lib_type <- str_split(bq.df.long$Sample,'\\_',simplify = T)[,2]
bq.df.long$Sample_uniq <- gsub("\\_\\d+$", "", bq.df.long$Sample)

bq.p <- get_boxplot(bq.df.long,'Base_quality','Base Quality',c(30,45),yintercept=F);bq.p

fastq_qc <- total_base.p+total_reads.p+gc.p+bq.p+plot_layout(guides = "collect",ncol = 2)&theme(legend.position='bottom')
ggsave("RNA.fastq_qc.png",fastq_qc,width=8.15, height=5.7*1.1,dpi = 300)





