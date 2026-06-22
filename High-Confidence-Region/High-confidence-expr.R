setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/High-Confidence-Region')
#构建高置信参考数据集----
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
#1.高置信基因集-----
my_theme <- theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 0),
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.x = element_text(color = "black",face = "bold",size = 16),
        axis.title.y = element_text(color = "black",face = "bold",size = 16),
        #axis.ticks.x = element_text(color = "black",face = "bold",size = 14),
        legend.title = element_text(face = "bold",size = 16),
        legend.text = element_text(face = "bold",size = 16,colour = "black"),
        legend.position = "bottom",legend.background = element_blank(),
        panel.grid.minor  = element_line(colour = NA),
        panel.grid.major.x   = element_line(colour = NA),
        panel.background = element_rect(fill="transparent",colour = NA)
  )+
  theme(strip.background.x = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.x = element_text(colour = "black",face = "bold",size = 16)) + 
  theme(strip.background.y = element_rect(fill = "white", colour = "white")) +
  theme(strip.text.y = element_text(colour = "black",face = "bold",size = 16)) + 
  theme(strip.placement = "inside") +
  theme(strip.switch.pad.grid = unit(1, "inch"))

DEGanalysis <- function(exprMat, group){
  library(edgeR)
  library(limma)
  
  dge <- DGEList(counts = exprMat)
  design <- model.matrix(~group)
  
  keep <- filterByExpr(dge, design,min.count = 0)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, design, plot=F)
  fit <- lmFit(v, design)
  
  fit <- eBayes(fit)
  result <- topTable(fit, coef=ncol(design), sort.by = 'logFC', number = Inf)
  result$gene = rownames(result)
  result$groupA =  levels(group)[1]
  result$groupB =  levels(group)[2]
  return(as.data.table(result))
}
get_deg_ref <- function(raw_count){
  df_count <- read.csv(raw_count)
  name=str_split(raw_count,'\\.',simplify = T)[,1]
  print(paste0(name,' Processing=========================='))
  if (name=='Stringtie'){
    df_count$Ensembl <- str_split(df_count$gene_id,'\\|',simplify = T)[,1]
    rownames(df_count) <- df_count$Ensembl
    df_count <- df_count[,!grepl('gene_id|Ensembl',colnames(df_count))]
  } else {
    df_count <- subset(df_count,Ensembl!='')
    rownames(df_count) <- df_count$Ensembl
    df_count <- df_count[,!grepl('RefSeq|Ensembl',colnames(df_count))]
  }

  batch_lst = c(ARG_PolyA='B1',
                BMK_RiboZero='B2',
                MGI_RiboZero='B3',
                NVG_PolyA='B4',
                SEQ_PolyA='B5',
                SEQ_RiboZero='B6',
                MGI_PolyA.T='B7',
                MGI_PolyA.M='B8')
  lab_lst = c(ARG='L1',BMK='L2',MGI='L3',NVG='L4',SEQ='L5')
  colnames(df_count) <- gsub('_0313','',colnames(df_count) )
  mgi_name <- read.csv('MGI_PolyA.M_name.csv',sep = ',')
  colnames(df_count) <- sapply(colnames(df_count),function(x){
    if (x %in% mgi_name$raw_name){
      return(subset(mgi_name,raw_name==x)$name_new)
    }else{
      return(x)
    }
  })
  mt_df <- data.frame(library=colnames(df_count),
                      sample=str_split(colnames(df_count),'\\_',simplify = T)[,3],
                      protocol=str_split(colnames(df_count),'\\_',simplify = T)[,2],
                      platform=str_split(colnames(df_count),'\\_',simplify = T)[,1],
                      lab=lapply(colnames(df_count),function(x){
                        aim = str_split(colnames(df_count),'\\_',simplify = T)[,1]
                        return(lab_lst[aim])}) %>% unlist() %>% as.character(),
                      batch=lapply(colnames(df_count),function(x){
                        aim = paste0(str_split(colnames(df_count),'\\_',simplify = T)[,1],'_',
                                     str_split(colnames(df_count),'\\_',simplify = T)[,2])
                        return(batch_lst[aim])}) %>% unlist() %>% as.character(),
                      rep=str_split(colnames(df_count),'\\_',simplify = T)[,4],
                      sample_rep=paste0(str_split(colnames(df_count),'\\_',simplify = T)[,3],'_',
                                        str_split(colnames(df_count),'\\_',simplify = T)[,4]))
  mt_df$protocol <- mt_df$protocol %>% gsub('.T|.M','',.)
  write.csv(mt_df,file = 'hcc1395_metadata.csv')
  write.csv(df_count,file = paste0(name,'.multibatch_count.csv'))
  #1. 计算可定量基因
  system(paste0('Rscript /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Reference_standard_protocol-main/scripts/char_detect-hcc1395.R',
                ' -i ',paste0(name,'.multibatch_count.csv'),' -m ','hcc1395_metadata.csv',' -o ','./'))
  
  #2. 取出可检测的基因表达矩阵
  dete_genes.df <- read.csv(sub('.csv','.detect_genelist.csv',paste0(name,'.multibatch_count.csv')))
  df.freq <-table(dete_genes.df$gene) %>% as.data.frame()
  df_count.dete <- df_count[subset(df.freq,Freq>1)$Var1,]
  df_count.dete[is.na(df_count.dete)]<- 0
  write.csv(df_count.dete,file = sub('.csv','.detectable.csv',paste0(name,'.multibatch_count.csv')))
  
  #3. 计算差异表达基因
  system(paste0('Rscript /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Reference_standard_protocol-main/scripts/deg_limma-hcc1395.R',
                ' -i ',sub('.csv','.detectable.csv',paste0(name,'.multibatch_count.csv')),' -m ','hcc1395_metadata.csv',' -o ','./'))
  
  #4. 去除PolyA和RiboZero两种不同方案之间显示显着差异的基因
  usample<-c("RiboZero","PolyA")
  sample_combn<-data.frame(
    sampleA=rep("RiboZero",3),
    sampleB= c("PolyA"))
  
  colnames(sample_combn)<-c("sampleA","sampleB")
  sample_combn$compare<-paste(sample_combn$sampleB,"/",sample_combn$sampleA,sep="")
  DEG.cal.P2R<-c()
  ubatch <- unique(mt_df$sample)
  for ( i in 1:length(ubatch)){
    meta_dat<-mt_df[mt_df$sample==ubatch[i],]
    expr<-df_count.dete[,match(meta_dat$library,colnames(df_count.dete))]
    
    j=1
    lst.library.sampleA=meta_dat$library[meta_dat$protocol==sample_combn$sampleA[j]]
    lst.library.sampleB=meta_dat$library[meta_dat$protocol==sample_combn$sampleB[j]]
    lst.library.forDEG <- c(lst.library.sampleA, lst.library.sampleB)
    
    y.DEG <- DEGanalysis(expr[,lst.library.forDEG],
                         factor(as.character(meta_dat$protocol[match(lst.library.forDEG,meta_dat$library)])))
    
    y.DEG$batch<-ubatch[i]
    y.DEG<-data.frame(y.DEG)
    y.DEG$compare<-paste(y.DEG$groupB,"/",y.DEG$groupA,sep="")
    
    y.DEG_rev<-data.frame(
      logFC= -y.DEG$logFC,
      AveExpr=y.DEG$AveExpr,
      t=y.DEG$t,
      P.Value=y.DEG$P.Value,
      adj.P.Val=y.DEG$adj.P.Val,
      B=y.DEG$B,
      gene=y.DEG$gene,
      groupA=y.DEG$groupB,
      groupB=y.DEG$groupA,
      # When limma does a variance analysis, the name of the numerator group of the FC must be sorted after the character order of the name of the denominator group.
      batch=y.DEG$batch,
      compare=rep(sample_combn$compare[j],nrow(y.DEG))
    )
    
    DEG.cal.P2R<-rbind(DEG.cal.P2R,y.DEG_rev)
    
    print(i)
  }
  
  DEG.cal.P2R<-data.frame(DEG.cal.P2R)
  DEG.cal.P2R$Type<-"non-DEG"
  DEG.cal.P2R$Type[intersect(which(DEG.cal.P2R$P.Value<0.05),which(DEG.cal.P2R$logFC>=1))]<-"up-regulate"
  DEG.cal.P2R$Type[intersect(which(DEG.cal.P2R$P.Value<0.05),which(DEG.cal.P2R$logFC<=(-1)))]<-"down-regulate"
  DEG.cal.P2R$Type %>% table()
  false_genes <- subset(DEG.cal.P2R,Type!='non-DEG')$gene %>% unique()
  
  deg_limma <- read.csv(sub('.csv','.detectable.csv',paste0(name,'.multibatch_count.csv')) %>% sub('.csv','.DEG_limma.csv',.))
  deg_limma.sub <- deg_limma[!deg_limma$gene %in% false_genes,]
  write.csv(deg_limma.sub,file = sub('.csv','.detectable.csv',
                                     paste0(name,'.multibatch_count.csv')) %>% 
              sub('.csv','.DEG_limma.csv',.) %>% sub('.csv','.rmFalseGenes.csv',.))
  #5. 进行Shapiro-Wilk检验
  batch_include <- data.frame(batch=deg_limma.sub$batch %>% unique())
  write.csv(batch_include,file = 'batch_included.csv',row.names = F)
  
  system(paste0('Rscript /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Reference_standard_protocol-main/scripts/normality-hcc1395.R',
                ' -i ',sub('.csv','.detectable.csv',
                           paste0(name,'.multibatch_count.csv')) %>% 
                  sub('.csv','.DEG_limma.csv',.) %>% sub('.csv','.rmFalseGenes.csv',.),' -p ','batch_included.csv',' -o ','./'))
  #6. 去除Shapiro-Wilk检验不合格的基因(这步不用做)
  #norm_df <- read.csv(sub('.csv','.detectable.csv',
  #                        paste0(name,'.multibatch_count.csv')) %>% 
  #                      sub('.csv','.DEG_limma.csv',.) %>% sub('.csv','.rmFalseGenes.csv',.) %>% sub('.csv','.normality_genelist.csv',.))
  #false_genes_2 <- subset(norm_df,log2fc_p>0.05)$gene %>% unique()
  #
  #deg_limma.sub.2 <- deg_limma.sub[!deg_limma.sub$gene %in% false_genes_2,]
  #write.csv(deg_limma.sub.2,file =   sub('.csv','.detectable.csv',
  #                                       paste0(name,'.multibatch_count.csv')) %>% 
  #            sub('.csv','.DEG_limma.csv',.) %>% sub('.csv','.rmFalseGenes.csv',.) %>% sub('.csv','.rmUnNorm.csv',.))
  
  #7. 生成高置信参考数据集------
  system(paste0('Rscript /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Reference_standard_protocol-main/scripts/char_ratio-hcc1395.R',
                ' -i ',sub('.csv','.detectable.csv',paste0(name,'.multibatch_count.csv')),
                ' -m ','hcc1395_metadata.csv',
                ' -g ',sub('.csv','.detect_genelist.csv',paste0(name,'.multibatch_count.csv')),
                ' -o ','./',
                ' -d ',sub('.csv','.detectable.csv',
                           paste0(name,'.multibatch_count.csv')) %>% 
                  sub('.csv','.DEG_limma.csv',.) %>% sub('.csv','.rmFalseGenes.csv',.)))
  
  #8. 识别高置信差异表达基因集----
  system(paste0('Rscript /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Reference_standard_protocol-main/scripts/char_refDEG-hcc1395.R',
                ' -r ',sub('.csv','.detectable.csv',paste0(name,'.multibatch_count.csv')) %>% sub('.csv','.ref_expr.csv',.),
                ' -o ','./',
                ' -d ',sub('.csv','.detectable.csv',
                           paste0(name,'.multibatch_count.csv')) %>% 
                  sub('.csv','.DEG_limma.csv',.) %>% sub('.csv','.rmFalseGenes.csv',.),
                ' -m ','hcc1395_metadata.csv'))
  
  print(paste0(name,' Finished=========================='))
  return(paste0(name,' Finished=========================='))
}

get_deg_ref("FeatureCounts.gene_quant.csv")

get_deg_ref("Stringtie.gene_quant.csv")

get_deg_ref("Salmon.gene_quant.csv")

get_deg_ref("HTseq.gene_quant.csv")


#合并四种工具检测得到的高置信参考数据集----
##要求至少两种工具支持，且表达类型一致，则保留，其他值均保留均值
dir(pattern = '*.RefData_DEGs.csv') 
fc_ref <- read.csv("FeatureCounts.multibatch_count.detectable.ref_expr.RefData_DEGs.csv")
fc_ref$tool <- 'FeatureCounts'

fc_ref[,c('n_up','n_non','n_down','final')]
fc_ref <- fc_ref %>%
  mutate(support_n = case_when(
    final == 'non-DEG' ~ n_non,
    final == 'up-regulate' ~ n_up,
    final == 'down-regulate' ~ n_down,
    TRUE ~ NA_real_
  ))

st_ref <- read.csv("Stringtie.multibatch_count.detectable.ref_expr.RefData_DEGs.csv")
st_ref$tool <- 'Stringtie'
st_ref <- st_ref %>%
  mutate(support_n = case_when(
    final == 'non-DEG' ~ n_non,
    final == 'up-regulate' ~ n_up,
    final == 'down-regulate' ~ n_down,
    TRUE ~ NA_real_
  ))

sm_ref <- read.csv("Salmon.multibatch_count.detectable.ref_expr.RefData_DEGs.csv")
sm_ref$tool <- 'Salmon'
sm_ref <- sm_ref %>%
  mutate(support_n = case_when(
    final == 'non-DEG' ~ n_non,
    final == 'up-regulate' ~ n_up,
    final == 'down-regulate' ~ n_down,
    TRUE ~ NA_real_
  ))

ht_ref <- read.csv("HTseq.multibatch_count.detectable.ref_expr.RefData_DEGs.csv")
ht_ref$tool <- 'HTseq'
ht_ref <- ht_ref %>%
  mutate(support_n = case_when(
    final == 'non-DEG' ~ n_non,
    final == 'up-regulate' ~ n_up,
    final == 'down-regulate' ~ n_down,
    TRUE ~ NA_real_
  ))

all_ref <- rbind(fc_ref,st_ref) %>% rbind(sm_ref,.) %>% rbind(ht_ref,.)

all_ref.sub <- all_ref[,c('gene_compare','fc','medianp','final','tool','support_n')]

all_ref.sub <- subset(all_ref.sub,grepl('^EN',all_ref.sub$gene_compare))

all_ref_cmb <- all_ref.sub %>% dplyr::group_by(gene_compare) %>%
  dplyr::summarise(fc_mean = mean(fc),
                   p_mean = mean(medianp),
                   freq = n(),
                   support_n=mean(support_n))

aim_genes <- subset(all_ref_cmb,freq >=2)$gene_compare

all_ref_cmb.new <- subset(all_ref_cmb,all_ref_cmb$gene_compare %in% aim_genes)

DEG_ref_cal<-data.frame(
  N_up=tapply(all_ref.sub$final,as.factor(all_ref.sub$gene_compare),function(x){length(which(x=="up-regulate"))}),
  N_non=tapply(all_ref.sub$final,as.factor(all_ref.sub$gene_compare),function(x){length(which(x=="non-DEG"))}),
  N_down=tapply(all_ref.sub$final,as.factor(all_ref.sub$gene_compare),function(x){length(which(x=="down-regulate"))}))

DEG_ref_cal$final<-"non-DEG"
DEG_ref_cal$final[intersect(which(DEG_ref_cal$N_up>=1),which(DEG_ref_cal$N_down>=1))]<-"conflicting"

DEG_ref_cal$final[intersect(which(DEG_ref_cal$N_up>=2),which(DEG_ref_cal$N_down==0))]<-"up-regulate"
DEG_ref_cal$final[intersect(which(DEG_ref_cal$N_down>=2),which(DEG_ref_cal$N_up==0))]<-"down-regulate"

DEG_ref_cal$gene_compare <- rownames(DEG_ref_cal)

all_ref.final <- dplyr::left_join(all_ref_cmb.new,DEG_ref_cal,by='gene_compare')
all_ref.final$`Confidence(%)` <- all_ref.final$support_n/8*100
all_ref.final <- all_ref.final[,c('gene_compare','fc_mean','p_mean','support_n','final','Confidence(%)')]

write.csv(all_ref.final,file = 'All.RefData_DEGs.csv',row.names = F)





#统计差异表达基因的数量--------------
#1. 注释到的基因--

for (i in dir(pattern = '*quant.csv')){
  name = sub('.gene_quant.csv','',i)
  print(name)
  df <- read.csv(i)
  print(nrow(df))
}

#2. 可检测的基因--
for (i in dir(pattern = '*.detect_genelist.csv')){
  name = str_split(i,'\\.',simplify = T)[,1]
  print(name)
  df <- read.csv(i)
  for (a in df$sample %>% unique()){
    print(a)
    print(subset(df,sample==a)$gene %>% length())
  }
}


#3.两个样本中都可检测的基因--
for (i in dir(pattern = '*.DEG_limma.csv')){
  name = str_split(i,'\\.',simplify = T)[,1]
  print(name)
  df <- read.csv(i)
  print(df$gene %>% unique() %>% length())
  
}

#4. 去除PolyA和RiboZero两种不同方案之间显示显着差异的基因--

for (i in dir(pattern = '*.rmFalseGenes.csv')){
  name = str_split(i,'\\.',simplify = T)[,1]
  print(name)
  df <- read.csv(i)
  print(df$gene %>% unique() %>% length())
  
}

#5. 差异表达的基因---

for (i in dir(pattern = '*.RefData_DEGs.csv')){
  name = str_split(i,'\\.',simplify = T)[,1]
    print(name)
    df <- read.csv(i)
    print(subset(df,final!='conflicting')$gene %>% unique() %>% length())

}







#统计基因表达水平的正太分布情况---
for (i in dir(pattern = '*.DEG_limma.csv')){
  system(paste0('Rscript /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Reference_standard_protocol-main/scripts/normality-hcc1395.R',
                ' -i ',i,' -p ','batch_included.csv',' -o ','./'))
  
  print(i)
  df <- read.csv(sub('.csv','.normality_genelist.csv',i))
  print(sum(df$log2fc_p > 0.05)/nrow(df)*100)
}

#计算定值引入的不确定度----
##1. 首先合并所有的DEG_limma.csv文件
df_limma <- read.csv(dir(pattern = '*.DEG_limma.csv')[1])
for (i in dir(pattern = '*.DEG_limma.csv')[-1]){
  print(i)
  df <- read.csv(i)
  df_limma <- rbind(df_limma,df)
}


write.csv(df_limma,file = 'All.DEG_limma.csv')

system(paste0('Rscript /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Reference_standard_protocol-main/scripts/uchar-hcc1395.R',
              ' -i All.DEG_limma.csv',' -p ','batch_included.csv',' -o ','./',' -r All.RefData_DEGs.csv'))

#uchar < 15 
ref.uchar <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/High-Confidence-Region/All.DEG_limma.Refdata_uchar.csv')
ref.uchar$uchar[is.na(ref.uchar$uchar)]=0
data.frame(uchar=c('0~15%','15~30%','30~40%'),
                       num=c(sum(ref.uchar$uchar < 15,na.rm = TRUE),
                             sum(ref.uchar$uchar >= 15 & ref.uchar$uchar < 30,na.rm = TRUE),
                             sum(ref.uchar$uchar >= 30,na.rm = TRUE)),
                       perc=c(sum(ref.uchar$uchar < 15,na.rm = TRUE)/nrow(ref.uchar),
                              sum(ref.uchar$uchar >= 15 & ref.uchar$uchar < 30,na.rm = TRUE)/nrow(ref.uchar),
                              sum(ref.uchar$uchar >= 30,na.rm = TRUE)/nrow(ref.uchar)))

uchar_df <- data.frame(uchar=c('0~20%','20%~'),
                       num=c(sum(ref.uchar$uchar < 20,na.rm = TRUE),
                             sum(ref.uchar$uchar >= 20,na.rm = TRUE)),
                       perc=c(sum(ref.uchar$uchar < 20,na.rm = TRUE)/nrow(ref.uchar),
                              sum(ref.uchar$uchar >= 20,na.rm = TRUE)/nrow(ref.uchar)))

#统计稳定性引入的不确定度----
get_u_1 <- function(u_lst,type){
  uhom_df <- data.frame(u_logFC=c('0~15%','15~30%','30~45%','45~60%','60%~'),
                        num=c(sum(u_lst < 0.15,na.rm = TRUE),
                              sum(u_lst >= 0.15 & u_lst < 0.3,na.rm = TRUE),
                              sum(u_lst >= 0.3 & u_lst < 0.45,na.rm = TRUE),
                              sum(u_lst >= 0.45 & u_lst < 0.6,na.rm = TRUE),
                              sum(u_lst >= 0.6,na.rm = TRUE)),
                        perc=c(sum(u_lst < 0.15,na.rm = TRUE)/length(u_lst),
                               sum(u_lst >= 0.15 & u_lst < 0.3,na.rm = TRUE)/length(u_lst),
                               sum(u_lst >= 0.3 & u_lst < 0.45,na.rm = TRUE)/length(u_lst),
                               sum(u_lst >= 0.45 & u_lst< 0.6,na.rm = TRUE)/length(u_lst),
                               sum(u_lst >= 0.6,na.rm = TRUE)/length(u_lst)))
  uhom_df$type <- type
  return(uhom_df)
}

get_u <- function(u_lst,type){
  uhom_df <- data.frame(u_logFC=c('0~20%','20%~'),
                        num=c(sum(u_lst < 0.2,na.rm = TRUE),
                              sum(u_lst >= 0.2,na.rm = TRUE)),
                        perc=c(sum(u_lst < 0.2,na.rm = TRUE)/length(u_lst),
                               sum(u_lst >= 0.2,na.rm = TRUE)/length(u_lst)
                        ))
  uhom_df$type <- type
  return(uhom_df)
}

ref_uncer <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/Refdata_add_uncertainty.csv')

#Long stability at -80°Ç
#ls_file='/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/stability_lm_detail_summary.csv'
#ls_df <- read.csv(ls_file,header = T)
get_u(ref_uncer$ulst,'ulst')

#Short stability at 4°C
ss_file='/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/DW4_stability/stability_dw4_detail_summary.csv'
ss_df <- read.csv(ss_file,header = T)
get_u(ss_df$usw4,'usw4')


#KD stability at -20°C
kd_file = '/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/KD_stability/stability_kd_detail_summary.csv'
kd_df <- read.csv(kd_file,header = T)
get_u(kd_df$ukd,'ukd')

#Short stability at 25°C
ss25_file = '/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/DW25_stability/stability_dw25_detail_summary.csv'
ss25_df <- read.csv(ss25_file,header = T)
get_u(ss25_df$usw25,'usw25')

#Short stability at CO2
ssco2_file = '/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/DWCO2_stability/stability_dwCO2_detail_summary.csv'
ssco2_df <- read.csv(ssco2_file,header = T)
get_u(ssco2_df$uco2,'uco2')

#Stability all
#ust_all <- merge(ref_uncer[ref_uncer$long_stability=='FALSE',c('compare','ulst')],
#                 ref_uncer[ref_uncer$dw4_stability=='FALSE',c('compare','usw4')]) %>% merge(.,ref_uncer[ref_uncer$dwCO2_stability=='FALSE',c('compare','uco2')])
ref_uncer_stab <- subset(ref_uncer,homoge=='FALSE' & long_stability == 'FALSE' & dw4_stability=='FALSE' & dwCO2_stability=='FALSE')
ref_uncer_stab$usw4 <- as.numeric(ref_uncer_stab$usw4)
ref_uncer_stab$usco2 <- as.numeric(ref_uncer_stab$usco2)
ref_uncer_stab$uts <- sqrt(ref_uncer_stab$ulst^2+ref_uncer_stab$usw4^2+ref_uncer_stab$usco2^2)
get_u(ref_uncer_stab$uts,'uts')


#combined uncertainty
table(ref_uncer_stab$uhom_logFC < ref_uncer_stab$uchar/3)
table(ref_uncer_stab$uts < ref_uncer_stab$uchar/3)
ref_uncer_stab$uc <- sqrt(ref_uncer_stab$uchar^2 + ref_uncer_stab$uhom_logFC^2 + ref_uncer_stab$uts^2)
get_u(ref_uncer_stab$uc,'uc')

ref_uncer_stab$U <- 2*ref_uncer_stab$uc
get_u(ref_uncer_stab$U,'U')

std_u <- rbind(get_u(ref_uncer$ulst,'ulst'),get_u(ss_df$usw4,'usw4')) %>% rbind(.,get_u(kd_df$ukd,'ukd')) %>% 
  rbind(.,get_u(ss25_df$usw25,'usw25')) %>% rbind(.,get_u(ssco2_df$uco2,'uco2')) %>% 
  rbind(.,get_u(ref_uncer_stab$uts,'uts')) %>% rbind(.,get_u(ref_uncer_stab$uc,'uc')) %>% rbind(.,get_u(ref_uncer_stab$U,'U'))


write.csv(std_u,file = '/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/All_uncertainty.csv',row.names = F)
write.csv(ref_uncer_stab,file = '/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/Refdata_add_uncertainty.pass.csv',row.names = F)

#给每个基因添加Confidence
ref_uncer_stab <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/Refdata_add_uncertainty.pass.csv')
ref_uncer_stab <- merge(ref_uncer_stab,all_ref.final[,c('gene_compare','Confidence(%)')])
Confidence_df <- data.frame(Confidence=c('0~90%','90%~'),
                            num=c(sum(ref_uncer_stab$`Confidence(%)` < 90,na.rm = TRUE),
                                  sum(ref_uncer_stab$`Confidence(%)` >= 90,na.rm = TRUE)),
                            perc=c(sum(ref_uncer_stab$`Confidence(%)` < 90,na.rm = TRUE)/nrow(ref_uncer_stab),
                                   sum(ref_uncer_stab$`Confidence(%)` >= 90,na.rm = TRUE)/nrow(ref_uncer_stab)))

#=------------
#与Quartet的标准参考集进行比较
hcc.ref <- read.csv('All.RefData_DEGs.csv',header = TRUE)
hcc.ref$gene <- str_split(hcc.ref$gene_compare,'\\ ',simplify = T)[,1]

gene_type <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/High-Confidence-Region/Homo_sapiens.GRCh38.109.gene_biotype.txt',
                      sep = '\t',header = F)
gene_type$gene <- str_split(gene_type$V9,pattern = "\\;",simplify = T)[,1]
gene_type$gene <- gsub('gene_id ','',gene_type$gene)

gene_type$gene_biotype <- str_extract(gene_type$V9, "gene_biotype\\s+([^;]+)") %>%
  str_remove("gene_biotype\\s+") %>%
  str_trim()

gene_type <- gene_type[,c('gene','gene_biotype')] %>% unique()

hcc.ref <- merge(hcc.ref,gene_type) 
#基因表达类型的分类---
final_counts <- hcc.ref %>%  # 请将 your_data 替换为你的数据框名称
  count(final)

# 绘制条形图
p_DEG <- ggplot(final_counts, aes(x = final, y = n, fill = final)) +
  geom_col(color = "black", width = 0.7) + # 添加黑色边框并调整宽度
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) + # 在柱子上方添加计数标签
  scale_fill_brewer(palette = "Set2") + # 使用 Set2 调色板
  labs(title = "Distribution of DEG Categories",
       x = "final Category",
       y = "Count") +
  my_theme+theme(axis.text.x = element_text(angle = 30, hjust = 1),
                 legend.position = "none")
p_DEG

#计算基因类型的分类----
gene_biotype_counts <- hcc.ref %>% 
  count(gene_biotype)

threshold <- 10
gene_biotype_counts <- gene_biotype_counts %>% 
  mutate(
    grouped_biotype = ifelse(n < threshold, "Other", as.character(gene_biotype))
  ) %>%
  group_by(grouped_biotype) %>%
  summarise(n = sum(n))

p_Biotype <- ggplot(gene_biotype_counts, aes(x = grouped_biotype, y = n, fill = grouped_biotype)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.5, size = 3.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Distribution of Gene Biotype",
       x = "Gene Biotype",
       y = "Count") +
  my_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

p_Biotype

#计算可信度分布----
p_confidence <- ggplot(hcc.ref, aes(x = Confidence...)) +
  geom_histogram(aes(y = ..density..), # 绘制密度直方图
                 bins = 20, # 设置分箱数
                 fill = "#69b3a2", 
                 color = "black", 
                 alpha = 0.7) +
  geom_density(color = "#e9a3c9", linewidth = 1) + # 叠加密度曲线
  labs(title = "Histogram of Confidence",
       x = "Confidence",
       y = "Density") +
  my_theme+
  theme(plot.title = element_text(hjust = 0.5))

p_confidence

library(patchwork)
p_HC_all <- p_DEG|p_Biotype|p_confidence

ggsave('HCC_ref_character.tiff',plot = p_HC_all,dpi = 300,width =8.15,height =5.7*2)
ggsave('HCC_ref_character.pdf',plot = p_HC_all,dpi = 300,width =8.15,height =5.7)


hcc.ref$gene <- str_split(hcc.ref$gene_compare,'\\ ',simplify = T)[,1]
hcc.ref$meanlogFC <- log2(hcc.ref$fc_mean+0.001)
quartet.ref <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/quartet-rseqc-report/exp2qcdt/data/ref_data_fc_value.csv',
                        header = TRUE)

#ADAR
subset(quartet.ref,gene=='ENSG00000197381')
subset(quartet.ref,gene=='ENSG00000185736')
subset(quartet.ref,gene=='ENSG00000160710')
subset(quartet.ref,gene=='ENSG00000205696')

#参考数据集的交集和差集----
library(ggVennDiagram)
library(ggplot2)
library("animation")
library(tidyr)
library(ggforce)
library(tidyverse)
hcc.gene <- unique(hcc.ref$gene)
quartet.gene <- unique(quartet.ref$gene)
combine_df <- data.frame(ID = c(hcc.gene,quartet.gene),
                         Lab = c(rep('HCC1395',length(hcc.gene)),
                                 rep('Quartet',length(quartet.gene))))

combine_df <- na.omit(combine_df)

venn_plot <- function(df,p_name,p_path){
  list_df = list()
  for (i in df[,2]){
    list_df[[i]] = df[df[,2]==i,1]
  }

  if (T){
    vp = ggVennDiagram(list_df,
                       label = "count",
                       edge_lty = "dashed",
                       edge_size = 0.5,
                       set_size = 5,
                       label_txtWidth = 1,
                       label_size = 5,
                       scale_size = 0.5)+scale_fill_distiller(palette = "RdBu")+
      theme_void()+labs(title = paste0("Reference datasets ","(",p_name,")"))+
      theme(legend.position = 'right',
            plot.margin=unit(c(1, 1, 1, 1),'cm'),
            plot.title = element_text(hjust = 0.5),
            text = element_text(size = 15),)+scale_x_continuous(expand = expansion(mult = .2))
  }

  ggsave(paste0(p_path,'/',p_name,'_','venn_plot.png'),plot = vp,dpi = 300,height =8.15*0.5,width =5.7)
  ggsave(paste0(p_path,'/',p_name,'_','venn_plot.pdf'),plot = vp,dpi = 300,height =8.15*0.5,width =5.7)
  
}

venn_plot(combine_df,'HCC1395.vs.Quartet',"HCC1395.vs.Quartet/")

common_gene.quartet.vs.hcc <- intersect(hcc.gene,quartet.gene)

quartet.ref_common <- subset(quartet.ref,quartet.ref$gene %in% common_gene.quartet.vs.hcc)[,c('gene','compare','meanlogFC')] %>% 
  tidyr::pivot_wider(
    names_from = compare,
    values_from = meanlogFC,
    values_fill = list(meanlogFC = NA)
    
  )

common_gene.quartet.vs.hcc.df <- merge(quartet.ref_common,
      subset(hcc.ref,hcc.ref$gene %in% common_gene.quartet.vs.hcc)[,c('gene','final','fc_mean')])
common_gene.quartet.vs.hcc.df$meanlogFC = log2(common_gene.quartet.vs.hcc.df$fc_mean)

# 使用ggplot2创建更精美的图形
library(patchwork)

# 创建绘图函数（ggplot2版本）
create_gg_corr_plot <- function(data, y_var, x_var, y_lab) {
  # 过滤NA值
  plot_data <- data %>% 
    select(all_of(c(x_var, y_var))) %>% 
    na.omit()
  
  if (nrow(plot_data) > 1) {
    # 计算统计量
    cor_test <- cor.test(plot_data[[x_var]], plot_data[[y_var]])
    cor_coef <- round(cor_test$estimate, 3)
    p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                      round(cor_test$p.value, 3))
    
    # 创建图形
    ggplot(plot_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
      geom_point(color = "steelblue", alpha = 0.7, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "red", alpha = 0.2) +
      labs(y = y_lab, x = "HCC1395/BL",
           
           subtitle = paste("r =", cor_coef, ", p =", p_value)) +
      my_theme
  }
}

# 创建三个图形
p1 <- create_gg_corr_plot(common_gene.quartet.vs.hcc.df, "D5/D6", "meanlogFC", "D5/D6")
p2 <- create_gg_corr_plot(common_gene.quartet.vs.hcc.df, "F7/D6", "meanlogFC", "F7/D6") 
p3 <- create_gg_corr_plot(common_gene.quartet.vs.hcc.df, "M8/D6", "meanlogFC", "M8/D6")

# 合并图形
p_HCC_quartet.cor <- (p1 / p2 / p3) 

ggsave(
  filename = "HCC1395.vs.Quartet/HCC_quartet.cor.pdf",
  plot = p_HCC_quartet.cor,
  width = 3,
  height = 6,
  dpi = 300
)


write.csv(common_gene.quartet.vs.hcc.df,file = "HCC1395.vs.Quartet/ccommon_gene.quartet.vs.hcc.df.csv")

consistency_df <- data.frame(compare=1,
                             consistency=1,
                             binom_p=0)

for (i in c("D5/D6", "F7/D6",'M8/D6')){
  # 首先处理NA值，只保留两列都有值的行
  df <- common_gene.quartet.vs.hcc.df
  df_filtered <- df[,c(i, "final")]
  df_filtered <- na.omit(df_filtered)
  colnames(df_filtered) <- c('A','B')
  
  # 计算一致性：两列值相同的行数
  consistent <- sum(df_filtered$A == df_filtered$B)
  
  # 计算总比较行数
  total <- nrow(df_filtered)
  
  # 计算一致性概率
  consistency_prob <- consistent / total
  
  # 进行二项检验
  # 原假设：一致性概率为0.5（随机猜测水平）
  binom_result <- binom.test(consistent, total, p = 0.5, alternative = "two.sided")
  
  compare=paste0(i,'.vs.','HCC1395')
  binom_p=binom_result$p.value

  df= data.frame(compare=compare,
                  consistency=consistency_prob,
                  binom_p=binom_p)
  consistency_df = rbind(consistency_df,df)
}


#与MAQC sample A/B比较----
maqc.ref <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/ElementVSIllumina/RNAseq/MAQC2_sampleA_sampleB/MAQC_taq_reference.csv',
                     sep = ',')
gene_ids <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Ensembl2RefSeq.txt',sep = '\t',header = F)
colnames(gene_ids)=c('gene','Gene.Name')
maqc.ref <- merge(maqc.ref,gene_ids,by='Gene.Name')

maqc.gene <- unique(maqc.ref$gene)
combine_df <- data.frame(ID = c(hcc.gene,maqc.gene),
                         Lab = c(rep('HCC1395',length(hcc.gene)),
                                 rep('MAQC',length(maqc.gene))))

combine_df <- na.omit(combine_df)

venn_plot(combine_df,'HCC1395.vs.MAQC',"HCC1395.vs.Quartet/")

common_gene.maqc.vs.hcc <- intersect(maqc.ref$gene,hcc.ref$gene)
common_gene.maqc.vs.hcc.df <- merge(maqc.ref[,c('gene','DEG')],
                                       subset(hcc.ref,hcc.ref$gene %in% common_gene.maqc.vs.hcc)[,c('gene','final')]) %>% unique()
write.csv(common_gene.maqc.vs.hcc.df,file = "HCC1395.vs.Quartet/common_gene.maqc.vs.hcc.df.csv")
df <- common_gene.maqc.vs.hcc.df
df_filtered <- na.omit(df[,c('DEG','final')])
colnames(df_filtered) <- c('A','B')

# 计算一致性：两列值相同的行数
consistent <- sum(df_filtered$A == df_filtered$B)

# 计算总比较行数
total <- nrow(df_filtered)

# 计算一致性概率
consistency_prob <- consistent / total

# 进行二项检验
# 原假设：一致性概率为0.5（随机猜测水平）
binom_result <- binom.test(consistent, total, p = 0.5, alternative = "two.sided")

compare=paste0('MAQC','.vs.','HCC1395')
binom_p=binom_result$p.value

df= data.frame(compare=compare,
               consistency=consistency_prob,
               binom_p=binom_p)
consistency_df = rbind(consistency_df,df)
consistency_df[-1,]

write.csv(consistency_df,file = "HCC1395.vs.Quartet/consistency_df.csv")
# 加载必要的包
library(ggplot2)
library(dplyr)
library(tibble)
library(scales)
library(ggtext)

# 设置Nature风格的颜色和主题
nature_blue <- "#1f77b4"
nature_red <- "#d62728"
nature_colors <- c(nature_blue, nature_red)

# 准备数据 - 添加显著性星号标记
consistency_df <- consistency_df[-1,] %>%
  mutate(
    significance = case_when(
      binom_p < 0.001 ~ "***",
      binom_p < 0.01 ~ "**",
      binom_p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    # 格式化比较名称
    compare = factor(compare, levels = c( "D5/D6.vs.HCC1395", "F7/D6.vs.HCC1395", 
                                         "M8/D6.vs.HCC1395", "MAQC.vs.HCC1395"))
  )

# 创建主图 - 一致性比例条形图
p <- ggplot(consistency_df, aes(x = compare, y = consistency, fill = consistency)) +
  geom_col(width = 0.7, alpha = 0.8) +
  geom_text(aes(label = paste0(round(consistency * 100, 1), "%", significance)), 
            vjust = -0.5, size = 4, fontface = "bold") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_fill_gradient(low = nature_red, high = nature_blue, limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.05), expand = c(0, 0)) +
  labs(
    title = "Comparison Consistency with HCC1395 Reference",
    subtitle = "Proportion of consistent gene expression classifications",
    x = "Comparison Group",
    y = "Consistency Proportion",
    caption = "*** p < 0.001; ** p < 0.01; * p < 0.05\nReference represents perfect consistency (100%)"
  ) +
  my_theme+  theme(
    axis.text.x = element_text(angle = 15, hjust = 1))

# 添加p值注释
p <- p + 
  annotate("text", x = 1.5, y = 0.95, 
           label = paste0("Lowest p-value: ", formatC(min(consistency_df$binom_p[-1]), format = "e", digits = 2)),
           size = 3.5, hjust = 0, color = "gray30")

# 打印图形
print(p)

# 保存图形
ggsave(
  filename = "HCC1395.vs.Quartet/consistency_comparison_nature_style.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

#差异基因柱状图
deg.df <- data.frame(Var1=1,
                     Freq=1,
                     Sample=1)
hcc.deg <- table(hcc.ref$final) %>% as.data.frame()
hcc.deg$Sample <- 'HCC1395'

deg.df <- rbind(deg.df,hcc.deg)

for (i in unique(quartet.ref$compare)){
  quartet.deg <- table(subset(quartet.ref,compare==i)$DEGtype)  %>% as.data.frame()
  quartet.deg$Sample <- paste0('Quartet ',i)
  deg.df <- rbind(deg.df,quartet.deg)
}

deg.df <- deg.df[-1,]
colnames(deg.df) <- c('DEG_type','Counts','Sample')

col_n <- c('#9370DB','#FFC1C1','#EEB4B4','#CD9B9B')
deg.p <- ggplot(deg.df,aes(x=DEG_type,y=Counts))+
  geom_bar(stat = "identity", position = position_dodge(),aes(fill=Sample),alpha = 0.7,color = "black")+
  scale_fill_manual(values = col_n)+
  my_theme;deg.p



#差异表达值分布图
fc.df <- data.frame(FC=c(hcc.ref$meanlogFC,quartet.ref$meanlogFC),
                    Sample=c(rep('HCC1395',length(hcc.ref$meanlogFC)),
                             rep('Quartet',length(quartet.ref$meanlogFC))))

fc.df <- na.omit(fc.df)

col_1 <- c('#9370DB','#FFC1C1')
fc.p <- ggplot(fc.df,aes(x=FC,fill=Sample))+
  geom_histogram(aes(y = after_stat(density)),position = "identity", bins = 100, alpha = 0.7, color = "black")+
  xlim(c(-10,10))+
  scale_fill_manual(values = col_1)+
  my_theme+
  labs( x = "log2(Fold Change)", y = "Frequency");fc.p

library(patchwork)
ref.compare <- deg.p+fc.p+plot_layout(guides = "collect") & theme(legend.position = "bottom",legend.background = element_blank())
ggsave('HCC1395.vs.Quartet/HCC1395.vs.Quartet.png',ref.compare,width=8.15*1.5, height=5.7*1.2,dpi = 300)


