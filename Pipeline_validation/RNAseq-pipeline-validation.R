#RNAseq 分析流程验证------
library(ggplot2)
library(dplyr)
library(patchwork)
#配色网站：https://liuzhaoze.github.io/posts/论文绘图配色/
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Pipeline_validation')
source('/Users/duanshumeng/生物信息/PGx_lab/毕业论文/Scripts/my_theme.R')
colors_1 <- c('#9BC985','#F7D58B','#B595BF','#797BB7')
colors_2 <- c('#62B197','#E18E6D')
dir()
#SNR------
SNR_df <- read.csv("Quartet_validation_SNR.csv")
head(SNR_df)
SNR_df$Quant <-  gsub('Salmon','RefSeq_Salmon',SNR_df$Quant ) %>% 
  gsub('FeatureCounts','Ensembl_STAR_FeatureCounts',.) %>% 
  gsub('HTseq','Ensembl_STAR_HTseq',.) %>% 
  gsub('Stringtie','Ensembl_STAR_Stringtie',.) 
snr.p <- ggplot(SNR_df,aes(Quant,SNR,color=Quant,fill=Quant))+
  geom_bar(stat="summary",fun=mean,position=position_dodge(width = 1))+ 
  geom_hline(yintercept = 12, linetype = "dashed", color = "darkred")+
  theme_bw()+
  scale_fill_manual(values = colors_1)+scale_color_manual(values = colors_1)+
  labs(x = NULL, y = 'SNR',size=12,face="bold") + coord_flip()+ 
  my_theme
#MCC & F1-----
mcc_df <- read.csv("Quartet_validation_MCC_F1.csv")
mcc_df$Quant <- gsub('Salmon\\|edgeR','RefSeq_Salmon_edgeR',mcc_df$quant) %>% gsub('Salmon\\|DEGseq','RefSeq_Salmon_DEGseq',.) %>%
  gsub('FeatureCounts\\|edgeR','Ensembl_STAR_FeatureCounts_edgeR',.) %>% gsub('FeatureCounts\\|DEGseq','Ensembl_STAR_FeatureCounts_DEGseq',.) %>%
  gsub('HTseq\\|edgeR','Ensembl_STAR_HTseq_edgeR',.) %>% gsub('HTseq\\|DEGseq','Ensembl_STAR_HTseq_DEGseq',.) %>% 
  gsub('Stringtie\\|edgeR','Ensembl_STAR_Stringtie_edgeR',.) %>% gsub('Stringtie\\|DEGseq','Ensembl_STAR_Stringtie_DEGseq',.)
head(mcc_df)
f1.p <- ggplot(mcc_df,aes(Quant,F1,color=ref,fill=ref))+
  geom_bar(stat="summary",fun=mean,position=position_dodge(width = 1))+ 
  theme_bw()+
  scale_fill_manual(values = colors_2)+scale_color_manual(values = colors_2)+
  labs(x = NULL, y = 'F1',size=12,face="bold")+ coord_flip() + 
  my_theme;f1.p

mcc.p <- ggplot(mcc_df,aes(Quant,MCC,color=ref,fill=ref))+
  geom_bar(stat="summary",fun=mean,position=position_dodge(width = 1))+ 
  theme_bw()+
  scale_fill_manual(values = colors_2)+scale_color_manual(values = colors_2)+
  geom_hline(yintercept = 0.54, linetype = "dashed", color = "darkred")+
  labs(x = NULL, y = 'MCC',size=12,face="bold") + coord_flip()+ 
  my_theme;mcc.p

val_p <- (snr.p/mcc.p)+plot_annotation(tag_levels = 'A')

ggsave('Validated.byQuartet.pdf',val_p,width=8.15, height=5.7,dpi = 300)
#ggsave('MCC.Quartet.png',mcc.p,width=8.15*1.2, height=5.7*0.8,dpi = 300)

#ggsave('SNR.quartet.png',snr.p,width=8.15*1.2, height=5.7*0.5,dpi = 300)
