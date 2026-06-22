#HCC1395-Ratio-based qualification evoluated by reference datasets
library(stringr)
library(dplyr)
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(patchwork)
my_theme <- theme_bw()+
  theme(plot.title = element_text(hjust=0.5,color = "black",face="plain"),
        axis.text.x =element_text(color = "black",size = 14,face = "bold",angle = 0),
        axis.text.y = element_text(color = "black",size = 14,face = "bold"),
        axis.title.x = element_text(color = "black",face = "bold",size = 16),
        axis.title.y = element_text(color = "black",face = "bold",size = 16),
        #axis.ticks.x = element_text(color = "black",face = "bold",size = 14),
        legend.title = element_text(face = "bold",size = 16),
        legend.text = element_text(face = "bold",size = 14,colour = "black"),
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

setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/High-Confidence-Region')

hcc_ref <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/High-Confidence-Region/All.RefData_DEGs.csv')

seq_lst = list(ARG_PolyA='DNBSEQ-T7',
               BMK_RiboZero='NovaSeq 6000',
               SEQ_PolyA='AVITI',
               SEQ_RiboZero='AVITI',
               NVG_PolyA='NovaSeq 6000',
               MGI_PolyA.M='MGISeq2000',
               MGI_PolyA.T='DNBSEQ-T7',
               MGI_RiboZero='DNBSEQ-T7')

#SNR: SNR_calculate.R
meta_df <- read.csv('hcc1395_metadata.csv')
meta_df <- meta_df[,-1] %>% unique()
meta_df$Batch <- meta_df$batch
#meta_df$batch <- paste0(meta_df$platform,'_',meta_df$protocol)
meta_df$batch <- str_split(meta_df$library,'\\_H',simplify = T)[,1]
snr_df <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/DEG_reference/SNR_FPKM.csv')
snr_df <- merge(snr_df,meta_df[,c('batch','Batch')],all.x = T) %>% unique()
snr_df$library <- str_split(snr_df$batch,'\\_',simplify = T)[,2]
snr_df$library <- gsub('.T','',snr_df$library) %>% gsub('.M','',.)
snr_df$sequencer <- lapply(snr_df$batch,function(x){
  print(x)
  seq_lst[[x]]
}) %>% unlist()

snr_df <- snr_df[order(snr_df$snr),]
snr_df$batch <- factor(snr_df$batch,levels = snr_df$batch)

snr.p = ggplot(snr_df,aes(x=batch,y=snr))+
  geom_bar(stat = "identity", position = position_dodge(),fill = "blue",color='black',,alpha = 0.6)+
  #geom_histogram(stat="count",position = "dodge",alpha = 0.7, color = "black")+
  #scale_fill_manual(values =col_1395)+
  my_theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs( x = "", y = "SNR");snr.p

# 按library分组计算中位数SNR用于排序
library_order <- snr_df %>%
  group_by(library) %>%
  summarise(median_snr = median(snr)) %>%
  arrange(median_snr) %>%
  pull(library)

# 绘制点图+文本标签
snr.p.lib <- ggplot(snr_df, aes(x = factor(library, levels = library_order), y = snr)) +
  geom_point(aes(color = sequencer, shape = sequencer), size = 6) +
  geom_text(aes(label = batch), vjust = -1.5, size = 4, color = "black") +
  labs(
       x = "Library Type",
       y = "SNR") +
  my_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
        legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(15, 16, 17, 18));snr.p.lib

sequencer_order <- snr_df %>%
  group_by(sequencer) %>%
  summarise(median_snr = median(snr)) %>%
  arrange(median_snr) %>%
  pull(sequencer)

# 绘制点图+文本标签
# : 库类型和测序仪的交互效应
snr.p.libseq <- subset(snr_df,sequencer!='MGISeq2000') %>%
  group_by(library, sequencer) %>%
  summarise(mean_SNR = mean(snr), .groups = "drop") %>%
  ggplot(aes(x = library, y = mean_SNR, color = sequencer, group = sequencer)) +
  geom_line(size = 1.2, alpha = 0.7) +
  geom_point(size = 5, aes(shape = sequencer)) +
  geom_text(aes(label = sprintf("%.2f", mean_SNR)), 
            vjust = -1.5, size = 4) +
  labs(x = "Library Type",
       y = "SNR",
       color = "Sequencer",
       shape = "Sequencer") +
  my_theme+
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 16, 17, 18));snr.p.libseq

#MCC
#MCC ----
get_DEGtype <- function(df){
  DEGtype <- apply(df,1,function(x){
    #print(x['adj.P.Val']%>% as.numeric())
    #print(x['logFC'] %>% as.numeric())
    if (x['P.Value'] %>% as.numeric() >= 0.05){
      return('non-DEG')
    }else{
      if (x['logFC'] %>% as.numeric() > 1){
        return('up-regulate')
      } else if (x['logFC'] %>% as.numeric()  < -1){
        return('down-regulate')
      } else {
        return('non-DEG')
      }
    }
  }) %>% unlist()
}

get_mcc <- function(tn,tp,fn,fp){
  print(paste(tn,tp,fn,fp,sep = ' '))
  
  sqrt_term = sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc = (tp*tn - fp*fn)/sqrt_term
  return(mcc)
}

if (FALSE) {get_f1 <- function(tn,tp,fn,fp){
  print(paste(tn,tp,fn,fp,sep = ' '))
  
  precision = tp/(tp+fp)
  recall = tp/(tp+fn)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  return(f1_score)
}
}
get_mcc_all <- function(test_expr,ref_expr,smp){
  TN_non_non = intersect(intersect(subset(test_expr,DEGtype=='non-DEG')$gene_id,ref_expr$gene),
                         subset(ref_expr,DEGtype=='non-DEG')$gene) %>% length() 
  FN_non_up = intersect(intersect(subset(test_expr,DEGtype=='non-DEG')$gene_id,ref_expr$gene),
                        subset(ref_expr,DEGtype=='up-regulate')$gene) %>% length() 
  
  FN_non_down = intersect(intersect(subset(test_expr,DEGtype=='non-DEG')$gene_id,ref_expr$gene),
                          subset(ref_expr,DEGtype=='down-regulate')$gene) %>% length()
  
  FP_up_non = intersect(intersect(subset(test_expr,DEGtype=='up-regulate')$gene_id,ref_expr$gene),
                        subset(ref_expr,DEGtype=='non-DEG')$gene) %>% length() 
  TP_up_up = intersect(intersect(subset(test_expr,DEGtype=='up-regulate')$gene_id,ref_expr$gene),
                       subset(ref_expr,DEGtype=='up-regulate')$gene) %>% length()
  FP_up_down = intersect(intersect(subset(test_expr,DEGtype=='up-regulate')$gene_id,ref_expr$gene),
                         subset(ref_expr,DEGtype=='down-regulate')$gene) %>% length()
  
  FP_down_non = intersect(intersect(subset(test_expr,DEGtype=='down-regulate')$gene_id,ref_expr$gene),
                          subset(ref_expr,DEGtype=='non-DEG')$gene) %>% length() 
  
  FP_down_up = intersect(intersect(subset(test_expr,DEGtype=='down-regulate')$gene_id,ref_expr$gene),
                         subset(ref_expr,DEGtype=='up-regulate')$gene) %>% length() 
  
  TP_down_down = intersect(intersect(subset(test_expr,DEGtype=='down-regulate')$gene_id,ref_expr$gene),
                           subset(ref_expr,DEGtype=='down-regulate')$gene) %>% length() 
  
  #FN_not_up = intersect(setdiff(test_expr$gene_id,ref_expr$gene),
  #                      subset(ref_expr,DEGtype=='up-regulate')$gene) %>% length() 
  
  #FN_not_down = intersect(setdiff(test_expr$gene_id,ref_expr$gene),
  #                       subset(ref_expr,DEGtype=='down-regulate')$gene) %>% length() 
  #grep("FN", ls(), value = TRUE)
  FN = FN_non_down + FN_non_up 
  #+ FN_not_down+FN_not_up
  #grep("TN", ls(), value = TRUE)
  TN = TN_non_non
  #grep('FP',ls(),value = TRUE)
  FP = FP_down_non + FP_down_up + FP_up_down+FP_up_non
  #grep('TP',ls(),value = TRUE)
  TP = TP_down_down+TP_up_up
  df = data.frame(TN=TN,
                  TP=TP,
                  FN=FN,
                  FP=FP,
                  Sample=smp)
  #mcc <- get_mcc(TN*0.01,TP*0.01,FN*0.01,FP*0.01)
  return(df)
  
}

MCC_df.all = data.frame(TN=1,
                    TP=1,
                    FN=1,
                    FP=1,
                    Sample=1,
                    MCC=1,
                    Tool=1,
                    PCC=1)

for (c in dir(pattern = "*.multibatch_count.detectable.DEG_limma.csv")){
  tool = str_split(c,'\\.',simplify = T)[,1]
  #print(tool)
  fc_expr <- read.csv(c)
  batchs <- fc_expr$batch %>% unique()
  mcc_df = data.frame(TN=1,
                      TP=1,
                      FN=1,
                      FP=1,
                      Sample=1,
                      PCC=1)
  for (i in batchs){
    print(i)
    test_expr <- subset(fc_expr,batch==i)
    test_expr$gene_id <- test_expr$gene
    test_expr$DEGtype <- test_expr$type
    hcc_ref$gene <- str_split(hcc_ref$gene_compare,'\\ ',simplify = T)[,1]
    hcc_ref$DEGtype <- hcc_ref$final
    mcc = get_mcc_all(test_expr,hcc_ref,i)
    hcc_ref_pcc <- merge(test_expr,hcc_ref,by='gene')
    hcc_ref_pcc$ref_logfc <- log2(hcc_ref_pcc$fc_mean+0.001)
    pcc = cor.test(hcc_ref_pcc$logfc, hcc_ref_pcc$ref_logfc,method = "pearson")$estimate
    mcc$PCC=pcc
    mcc_df <- rbind(mcc_df,mcc)

    
  }
  MCC.df <- mcc_df[-1,]
  MCC.df$MCC <- apply(MCC.df,1,function(x){
    #print(x)
    get_mcc(x['TN'] %>% as.numeric(),
            x['TP']%>% as.numeric(),
            x['FN']%>% as.numeric(),
            x['FP']%>% as.numeric())
  }) %>% unlist()
  MCC.df$Tool <- tool
  MCC_df.all <- rbind(MCC_df.all,MCC.df)
}

MCC_df.all <- MCC_df.all[-1,]



MCC_df.all <- merge(MCC_df.all,snr_df[,c('batch','Batch')],by.x = 'Sample',by.y = 'Batch')
MCC_df.all$batch <- as.character(MCC_df.all$batch)

MCC_df.all$library <- str_split(MCC_df.all$batch,'\\_',simplify = T)[,2]
MCC_df.all$library <- gsub('.T','',MCC_df.all$library) %>% gsub('.M','',.)
MCC_df.all$sequencer <- lapply(MCC_df.all$batch,function(x){
  print(x)
  seq_lst[[x]]
}) %>% unlist()

library(ggplot2)
library(dplyr)
library(patchwork)  # 用于组合多个图表
library(RColorBrewer)  # 用于颜色方案

# 创建技术组合标签
MCC_df.all$combo <- paste(MCC_df.all$Tool, MCC_df.all$library, MCC_df.all$sequencer, sep = " + ")

#MCC------
tool_order <- MCC_df.all %>%
  group_by(Tool) %>%
  summarise(median_MCC = median(MCC)) %>%
  arrange(median_MCC) %>%
  pull(Tool)

# 绘制箱线图+散点图
mcc.p.tool <- ggplot(MCC_df.all, aes(x = factor(Tool, levels = tool_order), y = MCC)) +
  geom_boxplot(aes(fill = Tool), alpha = 0.7, outlier.shape = NA,show.legend = FALSE) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.5) +
  labs( x = "Tool",
       y = "MCC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_brewer(palette = "Set1", name = "Batch")+
  stat_compare_means(method = "t.test",
                     #aes(label = "p.format"), #显示方式
                     label = "p.signif",
                     size = 5,label.y=c(0.93,0.96,0.97),
                     comparisons = list(c('Salmon','Stringtie'),c('Salmon','HTseq'),c('Salmon','FeatureCounts')))+my_theme;mcc.p.tool

# 计算每个批次的中位数MCC用于排序
batch_order <- MCC_df.all %>%
  group_by(batch) %>%
  summarise(median_MCC = median(MCC)) %>%
  arrange(median_MCC) %>%
  pull(batch)

# 绘制箱线图+散点图
mcc.p.batch <- ggplot(MCC_df.all, aes(x = factor(batch, levels = batch_order), y = MCC)) +
  geom_boxplot(aes(fill = batch), alpha = 0.7, outlier.shape = NA,show.legend = FALSE) +
  geom_jitter(width = 0.2, size = 3,alpha=0.5) +
  labs(x = "Batch",
       y = "MCC") +
  my_theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = 1:6);mcc.p.batch


lib_order <- MCC_df.all %>%
  group_by(library) %>%
  summarise(median_MCC = median(MCC)) %>%
  arrange(median_MCC) %>%
  pull(library)

# 创建绘图
mcc.p.lib <- ggplot(subset(MCC_df.all,sequencer !='MGISeq2000'), aes(x = factor(library, levels = lib_order), y = MCC)) +
  geom_boxplot(aes(fill = library), alpha = 0.6, width = 0.6, show.legend = FALSE) +
  geom_jitter(size = 3, width = 0.15, height = 0, alpha = 0.5) +
  labs(x = "Library Type",
       y = "MCC") +
  my_theme+
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # 不同形状表示测序仪
  ylim(min(MCC_df.all$MCC)*0.998, max(MCC_df.all$MCC)*1.001)+
  stat_compare_means(method = "t.test",
                     #aes(label = "p.format"), #显示方式
                     label = "p.signif",
                     size = 5,label.y=c(0.95),
                     comparisons = list(c('PolyA','RiboZero')));mcc.p.lib


# 计算每个测序仪的中位数PCC用于排序
seq_order <- MCC_df.all %>%
  group_by(sequencer) %>%
  summarise(median_MCC = median(MCC)) %>%
  arrange(median_MCC) %>%
  pull(sequencer)

seq_order <- c('NovaSeq 6000','DNBSEQ-T7','AVITI','MGISeq2000')
# 创建绘图
mcc.p.seq <- ggplot(subset(MCC_df.all,sequencer !='MGISeq2000'), aes(x = factor(sequencer, levels = seq_order), y = MCC)) +
  geom_boxplot(aes(fill = sequencer), alpha = 0.6, width = 0.6, show.legend = FALSE) +
  geom_jitter(size = 3, width = 0.15, height = 0, alpha = 0.5) +
  labs(x = "Sequencer",
       y = "MCC") +
  my_theme+
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(16, 17))+
  stat_compare_means(method = "t.test",
                     #aes(label = "p.format"), #显示方式
                     label = "p.signif",
                     size = 5,label.y=c(0.94,0.96,0.95),
                     comparisons = list(c('NovaSeq 6000','DNBSEQ-T7'),c('NovaSeq 6000','AVITI'),c('DNBSEQ-T7','AVITI')))
mcc.p.seq

# : 库类型和测序仪的交互效应
mcc.p.libseq <- subset(MCC_df.all,sequencer !='MGISeq2000') %>%
  group_by(library, sequencer) %>%
  summarise(mean_MCC = mean(MCC), .groups = "drop") %>%
  ggplot(aes(x = library, y = mean_MCC, color = sequencer, group = sequencer)) +
  geom_line(size = 1.2, alpha = 0.7) +
  geom_point(size = 5, aes(shape = sequencer)) +
  geom_text(aes(label = sprintf("%.4f", mean_MCC)), 
            vjust = -1.5, size = 4) +
  labs(x = "Library Type",
       y = "MCC",
       color = "Sequencer",
       shape = "Sequencer") +
  my_theme+
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 16, 17, 18))


# 计算每个组合的中位数PCC用于排序
combo_mcc <- MCC_df.all %>%
  group_by(combo, Tool) %>%  # 按组合和工具分组
  summarise(mean_MCC = mean(MCC),  # 计算平均MCC
            .groups = "drop") %>%  # 取消分组
  arrange(mean_MCC) %>%  # 按平均MCC排序
  mutate(combo = factor(combo, levels = unique(combo)))  # 将组合转换为因子以保持排序

mcc.p.order <- ggplot(combo_mcc, aes(x = combo, y = mean_MCC)) +
  geom_bar(aes(fill = Tool), stat = "identity", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", mean_MCC)), 
            vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Combination (Tool + Library + Sequencer)",
       y = "MCC",
       fill = "Tool") +  # 添加图例标题
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",  # 将图例放在底部
        panel.grid.major.x = element_blank()) +
  scale_fill_brewer(palette = "Set1") +
  ylim(0, max(combo_mcc$mean_MCC)*1.1)


#PCC------------
tool_order <- MCC_df.all %>%
  group_by(Tool) %>%
  summarise(median_PCC = median(PCC)) %>%
  arrange(median_PCC) %>%
  pull(Tool)

# 绘制箱线图+散点图
pcc.p.tool <- ggplot(MCC_df.all, aes(x = factor(Tool, levels = tool_order), y = PCC)) +
  geom_boxplot(aes(fill = Tool), alpha = 0.7, outlier.shape = NA,show.legend = FALSE) +
  geom_jitter( width = 0.2, size = 3, alpha = 0.5) +
  labs(x = "Tool",
       y = "PCC") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_brewer(palette = "Set1", name = "Batch")+my_theme+
  stat_compare_means(method = "t.test",
                     #aes(label = "p.format"), #显示方式
                     label = "p.signif",
                     size = 5,label.y=c(0.997,0.998,0.999),
                     comparisons = list(c('Salmon','Stringtie'),c('Salmon','HTseq'),c('Salmon','FeatureCounts')))+my_theme;mcc.p.tool
pcc.p.tool

# 计算每个批次的中位数MCC用于排序
batch_order <- MCC_df.all %>%
  group_by(batch) %>%
  summarise(median_PCC = median(PCC)) %>%
  arrange(median_PCC) %>%
  pull(batch)

# 绘制箱线图+散点图
pcc.p.batch <- ggplot(MCC_df.all, aes(x = factor(batch, levels = batch_order), y = PCC)) +
  geom_boxplot(aes(fill = batch), alpha = 0.7, outlier.shape = NA,show.legend = FALSE) +
  geom_jitter(width = 0.2, size = 3,alpha=0.5) +
  labs(x = "Batch",
       y = "PCC") +
  my_theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_shape_manual(values = 1:6);pcc.p.batch

lib_order <- MCC_df.all %>%
  group_by(library) %>%
  summarise(median_PCC = median(PCC)) %>%
  arrange(median_PCC) %>%
  pull(library)

# 创建绘图
pcc.p.lib <- ggplot(subset(MCC_df.all,sequencer != 'MGISeq2000'), aes(x = factor(library, levels = lib_order), y = PCC)) +
  geom_boxplot(aes(fill = library), alpha = 0.6, width = 0.6, show.legend = FALSE) +
  geom_jitter(
              size = 3, width = 0.15, height = 0, alpha = 0.5) +
  labs(x = "Library Type",
       y = "PCC") +
  my_theme+
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(15, 16, 17, 18)) +  # 不同形状表示测序仪
  ylim(min(MCC_df.all$PCC)*0.998, max(MCC_df.all$PCC)*1.001)+ stat_compare_means(method = "t.test",
                                                                                  #aes(label = "p.format"), #显示方式
                                                                                  label = "p.signif",
                                                                                  size = 5,label.y=c(0.997),
                                                                                  comparisons = list(c('PolyA','RiboZero')));pcc.p.lib

# 计算每个测序仪的中位数PCC用于排序
seq_order <- MCC_df.all %>%
  group_by(sequencer) %>%
  summarise(median_PCC = median(PCC)) %>%
  arrange(median_PCC) %>%
  pull(sequencer)

# 创建绘图
pcc.p.seq <- ggplot(subset(MCC_df.all,sequencer != 'MGISeq2000'), aes(x = factor(sequencer, levels = seq_order), y = PCC)) +
  geom_boxplot(aes(fill = sequencer), alpha = 0.6, width = 0.6, show.legend = FALSE) +
  geom_jitter(size = 3, width = 0.15, height = 0, alpha = 0.5) +
  labs(x = "Sequencer",
       y = "PCC") +
  my_theme+
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(16, 17))+  stat_compare_means(method = "t.test",
                                                              #aes(label = "p.format"), #显示方式
                                                              label = "p.signif",
                                                              size = 5,label.y=c(0.997,0.998,0.999),
                                                              comparisons = list(c('NovaSeq 6000','DNBSEQ-T7'),c('NovaSeq 6000','AVITI'),c('DNBSEQ-T7','AVITI')))
pcc.p.seq

# : 库类型和测序仪的交互效应
pcc.p.libseq <- subset(MCC_df.all,sequencer != 'MGISeq2000') %>%
  group_by(library, sequencer) %>%
  summarise(mean_PCC = mean(PCC), .groups = "drop") %>%
  ggplot(aes(x = library, y = mean_PCC, color = sequencer, group = sequencer)) +
  geom_line(size = 1.2, alpha = 0.7) +
  geom_point(size = 5, aes(shape = sequencer)) +
  geom_text(aes(label = sprintf("%.4f", mean_PCC)), 
            vjust = -1.5, size = 4) +
  labs(x = "Library Type",
       y = "PCC",
       color = "Sequencer",
       shape = "Sequencer") +
  my_theme+
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 16, 17, 18))


# 计算每个组合的中位数PCC用于排序
combo_pcc <- MCC_df.all %>%
  group_by(combo, Tool) %>%  # 按组合和工具分组
  summarise(mean_PCC = mean(PCC),  # 计算平均MCC
            .groups = "drop") %>%  # 取消分组
  arrange(mean_PCC) %>%  # 按平均MCC排序
  mutate(combo = factor(combo, levels = unique(combo)))  # 将组合转换为因子以保持排序

pcc.p.order <- ggplot(combo_pcc, aes(x = combo, y = mean_PCC)) +
  geom_bar(aes(fill = Tool), stat = "identity", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", mean_PCC)), 
            vjust = -0.5, size = 3.5, color = "black") +
  labs(x = "Combination (Tool + Library + Sequencer)",
       y = "PCC",
       fill = "Tool") +  # 添加图例标题
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "bottom",  # 将图例放在底部
        panel.grid.major.x = element_blank()) +
  scale_fill_brewer(palette = "Set1") +
  ylim(0, max(combo_pcc$mean_PCC)*1.1)



#计算MGI 5G数据量下的MCC----------------
system(paste0('Rscript /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Reference_standard_protocol-main/scripts/deg_limma-hcc1395.R',
              ' -i /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/DEG_reference/Salmon.gene_quant-MGI-5G.csv',
              ' -m ','/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/DEG_reference/MGI-5G-metadata.csv',' -o ','/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/DEG_reference/'))

mcc_df = data.frame(TN=1,
                    TP=1,
                    FN=1,
                    FP=1,
                    Sample=1,
                    PCC=1)
test_expr <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/DEG_reference/Salmon.gene_quant-MGI-5G.DEG_limma.csv')
test_expr$gene_id <- test_expr$gene
test_expr$DEGtype <- test_expr$type
hcc_ref$gene <- str_split(hcc_ref$gene_compare,'\\ ',simplify = T)[,1]
hcc_ref$DEGtype <- hcc_ref$final
mcc = get_mcc_all(test_expr,hcc_ref,'MGI-5G')
hcc_ref_pcc <- merge(test_expr,hcc_ref,by='gene')
hcc_ref_pcc$ref_logfc <- log2(hcc_ref_pcc$fc_mean+0.001)
pcc = cor.test(hcc_ref_pcc$logfc, hcc_ref_pcc$ref_logfc,method = "pearson")$estimate
mcc$PCC=pcc
mcc_df <- rbind(mcc_df,mcc)
MCC.df <- mcc_df[-1,]
MCC.df$MCC <- apply(MCC.df,1,function(x){
  #print(x)
  get_mcc(x['TN'] %>% as.numeric(),
          x['TP']%>% as.numeric(),
          x['FN']%>% as.numeric(),
          x['FP']%>% as.numeric())
}) %>% unlist()

#差异表达基因占比:

intersect(test_expr$gene,hcc_ref$gene) %>% length()
(intersect(test_expr$gene,hcc_ref$gene) %>% length())/(length(hcc_ref$gene))*100

for (c in dir(pattern = "*.multibatch_count.detectable.DEG_limma.csv")){
  print(c)
  test_expr <- read.csv(c)
  intersect(test_expr$gene,hcc_ref$gene) %>% length() %>% print()
  ((intersect(test_expr$gene,hcc_ref$gene) %>% length())/(length(hcc_ref$gene))*100) %>% print()
}

ggplot(MCC.df,aes(x=Sample,y=MCC,fill=Sample))+
  geom_bar(stat = "identity", position = position_dodge(),color = "black",,alpha = 0.6)+
  #geom_histogram(stat="count",position = "dodge",alpha = 0.7, color = "black")+
  #scale_fill_manual(values =col_1395)+
  my_theme+
  #theme(axis.text.x =element_text(color = "black",size = 12,face = "bold",angle = 45,hjust = 1))+
  labs( x = "", y = "PCC")


all.p.batch <- mcc.p.batch/pcc.p.batch+plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='bottom',
                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))
all.p.batch
ggsave('HCC1395.RNAseq_performance-batch.png',all.p.batch,width =8.15*1.5 ,height =5.7*2.5,dpi = 300)
ggsave('HCC1395.RNAseq_performance-batch.pdf',all.p.batch,width =8.15*1.5 ,height =5.7*2.5,dpi = 300)

all.p.all <- (snr.p+snr.p.libseq)/(mcc.p.batch|mcc.p.tool|mcc.p.libseq)/
  (pcc.p.batch|pcc.p.tool|pcc.p.libseq)/(mcc.p.order+pcc.p.order)+plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='bottom',
                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))

ggsave('HCC1395.RNAseq_performance-all.png',all.p.all,width =8.15*2.5 ,height =5.7*3.5,dpi = 300)

#仅比较测序仪和文库构建的影响----
all.p.libseq <- (mcc.p.lib+mcc.p.seq+plot_layout(widths = c(1,2)))/(pcc.p.lib+pcc.p.seq+plot_layout(widths = c(1,2)))/(mcc.p.libseq+pcc.p.libseq)+
  plot_annotation(tag_levels = 'a')+
  plot_layout(guides='collect') & theme(legend.position='bottom',
                                        plot.tag = element_text(color = "black",size = 18,face = "bold"))


ggsave('HCC1395.RNAseq_performance-libseq.png',all.p.libseq,width =8.15*1.5 ,height =5.7*2.5,dpi = 300)
ggsave('HCC1395.RNAseq_performance-libseq.pdf',all.p.libseq,width =8.15*1.5 ,height =5.7*2.5,dpi = 300)

#比较raw count和fc之间的差异
#df_count.dete <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/DEG_reference/Stringtie.gene_quant.FPKM.csv')
df_count.dete <- read.csv('FeatureCounts.gene_quant.csv')
mt_df <- read.csv('hcc1395_metadata.csv')
mt_df <- mt_df[,-1] %>% unique()
df_count.dete$gene <- df_count.dete$Ensembl
df_count.dete[df_count.dete == 0] <- NA
df_count.dete.melt <- df_count.dete %>% melt()
df_count.dete.melt <- merge(df_count.dete.melt[,c('gene','variable','value')],mt_df[,c('library','sample','sample_rep','batch')],by.x = 'variable',by.y = 'library') %>% unique()

hc_ref <- read.csv('All.RefData_DEGs.csv')
hc_ref$gene <- str_split(hc_ref$gene_compare,'\\ ',simplify = T)[,1]

#3. 计算差异表达基因
result <- data.frame()
df <- subset(df_count.dete.melt,df_count.dete.melt$gene %in% hc_ref$gene)
# 获取所有 unique 的 sample_rep 前缀（例如 HCC1395_1, HCC1395_2, HCC1395_3）
sample_prefixes <- unique(df$sample_rep)
# 筛选出HCC1395和HCC1395BL配对的样本
sample_pairs <- gsub("HCC1395", "HCC1395BL", sample_prefixes)

# 循环遍历每一对样本，计算Fold Change
for (sample_prefix in sample_prefixes) {
  
  # 找出当前 sample_prefix 对应的 HCC1395 和 HCC1395BL 的数据
  df_hcc1395 <- df[df$sample_rep == sample_prefix, ]
  df_hcc1395bl <- df[df$sample_rep == gsub("HCC1395", "HCC1395BL", sample_prefix), ]
  
  # 将 HCC1395 和 HCC1395BL 按照 gene(X) 合并
  df_merged <- merge(df_hcc1395, df_hcc1395bl, by = "gene", suffixes = c("_HCC1395", "_HCC1395BL"))
  
  # 计算Fold Change
  df_merged$log2fc <- df_merged$value_HCC1395 / df_merged$value_HCC1395BL
  
  # 把结果加入到结果数据框中
  result <- rbind(result, df_merged)
}


result.sub <- result[,c('gene','value_HCC1395','value_HCC1395BL','batch_HCC1395','log2fc')] %>% unique()
colnames(result.sub) <- c('gene','HCC1395','HCC1395BL','batch','logfc')
#write.csv(result.sub,file=gsub('.csv','.log2fc.csv',raw_count),row.names = F)

result.sub.melt <- result.sub[,c('gene','HCC1395','HCC1395BL','logfc')] %>% group_by(gene) %>%
  dplyr::summarise(HCC1395_cv=abs(sd(HCC1395)/mean(HCC1395)),
                   HCC1395BL_cv=abs(sd(HCC1395BL)/mean(HCC1395BL)),
                   logfc_cv=abs(sd(logfc)/mean(logfc))) %>% melt()


p1 <- ggplot(result.sub.melt,aes(x=variable,y=value))+
    #geom_violin(aes(fill = variable), width = 0.9,alpha=0.5)+
    geom_boxplot(aes(fill = variable), width = 0.2,alpha=0.7) +  # 箱线图展示分布
    #geom_jitter(aes(color = variable), width = 0.1, size = 3) +  # 散点展示个体值
    labs(
      x = 'Gene Quantification',
      y = "Standard Variation"
    ) +
    my_theme+ guides(fill=FALSE,color=FALSE)+
    stat_compare_means(method = "wilcox.test",
                       #aes(label = "p.format"), #显示方式
                       label = "p.signif",
                       size = 5,
                       comparisons = list(c('HCC1395_sd','logfc_sd'),c('HCC1395BL_sd','logfc_sd')));p1
