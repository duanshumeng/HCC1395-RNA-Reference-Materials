setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/High-Confidence-Region')
#构建高置信参考数据集----
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(patchwork)
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

#RNA层面 TNBC标志物差异性分析，确定是三阴性的----
#雌激素（ER）->ESR1，孕激素（ER）->PGR，人表皮生长因子受体（HER-2）->ERBB2
gene_df <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/基因组和转录组关联分析/Ensembl2RefSeq.txt',sep = '\t',header = F)
colnames(gene_df) <- c('gene','gene_name')
tnbc_gene <- subset(gene_df,gene_name %in% c('ESR1','PGR','ERBB2','ACTB'))

hcc_tpm <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/DEG_reference/Salmon.gene_quant.TPM.csv',sep = ',')
hcc_tpm.tnbc <- subset(hcc_tpm,RefSeq %in% c('ESR1','PGR','ERBB2','ACTB'))


hcc_seqc_tmp <- read.csv('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/RNA_editing/Salmon.gene_quant.TPM.SEQC2.csv',sep = ',')
hcc_seqc_tmp.tnbc <- subset(hcc_seqc_tmp,RefSeq %in% c('ESR1','PGR','ERBB2','ACTB'))


# 加载必要的包
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rstatix)


get_tnbc_expression <- function(hcc_tpm.tnbc){
  rownames(hcc_tpm.tnbc) <- hcc_tpm.tnbc$RefSeq
  
  hcc_tpm.tnbc <- hcc_tpm.tnbc[,grepl('HCC1395',colnames(hcc_tpm.tnbc))]
  sample_groups <- ifelse(grepl('HCC1395BL',colnames(hcc_tpm.tnbc)),'HCC1395BL','HCC1395')
  #sample_groups <- sample_groups[!sample_groups %in% c("RefSeq", "Ensembl")] # 移除非样本列
  
  # 为每个样本创建分组信息的数据框
  sample_info <- data.frame(
    sample = colnames(hcc_tpm.tnbc), # 排除第一列(RefSeq)和最后一列(Ensembl)
    group = sample_groups
  )
  
  # 将表达矩阵转换为长格式
  expression_long <- hcc_tpm.tnbc %>%
    rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "sample",
      values_to = "expression"
    ) %>%
    left_join(sample_info, by = "sample")
  
  # 定义要绘制的基因列表
  library(ggplot2)
  library(ggpubr)
  library(tidyverse)
  library(rstatix)
  genes_of_interest <- c("ERBB2", "ESR1", "PGR",'ACTB')
  
  # 创建函数来生成每个基因的图形
  create_gene_plot <- function(gene_name) {
    # 筛选特定基因的数据
    gene_data <- expression_long %>% filter(gene == gene_name)
    
    # 进行统计学检验
    stat_test <- gene_data %>%
      anova_test(expression ~ group) %>%
      adjust_pvalue(method = "bonferroni") %>%
      add_significance()
    
    # 进行两两比较
    pairwise_test <- gene_data %>%
      pairwise_t_test(expression ~ group, p.adjust.method = "bonferroni")
    
    # 创建图形
    p <- ggplot(gene_data, aes(x = group, y = expression, fill = group)) +
      geom_boxplot(
        width = 0.6,
        alpha = 0.8,
        outlier.shape = 16,
        outlier.size = 2,
        outlier.alpha = 0.5
      ) +
      geom_jitter(
        width = 0.15,
        size = 1.5,
        alpha = 0.6
      ) +
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 18,
        size = 3,
        color = "red"
      ) +
      labs(
        title = paste0(gene_name, " Expression"),
        x = "Cell Type",
        y = "TPM",
        subtitle = paste("ANOVA, p =", format(stat_test$p, scientific = TRUE, digits = 3))
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      scale_fill_brewer(palette = "Set2")
    
    # 添加统计学显著性标记
    if (nrow(pairwise_test %>% filter(p.adj.signif != "ns")) > 0) {
      p <- p + 
        stat_pvalue_manual(
          pairwise_test %>% filter(p.adj.signif != "ns"),
          y.position = max(gene_data$expression) * 1.05,
          step.increase = 0.1,
          label = "p.adj.signif"
        )
    }
    
    return(p)
  }
  p_lst = lst()
  # 为每个基因创建图形并保存
  for (gene in genes_of_interest) {
    p_lst[[gene]] = create_gene_plot(gene)
  }
  p_tnbc <- p_lst[[1]]|p_lst[[2]]|p_lst[[3]]|p_lst[[4]]
  return(p_tnbc)
}
p_tnbc_pgx <- get_tnbc_expression(hcc_tpm.tnbc)
p_tnbc_seqc2 <- get_tnbc_expression(hcc_seqc_tmp.tnbc)
ggsave('TNBC_gene_expression.PGx.tiff',p_tnbc_pgx,height =8.15*0.6,width =5.7*1.5,dpi = 300)
ggsave('TNBC_gene_expression.SEQC2.tiff',p_tnbc_seqc2,height =8.15*0.6,width =5.7*1.5,dpi = 300)

#DNA层面----
library(gggenes) # 用于绘制基因结构
library(ggrepel) # 避免标签重叠
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/High_confidence_SNVs&Indel/PGx/')
#参考https://mp.weixin.qq.com/s/BkHQmU4FPUaiouPdXlGb3g
library('maftools')
library('lolliplot')
#BRCA-----
BRCA_maf <- annovarToMaf("/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/High_confidence_SNVs&Indel/PGx/HCC1395BL_BRCA.variation.txt",
             Center = NULL,  refBuild = "hg38", 
             tsbCol = ' Tumor_Sample_Barcode',  
             table = "refGene",  
             ens2hugo = TRUE,  basename = "HCC1395BL.BRCA",  sep = "\t",  MAFobj = FALSE,  sampleAnno = NULL)
BRCA_maf$Tumor_Sample_Barcode <- 'HCC1395BL'
BRCA_laml = read.maf(maf = BRCA_maf)
pdf("BRCA1.pdf",, width=6, height=6)
BRCA1_p <- lollipopPlot(maf = BRCA_laml, gene = 'BRCA1', AACol = 'AAChange.refGene', showMutationRate = FALSE)
dev.off()
pdf("BRCA2.pdf",, width=6, height=6)
BRCA2_p <- lollipopPlot(maf = BRCA_laml, gene = 'BRCA2', AACol = 'AAChange.refGene', showMutationRate = FALSE)
dev.off()
#TP53-----
TP53_maf <- annovarToMaf("/Users/duanshumeng/生物信息/PGx_lab/HCC1395/SEQC2_analysis/hcc1395_analysis_pipeline/High_confidence_SNVs&Indel/PGx/HCC1395_TP53_variants.txt",
                         Center = NULL,  refBuild = "hg38", 
                         tsbCol = ' Tumor_Sample_Barcode',  
                         table = "refGene",  
                         ens2hugo = TRUE,  basename = "HCC1395.TP53",  sep = "\t",  MAFobj = FALSE,  sampleAnno = NULL)
TP53_maf$Tumor_Sample_Barcode <- 'HCC1395'
TP53_laml = read.maf(maf = TP53_maf)
pdf("TP53.pdf",, width=6, height=6)
TP53_p <- lollipopPlot(maf = TP53_laml, gene = 'TP53', AACol = 'AAChange.refGene', pointSize=1.5,
                       showMutationRate = FALSE)
dev.off()

setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/High-Confidence-Region')
#与Quartet的标准参考集进行比较--------
hcc.ref <- read.csv('All.RefData_DEGs.csv',header = TRUE)
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
                                       subset(hcc.ref,hcc.ref$gene %in% common_gene.quartet.vs.hcc)[,c('gene','meanlogFC')])
write.csv(common_gene.quartet.vs.hcc.df,file = "HCC1395.vs.Quartet/ccommon_gene.quartet.vs.hcc.df.csv")
consistency_df <- data.frame(compare=1,
                             consistency=1,
                             binom_p=0)

#logFC 分布--------------
# 2. 数据预处理
# 处理 HCC1395 数据
df_hcc <- data.frame(
  meanlogFC = hcc.ref$meanlogFC,
  Source = "HCC1395",
  CompareGroup = "HCC1395/HCC1395BL" # HCC只有一个组
)

# 处理 Quartet 数据
df_quartet <- data.frame(
  meanlogFC = quartet.ref$meanlogFC,
  Source = "Quartet",
  CompareGroup = quartet.ref$compare # 这里包含三个比较组
)

# 合并数据
plot_data <- rbind(df_hcc, df_quartet)

# 3. 绘图
compare_colors = c('#fb6a4b','#9e9bc7','#015932','#fdbb84','#4292c7')
names(compare_colors)=c("HCC1395/HCC1395BL","D5/D6","F7/D6","M8/D6","A/B")
p_logFC <- ggplot(plot_data, aes(x = meanlogFC, fill = CompareGroup, pattern = Source)) +
  # 使用密度图展现分布，也可以换成 geom_histogram_pattern
  geom_density(alpha = 0.5) +
  # 设置填充方式：HCC1395用无填充(none)，Quartet用斜线(stripe)
  # 设置颜色：你可以根据需要自定义四个组的颜色
  scale_fill_manual(values = compare_colors)+
  labs(
    title = "Distribution of meanlogFC",
    x = "meanlogFC",
    y = "Density",
    fill = "Comparison Group",
    pattern = "Data Source"
  ) +
  my_theme


# 加载必要的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# 数据预处理 - 转换为长格式
clean_data_long <- common_gene.quartet.vs.hcc.df %>%
  drop_na() %>%
  pivot_longer(
    cols = c("D5/D6", "F7/D6", "M8/D6"),
    names_to = "sample_type",
    values_to = "quartet_value"
  )

# 计算相关系数
correlations <- clean_data_long %>%
  group_by(sample_type) %>%
  summarise(
    # 使用broom包获取相关性检验的完整结果
    tidy_cor = list(broom::tidy(cor.test(
      quartet_value, 
      meanlogFC, 
      method = "pearson",
      use = "complete.obs"
    ))),
    .groups = "drop"
  ) %>%
  # 展开相关系数和p值
  unnest_wider(tidy_cor) %>%
  # 重命名列名
  rename(
    correlation = estimate,
    p_value = p.value,
    method = method
  ) %>%
  # 添加置信区间
  mutate(
    conf.low = conf.low,
    conf.high = conf.high
  ) %>%
  # 选择需要的列
  select(sample_type, correlation, p_value, method, conf.low, conf.high)



# 首先在correlations数据框中创建一个合并的标签
correlations <- correlations %>%
  mutate(
    # 合并r和p值到一个标签
    label = paste0(sample_type, ": r = ", 
                   round(correlation, 3), 
                   ", p = ", 
                   format.pval(p_value, digits = 2, eps = 0.001)),
    # 如果需要科学计数法
    label_sci = ifelse(p_value < 0.001, 
                       paste0(sample_type, ": r = ", round(correlation, 3), ", p < 0.001"),
                       paste0(sample_type, ": r = ", round(correlation, 3), ", p = ", round(p_value, 3)))
  )

# 创建Nature风格的单一散点图
p_hcc_quartet <- ggplot(clean_data_long, aes(x = meanlogFC, y = quartet_value, color = sample_type)) +
  # 添加散点
  geom_point(
    alpha = 0.8, 
    size = 1.5,
    shape = 21,
    aes(fill = sample_type),
    color = "white",
    stroke = 0.5
  ) +
  
  # 为每个样本类型添加回归线
  geom_smooth(
    method = "lm", 
    se = TRUE, 
    aes(fill = sample_type), 
    alpha = 0.1,
    size = 0.8
  ) +
  
  # 添加相关系数和p值标注（合并到一行）
  geom_text(
    data = correlations,
    aes(
      x = max(clean_data_long$meanlogFC, na.rm = TRUE) * 0.8,
      y = max(clean_data_long$quartet_value, na.rm = TRUE) * (0.9 - 0.15 * (as.numeric(factor(sample_type)) - 1)),
      label = label_sci,  # 使用合并的标签
      color = sample_type
    ),
    hjust = 1, 
    vjust = 0.5, 
    size = 3.2, 
    fontface = "bold", 
    show.legend = FALSE
  ) +
  
  # Nature风格的颜色配置
  scale_color_manual(
    name = "Quartet Samples",
    values = c("D5/D6" = "#1F77B4", "F7/D6" = "#FF7F0E", "M8/D6" = "#2CA02C"),
    labels = c("D5 vs D6", "F7 vs D6", "M8 vs D6")
  ) +
  scale_fill_manual(
    name = "Quartet Samples",
    values = c("D5/D6" = "#1F77B4", "F7/D6" = "#FF7F0E", "M8/D6" = "#2CA02C"),
    labels = c("D5 vs D6", "F7 vs D6", "M8 vs D6")
  ) +
  
  # 坐标轴标签和标题
  labs(
    title = "Correlation between HCC1395 and Quartet Samples",
    x = "HCC1395/BL log₂FC",
    y = "Quartet log₂FC",
    color = "Quartet Samples",
    fill = "Quartet Samples"
  ) +
  my_theme

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
common_gene.maqc.vs.hcc.df <- merge(maqc.ref[,c('gene','log2FC')],
                                    subset(hcc.ref,hcc.ref$gene %in% common_gene.maqc.vs.hcc)[,c('gene','meanlogFC')]) %>% unique()
write.csv(common_gene.maqc.vs.hcc.df,file = "HCC1395.vs.Quartet/common_gene.maqc.vs.hcc.df.csv")
df <- common_gene.maqc.vs.hcc.df
df_filtered <- na.omit(df[,c('DEG','final')])
colnames(df_filtered) <- c('A','B')


#logFC 分布--------------

# 处理 Quartet 数据
df_maqc <- data.frame(
  meanlogFC = maqc.ref$log2FC,
  Source = "MAQC",
  CompareGroup = 'A/B' # 这里包含三个比较组
)


# 合并数据
plot_data.maqc <- rbind(df_hcc, df_maqc)

# 3. 绘图
p_logFC.maqc <- ggplot(plot_data.maqc, aes(x = meanlogFC, fill = CompareGroup, pattern = Source)) +
  # 使用密度图展现分布，也可以换成 geom_histogram_pattern
  geom_density(alpha = 0.5) +
  # 设置颜色：你可以根据需要自定义四个组的颜色
  scale_fill_manual(values = compare_colors) +
  labs(
    title = "Distribution of meanlogFC",
    x = "meanlogFC",
    y = "Density",
    fill = "Comparison Group",
    pattern = "Data Source"
  ) +
  my_theme

p_compare_dist <- (p_logFC+p_logFC.maqc)+plot_layout(widths = c(1,1),guides = 'collect') & theme(legend.position = 'bottom')

ggsave("HCC1395.vs.Quartet/HCC1395.vs.Quartet.MAQC.distribution.pdf",p_compare_dist,width=8.15, height=5.7)


# 数据预处理 - 转换为长格式

colnames(common_gene.maqc.vs.hcc.df) <- c('gene','maqc_logFC','hcc1395_logFC')
# 计算相关系数
cor_test_result <- cor.test(
  common_gene.maqc.vs.hcc.df$hcc1395_logFC, 
  common_gene.maqc.vs.hcc.df$maqc_logFC,
  method = "pearson",  # 使用Pearson相关系数
  use = "complete.obs"  # 仅使用完整观测值
)

# 创建相关系数和p值的数据框
correlations.maqc <- data.frame(
  correlation = cor_test_result$estimate,
  p_value = cor_test_result$p.value,
  method = "Pearson"
)

# 创建Nature风格的单一散点图
p_hcc_maqc <- ggplot(common_gene.maqc.vs.hcc.df, aes(x = hcc1395_logFC, y = maqc_logFC)) +
  # 添加散点
  geom_point(alpha = 0.8, size = 1.5,shape = 21,fill="#1F77B4",
             color = "white",
             stroke = 1) +
  # 为每个样本类型添加回归线
  geom_smooth(method = "lm", se = TRUE, fill = "#1F77B4", alpha = 0.1) +
  # 添加相关系数标注
  # 添加相关系数和显著性标记
  annotate(
    "text",
    x = min(common_gene.maqc.vs.hcc.df$hcc1395_logFC, na.rm = TRUE) * 0.9,
    y = max(common_gene.maqc.vs.hcc.df$maqc_logFC, na.rm = TRUE) * 0.95,
    label = sprintf("r = %.3f\np = %.2e", 
                    correlations.maqc$correlation, 
                    correlations.maqc$p_value),
    color = "black",
    hjust = 0,
    vjust = 1,
    size = 4,
    fontface = "bold"
  ) +
  
  # 根据p值添加显著性星号标记
  annotate(
    "text",
    x = min(common_gene.maqc.vs.hcc.df$hcc1395_logFC, na.rm = TRUE) * 0.9,
    y = max(common_gene.maqc.vs.hcc.df$maqc_logFC, na.rm = TRUE) * 0.85,
    label = ifelse(correlations.maqc$p_value < 0.001, "***",
                   ifelse(correlations.maqc$p_value < 0.01, "**",
                          ifelse(correlations.maqc$p_value < 0.05, "*", "ns"))),
    color = ifelse(correlations.maqc$p_value < 0.05, "red", "gray"),
    size = 6,
    fontface = "bold",
    hjust = 0
  ) +
  # 坐标轴标签和标题
  labs(
    title = "",
    x = "HCC1395/BL logFC",
    y = "MAQC logFC"
  ) +my_theme


p_compare_dot <- (p_hcc_quartet +p_hcc_maqc)+plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

ggsave("HCC1395.vs.Quartet/HCC1395.vs.Quartet.MAQC.scatter.pdf",p_compare_dot,width=8.15, height=5.7)


if(F){
  
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
}
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


#三种标准物质的基因集注释--------

get_annotation <- function(result,gene_df,ref_type){
  core_genes <- result$gene_name %>% unique()
  background_genes <- gene_df$gene_name
  ego <- enrichGO(gene = core_genes, # 你的核心基因列表
                  universe = background_genes, # 背景基因集
                  OrgDb = org.Hs.eg.db,
                  keyType = 'SYMBOL', # 如果core_genes是Symbol
                  ont = "ALL", # 可选 "BP", "MF", "CC"
                  pAdjustMethod = "BH", # FDR校正
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE) # 将ENTREZID转成Symbol显示
  
  #write.csv(ego@result,file = paste0('genome_and_transcriptom_combine.GO.',ref_type,'.csv'),row.names = F)
  core_entrez <- bitr(core_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  bg_entrez <- bitr(background_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  #BiocManager::install("KEGG.db")
  kk <- enrichKEGG(gene = core_entrez,
                   organism = 'hsa', # 人类
                   universe = bg_entrez,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.2)
  
  return(list(ego,kk))
}

hcc.ref <- merge(hcc.ref,gene_df,by='gene')
hcc.ref.anno <- get_annotation(hcc.ref,gene_df)
hcc.ref.GO <- hcc.ref.anno[[1]]@result
hcc.ref.GO$Ref_type = 'HCC1395'
hcc.ref.KEGG <-  hcc.ref.anno[[2]]@result
hcc.ref.KEGG$Ref_type = 'HCC1395'

quartet.ref$gene_name <- quartet.ref$Gene_Name
quartet.ref.anno <- get_annotation(quartet.ref,gene_df)
quartet.ref.GO <- quartet.ref.anno[[1]]@result
quartet.ref.GO$Ref_type = 'Quartet'
quartet.ref.KEGG <-  quartet.ref.anno[[2]]@result
quartet.ref.KEGG$Ref_type = 'Quartet'

maqc.ref$gene_name <- maqc.ref$Gene.Name
maqc.ref.anno <- get_annotation(maqc.ref,gene_df)
maqc.ref.GO <- maqc.ref.anno[[1]]@result
maqc.ref.GO$Ref_type = 'MAQC'
maqc.ref.KEGG <-  maqc.ref.anno[[2]]@result
maqc.ref.KEGG$Ref_type = 'MAQC'

library(tidyverse)
library(pheatmap)
library(pheatmap)
library(viridis)
library(RColorBrewer)

aim_cols = c('ID','Description','p.adjust','Ref_type')
combined_go <- bind_rows(hcc.ref.GO %>%arrange(p.adjust) %>% head(15) %>% select(aim_cols), 
                           quartet.ref.GO %>%arrange(p.adjust) %>% head(15) %>% select(aim_cols), 
                           maqc.ref.GO %>%arrange(p.adjust) %>% head(15) %>% select(aim_cols))

combined_kegg <- bind_rows(hcc.ref.KEGG %>%arrange(p.adjust) %>% head(15) %>% select(aim_cols), 
                         quartet.ref.KEGG %>%arrange(p.adjust) %>% head(15) %>% select(aim_cols), 
                         maqc.ref.KEGG %>%arrange(p.adjust) %>% head(15) %>% select(aim_cols))


get_annotation_plot <- function(combined_kegg,anno_type){
  plot_data <- combined_kegg %>%
    select(Description, Ref_type, p.adjust) %>%
    mutate(
      `-log10(p.adjust)` = -log10(p.adjust)
    ) %>%
    # 去除可能存在的无限大值（由于p值为0导致）
    filter(is.finite(`-log10(p.adjust)`))
  
  # 5. 将长数据转换为宽数据，以便构建矩阵
  heatmap_matrix <- plot_data %>%
    select(Description, Ref_type, `-log10(p.adjust)`) %>%
    pivot_wider(
      names_from = Ref_type,
      values_from = `-log10(p.adjust)`,
      # 对于某些通路在某个组别可能不存在的情况，用0或NA填充
      values_fill = list(`-log10(p.adjust)` = 0)
    ) %>%
    # 将通路名设置为行名
    column_to_rownames(var = "Description") %>%
    # 转换为矩阵
    as.matrix()
  
  # 查看矩阵前几行
  head(heatmap_matrix)
  
  
  significant_pathways <- combined_kegg %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%  # 按p值排序
    pull(Description) %>%
    unique()
  
  # 从矩阵中筛选显著通路
  heatmap_matrix_sig <- heatmap_matrix[rownames(heatmap_matrix) %in% significant_pathways, ]
  
  # 确保列顺序正确
  column_order <- c("HCC1395", "Quartet", "MAQC")
  heatmap_matrix_sig <- heatmap_matrix_sig[, column_order]
  
  # 创建列注释
  annotation_col <- data.frame(
    CellLine = factor(colnames(heatmap_matrix_sig), levels = colnames(heatmap_matrix_sig))
  )
  rownames(annotation_col) <- colnames(heatmap_matrix_sig)
  
  # 定义Nature风格的配色方案
  # 使用viridis配色方案（Nature常用）
  nature_colors <- viridis(100)
  # 或者使用蓝-白-红方案但调整色调
  # nature_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  
  # 设置行注释颜色（如果需要）
  ann_colors <- list(
    CellLine = c(
      HCC1395 = "#1F77B4", 
      Quartet = "#FF7F0E", 
      MAQC = "#2CA02C"
    )
  )
  
  # 优化行名显示：缩短过长的通路名称
  shorten_pathway_names <- function(names) {
    # 移除"hsa"前缀和数字
    short_names <- gsub("hsa[0-9]+: ", "", names)
    # 进一步简化过长的名称
    short_names <- ifelse(nchar(short_names) > 40, 
                          paste0(substr(short_names, 1, 40), "..."), 
                          short_names)
    return(short_names)
  }
  
  # 应用缩短的名称
  short_rownames <- shorten_pathway_names(rownames(heatmap_matrix_sig))
  
  # 绘制Nature风格热图
  p <- pheatmap(
    mat = heatmap_matrix_sig,
    color = nature_colors,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    scale = "row",
    border_color = "grey80",  # 使用淡灰色边框区分单元格
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 7,         # 稍微减小行字体大小
    fontsize_col = 10,
    fontsize = 9,             # 设置基本字体大小
    angle_col = 45,           # 倾斜列名以便更好地显示
    main = anno_type,
    legend = TRUE,
    annotation_legend = TRUE,
    annotation_names_col = TRUE,
    labels_row = short_rownames,  # 使用缩短的行名
    gaps_col = c(1, 2),       # 在列之间添加间隔
    treeheight_row = 20,      # 减小行聚类树的高度
    treeheight_col = 0,       # 隐藏列聚类树
    silent = FALSE            # 显示绘图
  )
  return(p)
}

go_p <- get_annotation_plot(combined_go,'GO')
kegg_p <- get_annotation_plot(combined_kegg,'KEGG')

ggsave('HCC1395.vs.Quartet/GO_annotation_top15.pdf',go_p,height =8.15 ,width =5.7*1.5,dpi = 300)
ggsave('HCC1395.vs.Quartet/KEGG_annotation_top15.pdf',kegg_p,height =8.15 ,width =5.7*1.5,dpi = 300)

#挑选HCC1395细胞系的原因-------
1. 
