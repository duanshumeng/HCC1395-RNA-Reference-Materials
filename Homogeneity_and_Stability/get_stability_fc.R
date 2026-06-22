library(stringr)

#short-time stability
get_deg_ref <- function(raw_count,name_replace=NULL){
  df_count <- read.csv(raw_count)
  if (!is.null(name_replace)){
    name_df <- read.csv(name_replace,header = T)
    colnames(df_count) <- lapply(colnames(df_count),function(x){
      if (x %in% name_df$raw_name){
        return(subset(name_df,raw_name==x)$new_name)
      } else {
        return(x)
      }
      
    })
  }
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
  
  mt_df <- data.frame(library=colnames(df_count),
                      sample=str_split(colnames(df_count),'\\_',simplify = T)[,2],
                      protocol='PolyA',
                      platform='T7',
                      lab='MGI',
                      batch=str_split(colnames(df_count),'\\.',simplify = T)[,2],
                      rep=str_split(colnames(df_count),"_|\\.",simplify = T)[,3],
                      sample_rep=paste0(str_split(colnames(df_count),'\\_',simplify = T)[,2],'_',
                                        str_split(colnames(df_count),"_|\\.",simplify = T)[,3]))
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
  
  if (F){
    df_count.dete$gene <- rownames(df_count.dete)
    df_count.dete.melt <- df_count.dete %>% melt()
    df_count.dete.melt <- merge(df_count.dete.melt,mt_df[,c('library','sample','sample_rep','batch')],by.x = 'variable',by.y = 'library') %>% unique()
    
    
    #3. 计算差异表达基因
    result <- data.frame()
    df <- df_count.dete.melt
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
      df_merged$log2fc <- log2(df_merged$value_HCC1395 / df_merged$value_HCC1395BL+0.01)
      
      # 把结果加入到结果数据框中
      result <- rbind(result, df_merged)
    }
    
    
    result.sub <- result[,c('gene','value_HCC1395','value_HCC1395BL','batch_HCC1395','log2fc')] %>% unique()
    colnames(result.sub) <- c('gene','HCC1395','HCC1395BL','batch','logfc')
    write.csv(result.sub,file=gsub('.csv','.log2fc.csv',raw_count),row.names = F)
    
  }
  system(paste0('Rscript /Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Reference_standard_protocol-main/scripts/deg_limma-hcc1395.R',
                ' -i ',sub('.csv','.detectable.csv',paste0(name,'.multibatch_count.csv')),' -m ','hcc1395_metadata.csv',' -o ','./'))

  return(paste0(name,' Finished=========================='))
}
#Long stability-------

#KD stability----------
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/KD_stability')
get_deg_ref("FeatureCounts.gene_quant-KD.csv")
get_deg_ref("HTseq.gene_quant-KD.csv")
get_deg_ref("Salmon.gene_quant-KD.csv")
get_deg_ref("Stringtie.gene_quant-KD.csv")

#DW4°C stability----------------
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/DW4_stability')
get_deg_ref("FeatureCounts.gene_quant-DW4.csv","DW4_name_replace.csv")
get_deg_ref("HTseq.gene_quant-DW4.csv","DW4_name_replace.csv")
get_deg_ref("Salmon.gene_quant-DW4.csv","DW4_name_replace.csv")
get_deg_ref("Stringtie.gene_quant-DW4.csv","DW4_name_replace.csv")

#DW25°C stability----------------
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/DW25_stability')
get_deg_ref("FeatureCounts.gene_quant.csv","DW25_name_replace.csv")
get_deg_ref("HTseq.gene_quant.csv","DW25_name_replace.csv")
get_deg_ref("Salmon.gene_quant.csv","DW25_name_replace.csv")
get_deg_ref("Stringtie.gene_quant.csv","DW25_name_replace.csv")

#DW-CO2 stability----------------
setwd('/Users/duanshumeng/生物信息/PGx_lab/HCC1395/RNAseq/Data_analysis/Homogeneity_and_Stability/DWCO2_stability')
get_deg_ref("FeatureCounts.gene_quant.csv","DWCO2_name_replace.csv")
get_deg_ref("HTseq.gene_quant.csv","DWCO2_name_replace.csv")
get_deg_ref("Salmon.gene_quant.csv","DWCO2_name_replace.csv")
get_deg_ref("Stringtie.gene_quant.csv","DWCO2_name_replace.csv")




