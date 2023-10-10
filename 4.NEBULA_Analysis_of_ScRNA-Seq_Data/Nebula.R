# install.packages("devtools")
# library(devtools)
# install_github("lhe17/nebula")
# install.packages('randomcoloR')


library(nebula)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(showtext)
library(scales)
library(randomcoloR)


args<-commandArgs(TRUE)
work_dir <- args[1]

setwd(work_dir)
dir.create("outputs")

file_list <- list.files(pattern = ".csv$")
cluster <- unique(unlist(lapply(file_list,function(x){
  strsplit(x,'_')[[1]][1]
})))

print(cluster)
# ct <- 'Fibroblasts'
for (ct in cluster) {
  
  counts <- read.csv(paste0(ct,'_count.csv'), row.names = 1)
  count <- as.matrix(counts)
  batch <- read.csv(paste0(ct,'_batch.csv'), header = F,row.names = 1)
  label <- read.csv(paste0(ct,'_label.csv'), header = T,row.names = 1)
  output <- data.frame()
  
  for (mb in colnames(label)) {
    
    cat(ct,'-',mb,'\n')
    
    df = label[mb]
    colnames(df) <- 'V1'
    df = model.matrix(~V1 , data=df)
    
    re = nebula(count, batch$V2, pred=df)
    
    result <- re$summary
    result <- result[,c(8,2,6)]
    
    colnames(result) <- c('gene','logFC','pVal')
    rownames(result) <- result$gene
    
    result$adj.P <- p.adjust(result$pVal,method = 'BH')
    result$'-Log10(adj.P)' <- -log10(result$adj.P)
    
    result <- result[order(result$logFC,decreasing = T),]
    
    # 设定上调下调
    result$expression <- ifelse(result$logFC >= 1 & result$adj.P < 0.05,"Up-regulated",
                                ifelse(result$logFC <= -1 & result$adj.P < 0.05,"Down-regulated","NS."))
    
    # 删除logFC很高，但NS的代谢基�?    
    result <- result %>%
      filter(!(abs(result$logFC) > 10 & result$adj.P > 0.05))
    
    # 加一列标签，用以区分微生�?    
    result$label <- mb
    output <- rbind(output, result)
  }
  # 保存文件
  write.csv(output,paste0('outputs/',ct,'_DEG_among_microbe.csv'), quote = F) 
 
  plot_multiple_volcano <- function(output_sig,ct){
  # define a color palette
    dfcol<-data.frame(x=c(1:length(unique(output_sig$label))),
                      y=0,
                      label = unique(output_sig$label))
    mycol <- distinctColorPalette(length(unique(output_sig$label)))
    
    top10_sig <- data.frame()
    for (i in unique(output_sig$label)) {
      output_sub <- output_sig[output_sig$label==i,]
      top10 <- output_sub %>% top_n(10,abs(logFC))
      top10_sig <- rbind(top10_sig,top10)
    }
    
    # 多组火山�?    
    p1 <- ggplot()+
    geom_jitter(data = output_sig,
                aes(x = label, y = logFC, color = expression),
                size = 0.6,
                width = 0.3) +
    scale_y_continuous(limits = c(10,1e18),expand=c(0,0))+
    scale_color_manual(values = c("Down-regulated"='blue',"Up-regulated"='red')) +
    geom_text_repel(
      data=top10_sig,
      aes(x=label,y=logFC,label=gene),
      size=3
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 13,
                                color = "black",
                                face = "bold"),
      axis.line.y.left = element_line(color = "black",
                                      size = 1.2),
      axis.line.y = element_blank(),
      # axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1,0),
      legend.text = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      # axis.text.y = element_blank(),
    ) 

    p1

    p2 <- ggplot()+
  
      geom_jitter(data = output_sig,
                  aes(x = label, y = logFC, color = expression),
                  size = 0.6,
                  width =0.3) +
      geom_tile(data = dfcol,
                aes(x=x,y=y),
                height=3,
                color = "black",
                fill = mycol,
                alpha = 1,
                show.legend = F)+
      # labs(x="Cluster")+
      geom_text(data=dfcol,
                aes(x=x,y=y,label=label),
                size = 5,
                color ="white")+
      scale_y_continuous(limits = c(-10,10),expand=c(0,0))+
      scale_color_manual(values = c("Down-regulated"='blue',"Up-regulated"='red'))+
      geom_text_repel(
        data=top10_sig,
        aes(x=label,y=logFC,label=gene),
        force = 2,
        # max.overlaps = 10,
        size=3,
      )+
      theme_minimal()+
      theme(
        axis.title = element_text(size = 13,
                                  color = "black",
                                  face = "bold"),
        axis.line.y.left = element_line(color = "black",
                                        size = 1.2),
        axis.line.y = element_blank(),
        # axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
      )
    
    p2
    
    p3 <- ggplot()+
      # ylim(-1e18,-1e1) +
      scale_y_continuous(limits = c(-1e18,-10),expand=c(0,0))+
      scale_color_manual(values = c("Down-regulated"='blue',"Up-regulated"='red'))+
      geom_text_repel(
        data=top10_sig,
        aes(x=label,y=logFC,label=gene),
        force = 1.2,
        size=3,
      )+
      geom_jitter(data = output_sig,
                  aes(x = label, y = logFC, color = expression),
                  size = 0.6,
                  width = 0.3) +
      theme_minimal()+
      theme(
        axis.title = element_text(size = 13,
                                  color = "black",
                                  face = "bold"),
        axis.line.y.left = element_line(color = "black",
                                        size = 1.2),
        axis.line.y = element_blank(),
        # axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 10),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.text.y = element_blank()
      )
    p3

    p4 <- ggarrange(p1,p2,p3,heights=c(1/5, 3/5,1/5),ncol = 1, nrow = 3,common.legend = TRUE,legend="right",align = "v")
    p4
    print("p4 is or not null")
    is.null(p4)
    ggsave(paste0('outputs/',ct,"_multiple_volcano.pdf"),plot = p4,width = 5+2*length(unique(output_sig$label)),height=15,limitsize = F,device="pdf")
    ggsave(paste0('outputs/',ct,"_multiple_volcano.png"),plot = p4,width = 5+2*length(unique(output_sig$label)),height=15,limitsize = F,device="png")
  }

  output_sig <- output[output$expression != 'NS.',] # 删选显著的差异基因
  # 因子化microbiome列，便于排序
  output_sig$label <- factor(output_sig$label,levels = unique(output_sig$label))
  print(head(output_sig))
  if (nrow(output_sig)>1){
    print("start plot")
    plot_multiple_volcano(output_sig,ct)
  }
}
