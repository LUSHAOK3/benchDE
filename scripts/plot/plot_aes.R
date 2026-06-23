library(ggplot2)

nSample <- c(3,5,10,30,50)
plot.format=theme(plot.background=element_blank(),
                  panel.grid=element_blank(),
                  panel.background=element_blank(),
                  panel.border=element_rect(color="black",linewidth=0.5,fill=NA),
                  axis.line=element_blank(),
                  strip.text = element_text(color="black",size=7),
                  axis.ticks=element_line(color="black",linewidth=0.5),
                  axis.text=element_text(color="black",size=7),
                  axis.title=element_text(color="black",size=7),
                  plot.title=element_text(color="black",size=7),
                  legend.background=element_blank(),
                  legend.key=element_blank(),
                  legend.text=element_text(color="black",size=7),
                  legend.title=element_text(color="black",size=7))


cols <- c("#CAB2D6","#6A3D9A","#A6CEE3","#56b4e9","#e69f00","#11264f",
"#009e73","#d55e00","#FB9A99","#cc79a7","#E31A1C",
"#666666","#FFFF33","#FDBF6F","#FF7F00","#4DAF4A","#662506","#B15928")

methods <- c('DESeq2','DESeq','edgeR.lrt','edgeR.qlf','DSS','voom','NBPSeq','Wilcoxon',
          'T.test','NOISeq','ROTS','ABSSeq')
