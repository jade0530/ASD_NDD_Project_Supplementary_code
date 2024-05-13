library(tidyr)
library(dplyr)
library(readr)
library(ggplotify)
library(ggpubr)

# Create a separate legend with an empty graph for svg generation
# Dup/Del Signals
svg("legend.svg")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = "top",   fill = c("blue3", "#FFCC00", "#333333", "red2"), cex=1.2, legend = c("ASD", "ASD_ID_DD","ID_DD","SCZ"), box.lty=0)
dev.off()

# Hi-C Interaction 
# Enhacer HMM signals 
svg("legend_HMM.svg")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = "top",   fill = c("#ff0000", "#ff4500", "#32cd32", "#008000", "#006400", 
                             "#c2e105", "#ffff00", "#66cdaa", "#8a91d0", "#cd5c5c",
                             "#e9967a", "#bdb76b", "#808080", "#c0c0c0", "#ffffff"), cex=1.2, legend = c("Active TSS", "Flanking Active TSS, Transcr at gene 5' and 3'", "Strong transcription",
                                                                                                         "Weak transcription", "Genic enhancers", "Enhancers", "ZNF genes & repeats", "Heterochromatin"
                                                                                                         , 
                                                                                                         "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enh", "Bivalent Enhancer", "Repressed PolyComb", 
                                                                                                         "Weak Repressed PolyComb", "Quiescent/Low"), box.lty=0)
dev.off()

# Literature 
svg("interaction_legend.svg")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = "top",   fill = c("#10aeff", "#c210b0", "#1018c2"), cex=1.2, legend = c("Synapse NeuN samples", "PMID:30555922","PMID:31367015"), box.lty=0)
dev.off()

# Case/Control
svg("cnv_legend.svg")
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend(x = "top",   fill = c("#270f5c", "#d0021bff", "#5f8000"), cex=1.2, legend = c("Duplicated Case", "Deleted Case", "Control"), box.lty=0)
dev.off()
