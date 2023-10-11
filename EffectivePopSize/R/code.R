################################################################
#Author: Josephine Johnson
#
#


library(gdsfmt)
library(SNPRelate)
library(dplyr)
library(magrittr)
library(tidyr)
# Get the path of the plink files
bed <- "/Users/josephine.p.johnson/OneDrive - North Dakota University System/Documents/PhD/chapter1/LD_Ne/LD/vcf/plink.bed"
fam <- "/Users/josephine.p.johnson/OneDrive - North Dakota University System/Documents/PhD/chapter1/LD_Ne/LD/vcf/plink.fam"
bim <- "/Users/josephine.p.johnson/OneDrive - North Dakota University System/Documents/PhD/chapter1/LD_Ne/LD/vcf/plink.bim"
# Build the gdsdb
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn,out.gdsfn ="test.gds")
snpgdsSummary("test.gds")
genofile <- snpgdsOpen("test.gds")
# Getting the PCA object
pca <- snpgdsPCA(genofile, num.thread=2,autosome.only=FALSE)
# Calculate the percentage of variance explained by the first 10 PC
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
pc<-round(pc.percent)

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
pc.percent[1]
dim(pc.percent)
## Coloring by Batches

```{r, fig.width=10, fig.height=4}
table <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  EV5 = pca$eigenvect[,5],
                  stringsAsFactors = FALSE)
colnames(tab) <- c("ID", "EV1", "EV2", "EV3", "EV4", "EV5")
write.table(tab, file = "./Eigenvalues.txt", row.names = F, col.names = F, sep = "\t", quote = F)
#Load the metadata
metadata <- read.delim("/home/roberto/Sorghum/WGS/VCF_2.0/PCA/TERRA_metadata.txt")
#Join the eigenvalues with the metadata
taba <- tab %>%
  left_join(metadata, by="ID")
library(ggplot2)
library(wesanderson)
library(patchwork)
fill <- wes_palette("IsleofDogs1")
a <- ggplot(data = tab) +
  geom_point(mapping = aes(x = EV1,y = EV2), color = 'darkblue', show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA for SNPS in Depth 4",
  x = paste0("Principle Component 1 (",round(pc.percent[1])," %)"),
  y = paste0("Principle Component 2 (", round(pc.percent[2])," %)")) +
  theme_minimal()
p2 <- ggplot(taba, aes(EV1,EV3)) +
  geom_point(aes(colour = factor(Batch)), alpha =0.5, shape=19) +
  scale_color_manual(values=c("red", "blue", fill[1], fill[2], fill[3], fill[4])) +
  ylab("Eigenvalue 3 (3.29%)") +
  xlab("Eigenvalue 1 (4.99%)") +
  theme(legend.position="right", legend.direction="vertical",
        legend.title = element_blank()) +
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
S5 <- p1 + p2
ggsave("/workdir/Sorghum/Figures/FigS5A.png", plot=S5, device = "png", scale = 1, width = 10, height = 4, units = c("in"),dpi = 200, limitsize = TRUE)
```

## Coloring by Type (Sweet, cellulosic, grain)

```{r, fig.width=15, fig.height=6}
p <- ggplot(taba, aes(EV1, EV2)) +
  geom_point(aes(colour = factor(Type)), alpha =0.5, shape=19) +
  scale_color_manual(values=c("red", "blue", fill[1], fill[2], fill[3], fill[4])) +
  ylab("Eigenvalue 2 (3.90%)") +
  xlab("Eigenvalue 1 (4.99%)") +
  theme(legend.position="none", legend.direction="vertical",
        legend.title = element_blank()) +
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
p2 <- ggplot(taba, aes(EV1, EV3)) +
  geom_point(aes(colour = factor(Type)), alpha =0.5, shape=19) +
  scale_color_manual(values=c("red", "blue", fill[1], fill[2], fill[3], fill[4])) +
  ylab("Eigenvalue 3 (3.20%)") +
  xlab("Eigenvalue 1 (4.99%)") +
  theme(legend.position="right", legend.direction="vertical",
        legend.title = element_blank()) +
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
S5B <- p + p2
ggsave("/workdir/Sorghum/Figures/FigS5B.png", plot=S5B, device = "png", scale = 1, width = 10, height = 4, units = c("in"),dpi = 200, limitsize = TRUE)
```


## Coloring by races (All of them)

```{r, fig.width=15, fig.height=8}
pal <- wes_palette("Zissou1", 100, type = "continuous")
cal <- wes_palette("Cavalcanti1", 20, type = "continuous")
il <- wes_palette("GrandBudapest1", 20, type = "continuous")
p <- ggplot(taba, aes(EV1, EV2)) +
  geom_point(aes(colour = factor(Race)), alpha =0.6, shape=19) +
  scale_color_manual(values=c("red", pal[1], pal[10], "grey", cal[1], cal[2], cal[3], cal[20], cal[19], cal[18], cal[17], cal[16], il[7], il[5], il[3], "blue", "black", "white")) +
  ylab("Eigenvalue 2 (3.90%)") +
  xlab("Eigenvalue 1 (4.99%)") +
  theme(legend.position="none", legend.direction="vertical",
        legend.title = element_blank()) +
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
p2 <- ggplot(taba, aes(EV1, EV3)) +
  geom_point(aes(colour = factor(Race)), alpha =0.6, shape=19) +
  scale_color_manual(values=c("red", pal[1], pal[10], "grey", cal[1], cal[2], cal[3], cal[20], cal[19], cal[18], cal[17], cal[16], il[7], il[5], il[3], "blue", "black", "white")) +
  ylab("Eigenvalue 3 (3.20%)") +
  xlab("Eigenvalue 1 (4.99%)") +
  theme(legend.position="none", legend.direction="vertical",
        legend.title = element_blank()) +
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
p3 <- ggplot(taba, aes(EV1, EV4)) +
  geom_point(aes(colour = factor(Race)), alpha =0.6, shape=19) +
  scale_color_manual(values=c("red", pal[1], pal[10], "grey", cal[1], cal[2], cal[3], cal[20], cal[19], cal[18], cal[17], cal[16], il[7], il[5], il[3], "blue", "black", "white")) +
  ylab("Eigenvalue 4 (2.19%)") +
  xlab("Eigenvalue 1 (4.99%)") +
  theme(legend.position="none", legend.direction="vertical",
        legend.title = element_blank()) +
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
p4 <- ggplot(taba, aes(EV2, EV3)) +
  geom_point(aes(colour = factor(Race)), alpha =0.6, shape=19) +
  scale_color_manual(values=c("red", pal[1], pal[10], "grey", cal[1], cal[2], cal[3], cal[20], cal[19], cal[18], cal[17], cal[16], il[7], il[5], il[3], "blue", "black", "white")) +
  ylab("Eigenvalue 3 (3.20%)") +
  xlab("Eigenvalue 2 (3.90%)") +
  theme(legend.position="right", legend.direction="vertical",
        legend.title = element_blank()) +
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
Allraces <- p+ p2 + p3 + p4 + ncol(2)
ggsave("/workdir/Sorghum/Figures/FigS4A.png", plot=Allraces, device = "png", scale = 1, width = 15, height = 8, units = c("in"),dpi = 200, limitsize = TRUE)
```

## Figure 1A
Includes only the four races, complex and weedy types.

```{r, fig.height=3, fig.width=3.42}
wes <- wes_palette("Cavalcanti1")
fox <- wes_palette("FantasticFox1")
razas = c("Kafir", "Guinea", "Caudatum", "Durra", "Complex", "Weedy")
PCtable <-taba[which(taba$Race %in% razas),]
p <- ggplot(PCtable, aes(EV1, EV2)) +
  geom_point(aes(colour = factor(Race)), alpha =0.6, shape=19, size =1.3) +
  scale_color_manual(values=c(fox[3], "grey", wes[1], wes[2],"red", "black")) +
  ylab("PC 2 (3.90%)") +
  xlab("PC 1 (4.99%)") +
  theme(legend.position=c(0.89, 0.28), legend.direction="vertical",
        legend.title = element_blank(), legend.text=element_text(size=6),
        legend.key.height = unit(0.5, "cm"), legend.key = element_rect(colour = NA, fill = NA)) +
  #legend.text=element_text(size=X)
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "white"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
p
ggsave("/home/roberto/Sorghum/Figures/PCA_pres.tiff", plot= p, device = "tiff", scale = 1, width = 4.42, height = 4, units = c("in"),dpi = 500, limitsize = TRUE)
```

## Fig1A - 3D animated version

```{r}
set.seed(417)
library(plotly)
temp <- rnorm(100, mean=30, sd=5)
pressure <- rnorm(100)
dtime <- 1:100
wes <- wes_palette("Cavalcanti1")
fox <- wes_palette("FantasticFox1")
new <- wes_palette("IsleofDogs1")
plot_ly(x=temp, y=pressure, z=dtime, type="scatter3d", mode="markers", color=temp)
plot_ly(x=PCtable$EV1, y=PCtable$EV2, z=PCtable$EV4, type="scatter3d", mode="markers", color=as.character(PCtable$Race), colors = c(fox[3], new[6], wes[1], wes[2], fox[5], new[4]), size = 1.95  )
```



## Coloring by geographic Location


```{r, fig.width=15, fig.height=4}
p <- ggplot(taba, aes(EV1, EV2)) +
  geom_point(aes(colour = factor(Country_group)), alpha =0.6, shape=19) +
  #scale_color_manual(values=c("red", pal[1], pal[10], "grey", cal[1], cal[2], cal[3], cal[20], cal[19], cal[18], cal[17], cal[16], il[7], il[5], il[3], "blue", "black", "white")) +
  ylab("Eigenvalue 2 (3.90%)") +
  xlab("Eigenvalue 1 (4.99%)") +
  theme(legend.position="none", legend.direction="vertical",
        legend.title = element_blank()) +
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
p2 <- ggplot(taba, aes(EV1, EV4)) +
  geom_point(aes(colour = factor(Country_group)), alpha =0.6, shape=19) +
  #scale_color_manual(values=c("red", pal[1], pal[10], "grey", cal[1], cal[2], cal[3], cal[20], cal[19], cal[18], cal[17], cal[16], il[7], il[5], il[3], "blue", "black", "white")) +
  ylab("Eigenvalue 4 (2.19%)") +
  xlab("Eigenvalue 1 (4.99%)") +
  theme(legend.position="right", legend.direction="vertical",
        legend.title = element_blank()) +
  theme(axis.line.x = element_line(size=1, colour = "black"),
        axis.line.y = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text=element_text(family="Tahoma"),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        legend.key=element_rect(fill="white", colour="white"))
Geo <- p + p2
ggsave("/workdir/Sorghum/Figures/FigS4B.png", plot=Geo, device = "png", scale = 1, width = 15, height = 4, units = c("in"),dpi = 200, limitsize = TRUE)
