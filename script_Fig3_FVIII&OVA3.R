library("immunarch")
library("dplyr")
library("tidyverse")
library("ggplot2")
library("circlize")
library("ggseqlogo")

#Tous les 4FVIII commun avec OVA
DT <- repLoad(.path = "/mnt/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/_Analyse_corriges_4_FVIII&OVA")

DT <- repLoad(.path = "/Volumes/TOSHIBAT2/ngsTCRFVIII_Rscript/_Analyse_corriges_4_FVIII&OVA")


DT$data$RUN_1_1179_FVIII$Donor<-"1179"
DT$data$RUN_1_1179_FVIII$Antigen<-"FVIII"

DT$data$RUN_1_1179_OVA$Donor<-"1179"
DT$data$RUN_1_1179_OVA$Antigen<-"OVA"
#--------------------------------------------
DT$data$RUN_1_1195_FVIII$Donor<-"1195"
DT$data$RUN_1_1195_FVIII$Antigen<-"FVIII"

DT$data$RUN_1_1195_OVA$Donor<-"1195"
DT$data$RUN_1_1195_OVA$Antigen<-"OVA"
#-------------------------------------------
DT$data$RUN_1_1207_FVIII$Donor<-"1207"
DT$data$RUN_1_1207_FVIII$Antigen<-"FVIII"

DT$data$RUN_1_1207_OVA$Donor<-"1207"
DT$data$RUN_1_1207_OVA$Antigen<-"OVA"
#----------------------------------------------
DT$data$RUN_3_1180_FVIII$Donor<-"1180"
DT$data$RUN_3_1180_FVIII$Antigen<-"FVIII"

DT$data$RUN_3_1180_OVA$Donor<-"1180"
DT$data$RUN_3_1180_OVA$Antigen<-"OVA"



DT2 <- rbind(DT$data$RUN_1_1179_FVIII,DT$data$RUN_1_1179_OVA,DT$data$RUN_1_1195_FVIII,DT$data$RUN_1_1195_OVA,DT$data$RUN_1_1207_FVIII,DT$data$RUN_1_1207_OVA,DT$data$RUN_3_1180_FVIII,DT$data$RUN_3_1180_OVA)
view(DT2)
---------------------------------------------------------------------------------------------------------
#Figure 3A_FVIII&OVA repertoire distribution 4 donors 

Pp<-ggplot(DT2, aes(factor(Antigen),log10(Proportion)))
Pp+geom_violin(scale = "count",aes(factor(Antigen), fill = factor(Antigen))) + facet_grid(.~Donor)


Pp2<-ggplot(DT2, aes(factor(Antigen),log(Proportion)))
Pp2+geom_violin(scale = "count",aes(factor(Antigen), fill = factor(Antigen))) + facet_grid(.~Donor) + scale_fill_manual(values = c("#009966","#CC3300")) + theme_minimal() + geom_violin(draw_quantiles = c(0.5))

#essai largeur
Pp3<-ggplot(DT2, aes(factor(Antigen),log10(Proportion)))
Pp3+geom_violin(aes(factor(Antigen), fill = factor(Antigen))) + facet_grid(.~Donor) + scale_fill_manual(values = c("#009966","#CC3300")) + theme_minimal() 

#--------------------------------------------------------------
#50% rep 

FVIII_1 <- DT$data$RUN_1_1179_FVIII[1:16,] 
OVA_1 <- DT$data$RUN_1_1179_OVA[1:15,] 

FVIII_2 <- DT$data$RUN_1_1195_FVIII[1:6,] 
OVA_2 <- DT$data$RUN_1_1195_OVA[1:1,] 

FVIII_3 <- DT$data$RUN_1_1207_FVIII[1:15,] 
OVA_3 <- DT$data$RUN_1_1207_OVA[1:107,] 

FVIII_4 <- DT$data$RUN_3_1180_FVIII[1:15,] 
OVA_4 <- DT$data$RUN_3_1180_OVA[1:9,] 

DT50FVIII <- rbind(FVIII_1, FVIII_2,FVIII_3,FVIII_4)
DT50OVA <- rbind(OVA_1, OVA_2, OVA_3, OVA_4)
#-------------------------------------------------
# CD3 lenght on frequency (Clones in the 50% repertoire space)
#FVIII
DT50FVIII_f <- cbind(DT50FVIII, frequence = DT50FVIII$Clones/sum(DT50FVIII$Clones)) 

Lenght <- nchar(DT50FVIII$CDR3.aa)
ggplot(DT50FVIII_f, aes(x= Lenght , y= DT50FVIII_f$frequence))+geom_bar(stat="identity", width=0.7, fill="#009966")+xlim(5,25)+theme_minimal()+labs(subtitle = "Amino acid CDR3 length for FVIII")

#OVA
DT50OVA_f <- cbind(DT50OVA, frequence = DT50OVA$Clones/sum(DT50OVA$Clones)) 

Lenght <- nchar(DT50OVA$CDR3.aa)
ggplot(DT50OVA_f, aes(x= Lenght , y= DT50OVA_f$frequence))+geom_bar(stat="identity", width=0.7, fill="#CC3300")+xlim(5,25)+theme_minimal()+labs(subtitle = "Amino acid CDR3 length for Ovalbumin")

#----------------------------------------------------------------------------------------------------
#CD3 logo

#FVIII

FVIIIalign <- DT50FVIII_f
FVIIIalign$CDR3.aa <- FVIIIalign$V1

kmers_FVIII <- getKmers(FVIII50align, 19)
logoFVIII <- kmer_profile(kmers_FVIII, "prob")
vis(logoFVIII, .plot = "seq")

#OVA

OVA50align <-DT50OVA_f
OVA50align$CDR3.aa <- alignOVA$V1

Kmers_OVA <- getKmers(OVA50align,28 )
logoOVA <- kmer_profile(Kmers_OVA, "prob")
vis(logoOVA, .plot = "seq")
#-----------------------------------------------
#CD3 logo with 'ggseqlogo' try 

#FVIII

FVIII50align <- DT50FVIII_f
FVIII50align$CDR3.aa <- FVIIIalign$V1

#custom height 
aalist <- list('A', 'R', 'N', 'C','D','E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'  )
custom_mat <- matrix( ncol = 19, nrow=20, dimnames=list(c('A', 'R', 'N', 'C','D','E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'  )))
CD3split <- strsplit(FVIIIalign$V1, '')

for (position in (1:19)) {
  for (aa in (1:20)) {
    compteur <- 0
    for (sequence in 1:length(CD3split)) {
      if (CD3split[[sequence]][[position]] == aalist[[aa]]) {
        compteur <- compteur + 1
      }
    }
    custom_mat[aa, position] <- compteur
  }
}

custom_mat <- custom_mat/length(CD3split)


p1 <- ggseqlogo(custom_mat, method='custom', seq_type='aa') ; p1

#--------------------------------

#OVA

OVA50align <-DT50OVA_f
OVA50align$CDR3.aa <- alignOVA$V1

#custom height 
aalist <- list('A', 'R', 'N', 'C','D','E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'  )
custom_mat <- matrix( ncol = 28, nrow=20, dimnames=list(c('A', 'R', 'N', 'C','D','E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'  )))
CD3split <- strsplit(alignOVA$V1, '')

for (position in (1:28)) {
  for (aa in (1:20)) {
    compteur <- 0
    for (sequence in 1:length(CD3split)) {
      if (CD3split[[sequence]][[position]] == aalist[[aa]]) {
        compteur <- compteur + 1
      }
    }
    custom_mat[aa, position] <- compteur
  }
}

custom_mat2 <- custom_mat/length(CD3split)


p2 <- ggseqlogo(custom_mat2, method='custom', seq_type='aa') ; p2
#-------------------------------------------------------------


#VDJ pairing FVIII

#prep matrice for Chord diagram
listVF <- sort(unique(unlist(strsplit(unique(DT50FVIII_f$V.name), ", "))))
listJF <- sort(unique(unlist(strsplit(unique(DT50FVIII_f$J.name), ", "))))
VJFVIII <-matrix(0, length(listVF), length(listJF))
VJFVIII <- data.frame(VJFVIII)
colnames(VJFVIII) <- listJF
rownames(VJFVIII) <- listVF
for (v in (1:nrow(VJFVIII))){
  for (j in (1:ncol(VJFVIII))){
    for (a in (1:nrow(DT50FVIII_f))){
      if ((listVF[[v]] %in% unlist(strsplit(unique(DT50FVIII_f$V.name[a]), ", "))) & (listJF[[j]] %in% unlist(strsplit(unique(DT50FVIII_f$J.name[a]), ", ")))){
        VJFVIII[v, j] <- VJFVIII[v, j] + 1
      }
    }    
  }
}

col = c('TRBV10-3'="#99CC00",'TRBV11-2'="#CCFF00",'TRBV16'="#99FF00",
        'TRBV20-1'="#66CC66", 'TRBV25-1'="#33FF00",'TRBV3-1'="#00CC00",
        'TRBV6-3'="#33CC33",'TRBV7-2'="#006600", 'TRBV7-3'="#339900",
        'TRBV7-4'="#009966",'TRBV7-8'="#33CC00",'TRBV7-9'="#66FF33",
        'TRBV4-1'="#99FF66",'TRBV10-2'="#CCFF66",'TRBV11-3'="#99CC33",
        'TRBV28'="#66CC33", 'TRBV17'="#339900",'TRBV29-1'="#336600",
        'TRBV14'="#00FF00",'TRBV27'="#33FF33", 'TRBV5-4'="#669933",
        'TRBV12-3'="#666633",'TRBV21-1'="#666633",'TRBV24-1'="#666600",
        'TRBV10-1'="#999900", 'TRBV12-1'="#999933",'TRBV12-2'="#999966",
        'TRBV12-4'="#99FF33",'TRBV12-5'="#99CC66",'TRBV18'="#66FF66",
        'TRBV6-5'="#99FF99",'TRBV6-6'="#CCFF99",'TRBV9'="#CCCC99",
        'TRBV19'="#CCCC66",'TRBV13'="#CCCC33",'TRBV15'="#CCCC00",
        'TRBV30'="#99CC99", 'TRBV5-1'="#99FFCC",'TRBV5-2'="#66CC99",
        'TRBV23-1'="#669966",'TRBV5-3'="#339933", 'TRBV7-1'="#339966",
        'TRBV2'="#00CC66",'TRBJ1-1'="#009900", 'TRBJ1-2'="#006633",
        'TRBJ1-3'="#009933",'TRBJ1-4'="#33CC66",'TRBJ1-5'="#33CC99",
        'TRBJ1-6'="#336633",'TRBJ2-1'="#00CC33",'TRBJ2-2'="#CCFFCC", 
        'TRBJ2-3'="#33FF66",'TRBJ2-4'="#66FF99", 'TRBJ2-5'="#00FF66",
        'TRBJ2-6'="#00FF33",'TRBJ2-7'="#006600", 'TRBV4-2'="#99FF33")

VJFVIII <- as.matrix(VJFVIII)
dFVIII <- data.frame(from = rep(rownames(VJFVIII), times = ncol(VJFVIII)),
                   to = rep(colnames(VJFVIII), each = nrow(VJFVIII)),
                   value = as.vector(VJFVIII),
                   stringsAsFactors = FALSE)
chordDiagram(dFVIII,grid.col = col, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(dFVIII))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

#VDJ pairing OVA

#prep matrice for Chord diagram
listVO <- sort(unique(unlist(strsplit(unique(DT50OVA_f$V.name), ", "))))
listJO <-sort(unique(unlist(strsplit(unique(DT50OVA_f$J.name), ", "))))
VJOVA <-matrix(0, length(listVO), length(listJO))
VJOVA <- data.frame(VJOVA)
colnames(VJOVA) <- listJO
rownames(VJOVA) <- listVO
for (v in (1:nrow(VJOVA))){
  for (j in (1:ncol(VJOVA))){
    for (a in (1:nrow(DT50OVA_f))){
      if ((listVO[[v]] %in% unlist(strsplit(unique(DT50OVA_f$V.name[a]), ", "))) & (listJO[[j]] %in% unlist(strsplit(unique(DT50OVA_f$J.name[a]), ", ")))){
        VJOVA[v, j] <- VJOVA[v, j] + 1
      }
    }    
  }
}

col = c('TRBV10-1'="#FF0033",'TRBV10-2'="#FFCCCC",'TRBV11-2'="#FF3300",
        'TRBV11-3'="#990033", 'TRBV12-1'="#FF9966",'TRBV12-2'="#FF6633",
        'TRBV12-3'="#CC3300",'TRBV12-4'="#CC6633", 'TRBV19'="#990000",
        'TRBV20-1'="#CC3333",'TRBV3-1'="#FF3333",'TRBV4-2'="#CC3366",
        'TRBV5-2'="#CC6666",'TRBV5-4'="#FFCCCC",'TRBV6-5'="#FF99CC",
        'TRBV6-6'="#CC6699", 'TRBV7-1'="#993366",'TRBV7-3'="#CC0033",
        'TRBV7-4'="#FF0033",'TRBV26'="#FF3366", 'TRBV27'="#FF0000",
        'TRBV12-5'="#CC0000",'TRBV15'="#990000",'TRBV2'="#660000",
        'TRBV28'="#993333", 'TRBV4-1'="#993300",'TRBV14'="#663333",
        'TRBV18'="#CC3333",'TRBV5-1'="#FF3333",'TRBV6-8'="#990033",
        'TRBV9'="#FF6666",'TRBV13'="#CC6666",'TRBV4-3'="#FF9999",
        'TRBV24-1'="#FFCCFF",'TRBV10-3'="#660033",'TRBV25-1'="#FF6699",
        'TRBV6-3'="#FFCCCC", 'TRBV7-2'="#CC6666",'TRBJ1-1'="#FF3300",
        'TRBJ1-2'="#FF6633",'TRBJ1-3'="#CC3300",'TRBJ1-4'="#660000",
        'TRBJ1-5'="#990000",'TRBJ1-6'="#CC0000",'TRBJ2-1'="#FF0000",
        'TRBJ2-2'="#993300",'TRBJ2-3'="#CC0033",'TRBJ2-4'="#FF0033",
        'TRBJ2-5'="#990033",'TRBJ2-6'="#FF0000",'TRBJ2-7'="#FF6633")

VJOVA <- as.matrix(VJOVA)
dOVA <- data.frame(from = rep(rownames(VJOVA), times = ncol(VJOVA)),
                     to = rep(colnames(VJOVA), each = nrow(VJOVA)),
                     value = as.vector(VJOVA),
                     stringsAsFactors = FALSE)
chordDiagram(dOVA,grid.col = col, annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(dOVA))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important



#different font characters 

circos.track(track.index = 1, panel.fun = function(x, y) {
circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
            facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.5)

}, bg.border = NA)


#---------------------------------------------------------------------

#Gini-Simpson index (50% rep)

gini.simp_FVIII50 <- repDiversity(DT50FVIII_f, "gini.simp")
mean(gini.simp_FVIII50)

gini.simp_OVA50 <- repDiversity(DT50OVA_f, "gini.simp")
mean(gini.simp_OVA50)

#repertoire total diversity (4donors FVIII_OVA)

gini.simp_FVIII <- repDiversity(DT2[DT2$Antigen == 'FVIII',], "gini.simp")
mean(gini.simp_FVIII)

gini.simp_OVA <- repDiversity(DT2[DT2$Antigen == 'OVA',], "gini.simp")
mean(gini.simp_OVA)

#donor by donor rep totale

DT1179 <-DT2[DT2$Antigen == 'FVIII' & DT2$Donor == '1179',]
gini.simp_FVIII1179 <- repDiversity(DT1179, "gini.simp")
mean(gini.simp_FVIII1179)

DT1179 <-DT2[DT2$Antigen == 'OVA' & DT2$Donor == '1179',]
gini.simp_OVA1179 <- repDiversity(DT1179, "gini.simp")
mean(gini.simp_OVA1179)


DT1180 <-DT2[DT2$Antigen == 'FVIII' & DT2$Donor == '1180',]
gini.simp_FVIII1180 <- repDiversity(DT1180, "gini.simp")
mean(gini.simp_FVIII1180)

DT1180 <-DT2[DT2$Antigen == 'OVA' & DT2$Donor == '1180',]
gini.simp_OVA1180 <- repDiversity(DT1180, "gini.simp")
mean(gini.simp_OVA1180)

DT1195 <-DT2[DT2$Antigen == 'FVIII' & DT2$Donor == '1195',]
gini.simp_FVIII1195 <- repDiversity(DT1195, "gini.simp")
mean(gini.simp_FVIII1195)

DT1195 <-DT2[DT2$Antigen == 'OVA' & DT2$Donor == '1195',]
gini.simp_OVA1195 <- repDiversity(DT1195, "gini.simp")
mean(gini.simp_OVA1195)

DT1207 <-DT2[DT2$Antigen == 'FVIII' & DT2$Donor == '1207',]
gini.simp_FVIII1207 <- repDiversity(DT1207, "gini.simp")
mean(gini.simp_FVIII1207)

DT1207 <-DT2[DT2$Antigen == 'OVA' & DT2$Donor == '1207',]
gini.simp_OVA1207 <- repDiversity(DT1207, "gini.simp")
mean(gini.simp_OVA1207)

#plotting (useless)
diversity <- barplot(height = c(mean(gini.simp_OVA1207), mean(gini.simp_FVIII1207), mean(gini.simp_OVA1195), mean(gini.simp_FVIII1195), mean(gini.simp_OVA1180), mean(gini.simp_FVIII1180), mean(gini.simp_OVA1179), mean(gini.simp_FVIII1179)))
                     
Diversity <- c(mean(gini.simp_FVIII50),mean(gini.simp_OVA50))


barplot(table(Diversity),
        main="Diversity Gini.Simp",
        xlab="Antigen",
        ylab="Diversity",
        border="red",
        col="blue",
        density=10)                     
                   
#-----------------------------------------------------------------------------------------------
#gene usage (50%)

FVIII_gene <- geneUsage(DT50FVIII_f, "hs.trbv", .norm = T, .ambig = "maj")
vis(ova_gene) 

FVIII_gene_J <- geneUsage(DT50FVIII_f, "hs.trbj", .norm = T, .ambig = "maj")
vis(ova_gene_J)


ova_gene <- geneUsage(DT50OVA_f, "hs.trbv", .norm = T, .ambig = "maj")
vis(ova_gene)

ova_gene_J <- geneUsage(DT50OVA_f, "hs.trbj", .norm = T, .ambig = "maj")
vis(ova_gene_J) 

#-------------------------------------------------------------------------------------
#gene usage total (4 donors)

FVIII_genetot <- rbind(DT$data$RUN_1_1179_FVIII, DT$data$RUN_1_1195_FVIII, DT$data$RUN_1_1207_FVIII, DT$data$RUN_3_1180_FVIII) 


FVIII_genetot_V <- geneUsage(FVIII_genetot, "hs.trbv", .norm = T, .ambig = "maj")
vis(FVIII_genetot_V) 

FVIII_genetot_J <- geneUsage(FVIII_genetot, "hs.trbj", .norm = T, .ambig = "maj")
vis(FVIII_genetot_J)
#--------------------------------
OVA_genetot <- rbind(DT$data$RUN_1_1179_OVA, DT$data$RUN_1_1195_OVA, DT$data$RUN_1_1207_OVA, DT$data$RUN_3_1180_OVA)

ova_genetot_V <- geneUsage(OVA_genetot, "hs.trbv", .norm = T, .ambig = "maj")
vis(ova_genetot_V)

ova_genetot_J <- geneUsage(OVA_genetot, "hs.trbj", .norm = T, .ambig = "maj") 
vis(ova_genetot_J) 

#---------------------------------------------------------------------------------------
public_rep_tot <- pubRep(
 DT$data,
  .col = "aa",
  .quant = c("prop"),
  .coding = TRUE,
  .min.samples = 1,
  .max.samples = NA,
  .verbose = TRUE
)

write.csv(public_rep_tot, file = "public_repOVA_FVIII.csv")
write.table(public_rep_tot, file = "public_repOVA_FVIII.txt", sep = "\t",
            row.names = TRUE, col.names = NA)


#---------------------------------------------------------------------------------------
FVIII_common <- intersect(DT$data$RUN_1_1179_FVIII$CDR3.aa, DT$data$RUN_1_1195_FVIII$CDR3.aa, DT$data$RUN_1_1207_FVIII$CDR3.aa, DT$data$RUN_3_1180_FVIII$CDR3.aa)


public_rep_FVIII <- pubRep(
  FVIII_common,
  .col = "aa",
  .quant = c("prop"),
  .coding = TRUE,
  .min.samples = 1,
  .max.samples = NA,
  .verbose = TRUE
)

#-----------------------------------------------------------------------------
#finding public sequences (4donors FVIII vs 5donors FVIII&POOL)
DTpool <- repLoad(.path = "/mnt/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/_Analyse_fichiers_corriges_FVIII&POOL_unique")


DTpool$data$RUN_2_1238_FVIII$Donor <-"1238"
DTpool$data$RUN_2_1238_FVIII$Antigen<-"FVIII"

DTpool$data$RUN_2_1238_FVIII_POOL$Donor <-"1238"
DTpool$data$RUN_2_1238_FVIII_POOL$Antigen<-"POOL"
#--------------------------------------------
DTpool$data$RUN_2_1245_FVIII$Donor <-"1245"
DTpool$data$RUN_2_1245_FVIII$Antigen<-"FVIII"

DTpool$data$RUN_2_1245_FVIII_POOL$Donor <-"1245"
DTpool$data$RUN_2_1245_FVIII_POOL$Antigen<-"POOL"
#-------------------------------------------
DTpool$data$RUN_2_1255_FVIII$Donor <-"1255"
DTpool$data$RUN_2_1255_FVIII$Antigen<-"FVIII"

DTpool$data$RUN_2_1255_FVIII_POOL$Donor <-"1255"
DTpool$data$RUN_2_1255_FVIII_POOL$Antigen<-"POOL"
#---------------------------------------------
DTpool$data$RUN_2_1285_FVIII$Donor <-"1285"
DTpool$data$RUN_2_1285_FVIII$Antigen<-"FVIII"

DTpool$data$RUN_2_1285_FVIII_POOL$Donor <-"1285"
DTpool$data$RUN_2_1285_FVIII_POOL$Antigen<-"POOL"
#----------------------------------------------
DTpool$data$RUN_3_1327_FVIII$Donor <-"1327"
DTpool$data$RUN_3_1327_FVIII$Antigen<-"FVIII"

DTpool$data$RUN_3_1327_FVIII_POOL$Donor <-"1327"
DTpool$data$RUN_3_1327_FVIII_POOL$Antigen<-"POOL"


DTpool5 <- rbind(DTpool$data$RUN_2_1238_FVIII_POOL,DTpool$data$RUN_2_1245_FVIII_POOL, DTpool$data$RUN_2_1255_FVIII_POOL, DTpool$data$RUN_2_1285_FVIII_POOL, DTpool$data$RUN_3_1327_FVIII_POOL)

DTFVIII5 <- rbind(DTpool$data$RUN_2_1238_FVIII, DTpool$data$RUN_2_1245_FVIII, DTpool$data$RUN_2_1255_FVIII, DTpool$data$RUN_2_1285_FVIII, DTpool$data$RUN_3_1327_FVIII)





intersect(DT50FVIII_f$CDR3.aa, DTFVIII5$CDR3.aa) #34 public clones
intersect(DT50FVIII_f$CDR3.aa, DTOVA_f$CDR3.aa)  #7 clones 
