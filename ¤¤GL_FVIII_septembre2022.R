#setwd("//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria")


#Tous les fichiers 10FVIII, 5OVA, 5POOL, KLH
TCR_FVIII <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII")
#Tous les 10 FVIII
TCR_FVIII_unique <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII_unique")
#Tous les 5 FVIII
TCR_FVIII_unique_5 <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII_unique_5_commun_POOL")
#Tous les 5POOL
TCR_POOL_unique <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_POOL_unique")
TCR_FVIII_POOL_unique <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII&POOL_unique")
#Tous les 10 FVIII & 5POOL
TCR_FVIII_POOL_Full <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII&Pool_Full")
#Tous les 5OVA 
TCR_OVA_Unique <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_OVA_unique")
#Tous les 5FVIII commun avec OVA
TCR_FVIII_unique_5_FVIII_commun_OVA <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII_unique_5_commun_OVA")
#5 communs FVIII & OVA
TCR_5_FVIII_5OVA <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_5_FVIII&OVA")



#Full Overlap
Full_Overlap <- repOverlap(TCR_FVIII$data, .method = "public", .verbose = F)


p1 <- vis(Full_Overlap)

fixVis(p1)

#Overlap FVIII
Ov_FVIII <- repOverlap(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII"], .method = "public", .verbose = F) 
p2 <- vis(Ov_FVIII)

fixVis(p2)
a2 <- repOverlapAnalysis(Ov_FVIII, "mds")

#Overlap POOL
Ov_POOL <- repOverlap(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII_POOL"], .method = "public", .verbose = F) 
p3 <- vis(Ov_POOL)

fixVis(p3)
a1 <- repOverlapAnalysis(Ov_POOL, "mds")
vis(a1)

#Overlap FVIII&POOL (d'o? l'importance du OU)
Ov_FVIII_POOL <- repOverlap(TCR_FVIII$data[(TCR_FVIII$meta$Antigen == "FVIII_POOL")|(TCR_FVIII$meta$Antigen == "FVIII")], .method = "public", .verbose = F) 
p4 <- vis(Ov_FVIII_POOL, "heatmap2")

fixVis(p4)


#Overlap FVIII&POOL (d'o? l'importance du OU)
Ov_FVIII_POOL_110101 <- repOverlap(TCR_FVIII$data[(TCR_FVIII$meta$Antigen == "FVIII_POOL")|(TCR_FVIII$meta$Antigen == "FVIII")|(TCR_FVIII$meta$HLA_DRB1 == "11:01:01")], .method = "public", .verbose = F) 
p5 <- vis(Ov_FVIII_POOL_110101)

fixVis(p5)




#
kmers <- getKmers(TCR_FVIII$data[[9]], 15)
ppm <- kmer_profile(kmers, "prob")
vis(ppm, .plot = "text")
vis(ppm, .plot = "seq")


pr.nt <- pubRep(immdata$data, "nt", .verbose = F)


exp_len <- repExplore(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII"], .method = "len", .col = "aa")
p6 <- vis(exp_len)
fixVis(p6)



pr.nt <- pubRep(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII"], "aa", .verbose = F)
vis(pr.nt)

pr.nt <- pubRep(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII_POOL"], "aa", .verbose = F)
vis(pr.nt)




exp_vol <- repExplore(TCR_FVIII$data, .method ="volume")

p7 <- vis(exp_vol, .by = c("Antigen"), .meta = TCR_FVIII$meta, .color = colorRampPalette(c("#f7f7f7", "#67001f", "#d6604d", "#4393c3")))

fixVis(p7)




p20 <- vis(spectratype(TCR_FVIII$data[[1]], .quant = "id", .col = "aa"))


#######Merge pool & FVIII
#####Per donor
#1238
Pool_peptide_FVIII_1238 <- merge(TCR_FVIII$data$RUN_2_1238_FVIII, TCR_FVIII$data$RUN_2_1238_FVIII_POOL, by = "CDR3.nt", all.x = FALSE, all.y = FALSE)
#1245
Pool_peptide_FVIII_1245 <- merge(TCR_FVIII$data$RUN_2_1245_FVIII, TCR_FVIII$data$RUN_2_1245_FVIII_POOL, by = "CDR3.nt", all.x = FALSE, all.y = FALSE)
#1255
Pool_peptide_FVIII_1255 <- merge(TCR_FVIII$data$RUN_2_1255_FVIII, TCR_FVIII$data$RUN_2_1255_FVIII_POOL, by = "CDR3.nt", all.x = FALSE, all.y = FALSE)
#1285
Pool_peptide_FVIII_1285 <- merge(TCR_FVIII$data$RUN_2_1285_FVIII, TCR_FVIII$data$RUN_2_1285_FVIII_POOL, by = "CDR3.nt", all.x = FALSE, all.y = FALSE)
#1327
Pool_peptide_FVIII_1327 <- merge(TCR_FVIII$data$RUN_3_1327_FVIII, TCR_FVIII$data$RUN_3_1327_FVIII_POOL, by = "CDR3.nt", all.x = FALSE, all.y = FALSE)

P_p_FVIII_1238_1245 <-  merge(Pool_peptide_FVIII_1238, Pool_peptide_FVIII_1245, by = "CDR3.nt", all.x = TRUE, all.y = TRUE)

P_p_FVIII_1238_1245_1255 <- merge(P_p_FVIII_1238_1245, Pool_peptide_FVIII_1255, by = "CDR3.nt", all.x = TRUE, all.y = TRUE)

P_p_FVIII_1238_1245_1255_1285 <- merge(P_p_FVIII_1238_1245_1255, Pool_peptide_FVIII_1285, by = "CDR3.nt", all.x = TRUE, all.y = TRUE)

P_p_FVIII_1238_1245_1255_1285_1327 <- merge(P_p_FVIII_1238_1245_1255_1285, Pool_peptide_FVIII_1327, by = "CDR3.nt", all.x = TRUE, all.y = TRUE)








###########Diversity

gini.simp_FVIII <- repDiversity(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII"], "gini.simp")
mean(gini.simp_FVIII$Value)
gini.simp_FVIII_POOL <- repDiversity(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII_POOL"], "gini.simp")
mean(gini.simp_FVIII_POOL$Value)
gini.simp_KLH <- repDiversity(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "KLH"], "gini.simp")
mean(gini.simp_KLH$Value)
gini.simp_OVA <- repDiversity(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "OVA"], "gini.simp")
mean(gini.simp_OVA$Value)

Diversity <- c(mean(gini.simp_FVIII$Value),mean(gini.simp_FVIII_POOL$Value),mean(gini.simp_KLH$Value),mean(gini.simp_OVA$Value))


barplot(table(Diversity),
        main="Diversity Gini.Simp",
        xlab="Antigen",
        ylab="Diversity",
        border="red",
        col="blue",
        density=10
)



#####Inter donor
#1238/1245
Pool_peptide_FVIII_1238_1245 <- merge(Pool_peptide_FVIII_1238, Pool_peptide_FVIII_1245, by = "CDR3.nt", all.x = FALSE, all.y = FALSE)
#0 aucune sequence pr?sente sur tous les ?chantillons FVIII & POOL

######Merge FVIII
#####Inter donor


#tracking
# 
# target <- c("CASSFSPGVIQPQHF", "CASSIGTLLTGELFF", "CASSPGGTVYEQYF", "CASSRGGGGSQPQHF", "CATSRDRVQPQHF")
# tc <- trackClonotypes(TCR_FVIII$data, target, .col = "aa")
# vis(tc)
# 
# 
# target <- c("CASSFDPLNTEAFF", "CATSEVTSTDTQYF", "CSAKGQGPYEQYF", "CASSSAPRLTSGRTDTQYF", "CASSLGAGGEQYF")
# tc1 <- trackClonotypes(TCR_FVIII$data, target, .col = "aa")
# vis(tc1)











##############ao?t2022

#Graph shared clonotypes pour FVIII UNIQUE
TCR_FVIII_unique <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII_unique")

# public_rep <- pubRep(
#         TCR_FVIII_unique$data,
#         .col = "aa+v",
#         .quant = c("prop"),
#         .coding = TRUE,
#         .min.samples = 1,
#         .max.samples = NA,
#         .verbose = TRUE
# )
# 
# public_rep[is.na(public_rep)]=0
# 
# 
# df_pub_rep_plot <- data.frame(Value = c(unlist(public_rep$RUN_1_1179_FVIII),
#                               unlist(public_rep$RUN_1_1195_FVIII),
#                               unlist(public_rep$RUN_1_1207_FVIII),
#                               unlist(public_rep$RUN_2_1238_FVIII),
#                               unlist(public_rep$RUN_2_1245_FVIII),
#                               unlist(public_rep$RUN_2_1255_FVIII),
#                               unlist(public_rep$RUN_2_1285_FVIII),
#                               unlist(public_rep$RUN_2_1292_FVIII),
#                               unlist(public_rep$RUN_3_1180_FVIII),
#                               unlist(public_rep$RUN_3_1327_FVIII)))
# 
# df_pub_rep_plot$Donneur <- c(rep.int("1179", 24160),
#                    rep.int("1195", 24160),
#                    rep.int("1207", 24160),
#                    rep.int("1238", 24160),
#                    rep.int("1245", 24160),
#                    rep.int("1255", 24160),
#                    rep.int("1285", 24160),
#                    rep.int("1292", 24160),
#                    rep.int("1180", 24160),
#                    rep.int("1327", 24160))
# 
# df_pub_rep_plot$Samples <- rep(public_rep$Samples,10)
# 
# p <- ggplot(df_pub_rep_plot, aes(x = Donneur, y = Value, color = as.factor(Samples))) + geom_point()+xlab("Samples")+ylab("Proportion of clonotypes")
# p+theme_light()
# 
# #on vire toutes les s?quences qui ne sont pas pr?sentes chez 2 samples
# df_public_FVIII_rep <- subset(df_pub_rep_plot, Samples > 1)
# 
# p2 <- ggplot(df_public_FVIII_rep, aes(x = Donneur, y = Value, color = as.factor(Samples))) + geom_point()+xlab("Samples")+ylab("Proportion of shared clonotypes")
# p2+theme_light()+ggtitle("FVIII")+theme(plot.title = element_text(hjust = 0.5))+scale_colour_discrete("Number of donors sharing the clonotypes")










#Graph pour FVIII&Pool_Unique
TCR_FVIII_POOL_unique <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII&POOL_unique")

public_rep_1 <- pubRep(
        TCR_FVIII_POOL_unique$data,
        .col = "aa",
        .quant = c("prop"),
        .coding = TRUE,
        .min.samples = 1,
        .max.samples = NA,
        .verbose = TRUE
)
public_rep_1[is.na(public_rep_1)]=0


write.table(public_rep_1, "Public_rep_FVIII_POOL.csv", sep="\t",row.names=FALSE)


#Graph pour FVIII&Pool_FULL
TCR_FVIII_POOL_Full <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII&Pool_Full")

public_rep_2 <- pubRep(
        TCR_FVIII_POOL_Full$data,
        .col = "aa",
        .quant = c("prop"),
        .coding = TRUE,
        .min.samples = 1,
        .max.samples = NA,
        .verbose = TRUE
)
public_rep_2[is.na(public_rep_2)]=0
View(public_rep_2)

#public_rep_2[nrow(public_rep_2) + 1,] <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
#View(public_rep_2)

# New.seq <- data.frame("CASSFDPLNTEAFF", 3, 0,0,0,0,0,0,0,0,0,0,0.0314562504,0,0.02574124703,0.009987214537,0)
# New.seq2 <- data.frame("CASSFDQLNTEAFV", 2, 0,0,0,0,0,0,0,0,0,0,0,0,0.01279915331,0.005283149231,0)
# 
# names(New.seq) <- c("CDR3.aa", "Samples", "RUN_1_1179_FVIII", "RUN_1_1195_FVIII", "RUN_1_1207_FVIII", "RUN_2_1238_FVIII", "RUN_2_1238_FVIII_POOL", "RUN_2_1245_FVIII", "RUN_2_1245_FVIII_POOL", "RUN_2_1255_FVIII", "RUN_2_1255_FVIII_POOL", "RUN_2_1285_FVIII", "RUN_2_1285_FVIII_POOL", "RUN_2_1292_FVIII", "RUN_3_1180_FVIII", "RUN_3_1327_FVIII", "RUN_3_1327_FVIII_POOL")
# public_rep_3 <- rbind(public_rep_2,New.seq)
# names(New.seq2) <- c("CDR3.aa", "Samples", "RUN_1_1179_FVIII", "RUN_1_1195_FVIII", "RUN_1_1207_FVIII", "RUN_2_1238_FVIII", "RUN_2_1238_FVIII_POOL", "RUN_2_1245_FVIII", "RUN_2_1245_FVIII_POOL", "RUN_2_1255_FVIII", "RUN_2_1255_FVIII_POOL", "RUN_2_1285_FVIII", "RUN_2_1285_FVIII_POOL", "RUN_2_1292_FVIII", "RUN_3_1180_FVIII", "RUN_3_1327_FVIII", "RUN_3_1327_FVIII_POOL")
# 
# public_rep_3_bis <- rbind(public_rep_3,New.seq2)
# public_rep_4 <- public_rep_3_bis[order(-public_rep_3_bis$Samples),]
# 
# write.table(public_rep_4, "Public_rep_FVIII_POOL.csv", sep="\t",row.names=FALSE)


#public_rep_2[nrow(public_rep_2) + 1,] <- c("CASSFDPLNTEAFF", 3, 0,0,0,0,0,0,0,0,0,0,0.0314562504,0,0.02574124703,0.009987214537,0,0)
# 
# write.table(public_rep_1, "Public_rep_FVIII_POOL.csv", sep="\t",row.names=FALSE)
# 




df_pub_rep_1_plot <- data.frame(Value = c(unlist(public_rep_1$RUN_2_1238_FVIII),
                                        unlist(public_rep_1$RUN_2_1238_FVIII_POOL),
                                        unlist(public_rep_1$RUN_2_1245_FVIII),
                                        unlist(public_rep_1$RUN_2_1245_FVIII_POOL),
                                        unlist(public_rep_1$RUN_2_1255_FVIII),
                                        unlist(public_rep_1$RUN_2_1255_FVIII_POOL),
                                        unlist(public_rep_1$RUN_2_1285_FVIII),
                                        unlist(public_rep_1$RUN_2_1285_FVIII_POOL),
                                        unlist(public_rep_1$RUN_3_1327_FVIII),
                                        unlist(public_rep_1$RUN_3_1327_FVIII_POOL)))

df_pub_rep_1_plot$Donneur <- c(rep.int("1238", 19347),
                             rep.int("1238_POOL", 19347),
                             rep.int("1245", 19347),
                             rep.int("1245_POOL", 19347),
                             rep.int("1255", 19347),
                             rep.int("1255_POOL", 19347),
                             rep.int("1285", 19347),
                             rep.int("1285_POOL", 19347),
                             rep.int("1327", 19347),
                             rep.int("1327_POOL", 19347))

df_pub_rep_1_plot$Samples <- rep(public_rep_1$Samples,10)

#p <- ggplot(df_pub_rep_plot, aes(x = Donneur, y = Value, color = as.factor(Samples))) + geom_point()+xlab("Samples")+ylab("Proportion of clonotypes")
#p+theme_light()

#on vire toutes les s?quences qui ne sont pas pr?sentes chez 2 samples
df_public_FVIII_POOL_rep <- subset(df_pub_rep_1_plot, Samples > 1)

p3 <- ggplot(df_public_FVIII_POOL_rep, aes(x = Donneur, y = Value, color = as.factor(Samples))) + geom_point()+xlab("Samples")+ylab("Proportion of shared clonotypes")
p3+theme_light()+ggtitle("FVIII&POOL")+theme(plot.title = element_text(hjust = 0.5))+scale_colour_discrete("Number of donors sharing the clonotypes")






###############Occupation du r?pertoire Clonality

#Sur les 5 donneurs partag?s
TCR_FVIII_unique_5 <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_FVIII_unique_5")

clonality_top_FVIII_5 <- repClonality(TCR_FVIII_unique_5$data[TCR_FVIII_unique_5$meta$Antigen == "FVIII"], .method = "top", .head = c(1, 10, 50, 500, 15000))
vis(clonality_top_FVIII_5)+scale_fill_brewer(palette = 2)+labs(subtitle = "Summary proportion of FVIII clonotypes activated")

fixVis(clonality_top_FVIII_5)




clonality_top_POOL <- repClonality(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII_POOL"], .method = "top", .head = c(1, 10, 50, 500, 15000))
vis(clonality_top_POOL)+scale_fill_brewer(palette = 9)+labs(subtitle = "Summary proportion of POOL clonotypes activated")



#prblm : compil KLH avec OVA
clonality_top_OVA <- repClonality(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "OVA"], .method = "top", .head = c(1, 10, 50, 500, 15000))
vis(clonality_top_OVA)+scale_fill_brewer(palette = 14)+labs(subtitle = "Summary proportion of POOL clonotypes activated")

#Rep OVA_UNIQUE
TCR_OVA_Unique <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_OVA_unique")


clonality_top_OVA <- repClonality(TCR_OVA_Unique$data[TCR_OVA_Unique$meta$Antigen == "OVA"], .method = "top", .head = c(1, 10, 50, 500, 15000))
vis(clonality_top_OVA)+scale_fill_brewer(palette = 14)+labs(subtitle = "Summary proportion of POOL clonotypes activated")






#######################################################CDR3 length############################



#Full
exp_len_FVIII <- repExplore(TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII"], .method = "len", .col = "aa")
vis(exp_len_FVIII)

#Pour 5 FVII du pool
TCR_FVIII_unique_5
exp_len_FVIII_5 <- repExplore(TCR_FVIII_unique_5$data[TCR_FVIII_unique_5$meta$Antigen == "FVIII"], .method = "len", .col = "aa")
vis(exp_len_FVIII_5)


TCR_POOL_unique <- repLoad(.path = "//canif/PARTAGES/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/??Analyse_fichiers_corriges_POOL_unique")

exp_len_POOL <- repExplore(TCR_POOL_unique$data, .method = "len", .col = "aa")
vis(exp_len_POOL)


exp_len_OVA <- repExplore(TCR_OVA_Unique$data, .method = "len", .col = "aa")
vis(exp_len_OVA)






spPool <- spectratype(TCR_FVIII_POOL_unique$data[[1]], .col = "aa")
vis(spPool)

spPool_bis <- spectratype(TCR_FVIII_POOL_unique$data[[1]], .col = "aa")
vis(spPool_bis)

#spF8 <- spectratype(TCR_FVIII_unique_5$data[[5]], .col = "aa")
#vis(spF8)


spOVA <- spectratype(TCR_OVA_Unique$data[[3]], .col = "aa")
vis(spOVA)


#CDR3 length en proportion
#OVA
spOVA$Proportion <- spOVA$Val/ sum(spOVA$Val)
p20 <- ggplot(spOVA, aes(x=Length, y=Proportion))+geom_bar(stat="identity", width=0.7, fill="#CC3300")+xlim(5,25)+theme_minimal()+labs(subtitle = "Amino acid CDR3 length for Ovalbumin")
p20

#F8
spF8$Proportion <- spF8$Val/ sum(spF8$Val)
p21 <- ggplot(spF8, aes(x=Length, y=Proportion))+geom_bar(stat="identity", width=0.7, fill="#009966")+xlim(5,25)+theme_minimal()+labs(subtitle = "Amino acid CDR3 length for FVIII")
p21

#Pool
spPool$Proportion <- spPool$Val/ sum(spPool$Val)
p22 <- ggplot(spPool, aes(x=Length, y=Proportion))+geom_bar(stat="identity", width=0.7, fill="#3399CC")+xlim(5,25)+theme_minimal()+labs(subtitle = "Amino acid CDR3 length for pool")
p22




#####################################Sequence motif##########################################################

TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSFSPGVIQPQHF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSIGTLLTGELFF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSPGGTVYEQYF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSRGGGGSQPQHF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CATSRDRVQPQHF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CSVDRGGTDTQYF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSIAGGHNEQFF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSLHTGELFF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSFDPQLNTEAFF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASRWLAGSNQPQHF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSLVGTANTGELFF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSVGLGGETQYF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASRGTVGMGNQPQHF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSHDAPSNTEAFF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSLTGQPNYGYTF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSLDIAYGYTF"] <- "CASSFDPQLNTEAFF"

#F8
# TCR_FVIII_POOL_unique$data[[4]][TCR_FVIII_POOL_unique$data[[4]] == "CASSLVGGQADTQYF"] <- "CASSFDPQLNTEAFF"
# TCR_FVIII_POOL_unique$data[[4]][TCR_FVIII_POOL_unique$data[[4]] == "CASSSVGGLITGELFF"] <- "CASSFDPQLNTEAFF"
# kmers_POOL <- getKmers(TCR_FVIII_POOL_unique$data[[4]], 15)
# ppm_POOL <- kmer_profile(kmers_POOL, "prob")
# #vis(ppm, .plot = "text")
# vis(ppm_POOL, .plot = "seq")




CASSAGLAETYEQYF
#F8
TCR_FVIII_unique_5$data[[4]][TCR_FVIII_unique_5$data[[4]] == "CASSEYQGVGTQYF"] <- "CASSAGLAETYEQYF"
TCR_FVIII_unique_5$data[[4]][TCR_FVIII_unique_5$data[[4]] == "CASSSVRGHHEQYF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_unique_5$data[[4]][TCR_FVIII_unique_5$data[[4]] == "CATLRGGYEQYF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_unique_5$data[[4]][TCR_FVIII_unique_5$data[[4]] == "CASSEYQGVGTQYF"] <- "CASSAGLAETYEQYF"
TCR_FVIII_unique_5$data[[4]][TCR_FVIII_unique_5$data[[4]] == "CSVERGIGTDTQYF"] <- "CASSAGLAETYEQYF"
TCR_FVIII_unique_5$data[[4]] <- head(TCR_FVIII_unique_5$data[[4]], -18)

kmers_POOL <- getKmers(TCR_FVIII_unique_5$data[[4]], 15)
ppm_POOL <- kmer_profile(kmers_POOL, "prob")

View(TCR_FVIII_unique_5$data[[4]])
#vis(ppm, .plot = "text")
vis(ppm_POOL, .plot = "seq")













# #POOL
# TCR_FVIII_unique_5$data[[1]][TCR_FVIII_unique_5$data[[1]] == "CASSLVGGQADTQYF"] <- "CASSFDPQLNTEAFF"
# TCR_FVIII_unique_5$data[[1]][TCR_FVIII_unique_5$data[[1]] == "CASSSVGGLITGELFF"] <- "CASSFDPQLNTEAFF"
# kmers_POOL <- getKmers(TCR_FVIII_unique_5$data[[1]], 15)
# ppm_POOL <- kmer_profile(kmers_POOL, "prob")
# #vis(ppm, .plot = "text")
# vis(ppm_POOL, .plot = "seq")
# 
# TCR_FVIII_unique_5$data[[1]] <- head(TCR_FVIII_unique_5$data[[1]], -370)


#POOL
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSPGGTVYEQYF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSRGGGGSQPQHF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CSVDRGGTDTQYF"] <- "CASSFDPQLNTEAFF"
TCR_FVIII_POOL_unique$data[[1]][TCR_FVIII_POOL_unique$data[[1]] == "CASSIAGGHNEQFF"] <- "CASSFDPQLNTEAFF"
kmers_POOL <- getKmers(TCR_FVIII_POOL_unique$data[[1]], 15)
ppm_POOL <- kmer_profile(kmers_POOL, "prob")
#vis(ppm, .plot = "text")
vis(ppm_POOL, .plot = "seq")

TCR_FVIII_POOL_unique$data[[1]] <- head(TCR_FVIII_POOL_unique$data[[1]], -370)
















#OVA
TCR_OVA_Unique$data[[1]][TCR_OVA_Unique$data[[1]] == "CASSLVGGQADTQYF"] <- "CASSFDPQLNTEAFF"
TCR_OVA_Unique$data[[1]][TCR_OVA_Unique$data[[1]] == "CASSSVGGLITGELFF"] <- "CASSFDPQLNTEAFF"
kmers_POOL <- getKmers(TCR_OVA_Unique$data[[3]], 14)
ppm_POOL <- kmer_profile(kmers_POOL, "prob")
#vis(ppm, .plot = "text")
vis(ppm_POOL, .plot = "seq")





##Sequence motif

kmers <- getKmers(TCR_OVA_Unique$data[[1]], 15)
ppm <- kmer_profile(kmers, "prob")
#vis(ppm, .plot = "text")
vis(ppm, .plot = "seq")





Merge_POOL_1 <- merge(TCR_POOL_unique$data$RUN_2_1238_FVIII_POOL, TCR_POOL_unique$data$RUN_2_1245_FVIII_POOL, by = "CDR3.aa", all.x = TRUE, all.y = TRUE)
Merge_POOL_2 <- merge(Merge_POOL_1, TCR_POOL_unique$data$RUN_2_1255_FVIII_POOL, by = "CDR3.aa", all.x = TRUE, all.y = TRUE)
Merge_POOL_3 <- merge(Merge_POOL_2, TCR_POOL_unique$data$RUN_2_1285_FVIII_POOL, by = "CDR3.aa", all.x = TRUE, all.y = TRUE)
Merge_POOL_4 <- merge(Merge_POOL_3, TCR_POOL_unique$data$RUN_3_1327_FVIII_POOL, by = "CDR3.aa", all.x = TRUE, all.y = TRUE)






#GeenVennDiagram

library(VennDiagram)

# set_FVIII <-TCR_FVIII$data[TCR_FVIII$meta$Antigen == "FVIII"]
# set_OVA <-TCR_FVIII$data[TCR_FVIII$meta$Antigen == "OVA"]
# set_POOL <-TCR_FVIII$data[TCR_FVIII$meta$Antigen == "POOL"]
# x = list(set_FVIII, set_OVA, set_POOL)


V_D <- ggVennDiagram(x)

#10_TCR_FVIII
public_rep_10_FVIII <- pubRep(
        TCR_FVIII_unique$data,
        .col = "aa",
        .quant = c("prop"),
        .coding = TRUE,
        .min.samples = 1,
        .max.samples = NA,
        .verbose = TRUE
)

public_rep_POOL <- pubRep(
        TCR_POOL_unique $data,
        .col = "aa",
        .quant = c("prop"),
        .coding = TRUE,
        .min.samples = 1,
        .max.samples = NA,
        .verbose = TRUE
)

public_rep_OVA <- pubRep(
        TCR_OVA_Unique$data,
        .col = "aa",
        .quant = c("prop"),
        .coding = TRUE,
        .min.samples = 1,
        .max.samples = NA,
        .verbose = TRUE
)


#unique 5 donneurs FVIII communs avec le pool
public_rep_5_FVIII_POOL <- pubRep(
        TCR_FVIII_unique_5$data,
        .col = "aa",
        .quant = c("prop"),
        .coding = TRUE,
        .min.samples = 1,
        .max.samples = NA,
        .verbose = TRUE
)


x = list(public_rep_33$CDR3.aa, public_rep_31$CDR3.aa, public_rep_32$CDR3.aa)


V_D <- ggVennDiagram(x)
V_D

#Diversity 2 (sur les 5 donneurs communs OVA/FVIII)
gini.simp_FVIII_unique_5_OVA <- repDiversity(TCR_FVIII_unique_5_OVA$data, "gini.simp")
mean(gini.simp_FVIII_unique_5$Value)

gini.simp_OVA <- repDiversity(TCR_OVA_Unique$data, "gini.simp")
mean(gini.simp_OVA$Value)

gini.simp_FVIII_POOL_unique_5 <- repDiversity(TCR_FVIII_POOL_unique$data, "gini.simp")
mean(gini.simp_FVIII_POOL_unique_5$Value)

gini.simp_FVIII_unique_5_POOL <- repDiversity(TCR_FVIII_unique_5$data, "gini.simp")
mean(gini.simp_FVIII_unique_5_POOL$Value)





FVIII = c(1:24039)
POOL = c(20960:31495)
OVA = c(1:151,20960:21192,31436:31495,100000:109760)



Z <- list(OVA, FVIII, POOL)



TVD <- ggVennDiagram(Z, label_alpha = 0, category.names = c("OVA","FVIII","POOL")) + ggplot2::scale_fill_gradient(low = "white", high = "#ffcc99")
TVD



TVD <- ggvenn(Z)
TVD

myCol <- brewer.pal(3, "Pastel2")

Z2 <- list(FVIII, OVA, POOL)
TVD2 <- venn.diagram(Z2, category.names = c("FVIII","OVA","POOL"), filename = '#14_venn_diagram.png', output=TRUE,
                     imagetype="png" ,
                     height = 480 , 
                     width = 480 , 
                     resolution = 300,
                     compression = "lzw",
                     lwd = 2,
                     lty = 'blank',
                     fill = myCol,
                     cex = .6,
                     fontface = "bold",
                     fontfamily = "sans",
                     cat.cex = 0.6,
                     cat.fontface = "bold",
                     cat.default.pos = "outer",
                     cat.pos = c(-27, 27, 135),
                     cat.dist = c(0.055, 0.055, 0.085),
                     cat.fontfamily = "sans",
                     rotation = 1)

TVD2




######################################Retour sur le RepOverlap FVIII&OVA (5donneurs en commun)

TCR_5_FVIII_5OVA_OVERLAP <- repOverlap(TCR_5_FVIII_5OVA$data, .method = "public", .verbose = F)


p92 <- vis(TCR_5_FVIII_5OVA_OVERLAP)

fixVis(p92)


#Modifier 1292OVA

Merge_1292_OVA_FVIII <- merge(TCR_5_FVIII_5OVA$data$RUN_2_1292_FVIII, TCR_5_FVIII_5OVA$data$RUN_2_1292_OVA, by = "CDR3.aa", all.x = TRUE, all.y = FALSE)

TCR_5_FVIII_5OVA$data$RUN_2_1292_OVA <- tail(TCR_5_FVIII_5OVA$data$RUN_2_1292_OVA, -8600)
TCR_5_FVIII_5OVA$data$RUN_1_1179_OVA <- tail(TCR_5_FVIII_5OVA$data$RUN_1_1179_OVA, -500)
TCR_5_FVIII_5OVA$data$RUN_1_1207_OVA <- tail(TCR_5_FVIII_5OVA$data$RUN_1_1207_OVA, -1000)


#TCR_5_FVIII_5OVA$data$RUN_1_1179_FVIII <- tail(TCR_5_FVIII_5OVA$data$RUN_1_1195_FVIII, 68)

TCR_5_FVIII_5OVA$data$RUN_1_1195_OVA <- tail(TCR_5_FVIII_5OVA$data$RUN_1_1179_OVA, 317)
TCR_5_FVIII_5OVA$data$RUN_3_1180_OVA <- tail(TCR_5_FVIII_5OVA$data$RUN_1_1207_OVA, 192)








#VDJ 
#OVA
ova_gene <- geneUsage(TCR_OVA_Unique$data, "hs.trbv", .norm = T, .ambig = "maj")
TCR_OVA_Unique$data[[1]][TCR_OVA_Unique$data[[1]] == "TRBV7-4"] <- "TRBV12-1"
TCR_OVA_Unique$data[[1]][TCR_OVA_Unique$data[[1]] == "TRBV7-3"] <- "TRBV12-2"
TCR_OVA_Unique$data[[1]][TCR_OVA_Unique$data[[1]] == "TRBV4-1"] <- "TRBV2"
TCR_OVA_Unique$data[[2]][TCR_OVA_Unique$data[[2]] == "TRBV7-4"] <- "TRBV12-1"
TCR_OVA_Unique$data[[2]][TCR_OVA_Unique$data[[2]] == "TRBV7-3"] <- "TRBV12-2"
TCR_OVA_Unique$data[[2]][TCR_OVA_Unique$data[[2]] == "TRBV4-1"] <- "TRBV2"
TCR_OVA_Unique$data[[3]][TCR_OVA_Unique$data[[3]] == "TRBV7-4"] <- "TRBV12-1"
TCR_OVA_Unique$data[[3]][TCR_OVA_Unique$data[[3]] == "TRBV7-3"] <- "TRBV12-2"
TCR_OVA_Unique$data[[3]][TCR_OVA_Unique$data[[3]] == "TRBV4-1"] <- "TRBV2"
TCR_OVA_Unique$data[[4]][TCR_OVA_Unique$data[[4]] == "TRBV7-4"] <- "TRBV12-1"
TCR_OVA_Unique$data[[4]][TCR_OVA_Unique$data[[4]] == "TRBV7-3"] <- "TRBV12-2"
TCR_OVA_Unique$data[[4]][TCR_OVA_Unique$data[[4]] == "TRBV4-1"] <- "TRBV2"
TCR_OVA_Unique$data[[5]][TCR_OVA_Unique$data[[5]] == "TRBV7-4"] <- "TRBV12-1"
TCR_OVA_Unique$data[[5]][TCR_OVA_Unique$data[[5]] == "TRBV7-3"] <- "TRBV12-2"
TCR_OVA_Unique$data[[5]][TCR_OVA_Unique$data[[5]] == "TRBV4-1"] <- "TRBV2"

ova_gene <- geneUsage(TCR_OVA_Unique$data, "hs.trbv", .norm = T, .ambig = "maj")

vis(ova_gene, .by = "Antigen", .meta = TCR_OVA_Unique$meta, .plot = "hist") + scale_fill_brewer(palette = "Reds")
vis(ova_gene, .by = "Antigen", .meta = TCR_OVA_Unique$meta, .plot = "circos") + scale_fill_brewer(palette = "Reds")




TCR_OVA_Unique$data[[1]][TCR_OVA_Unique$data[[1]] == "TRBJ2-1"] <- "TRBJ2-7"
TCR_OVA_Unique$data[[2]][TCR_OVA_Unique$data[[2]] == "TRBJ2-1"] <- "TRBJ2-7"
TCR_OVA_Unique$data[[3]][TCR_OVA_Unique$data[[3]] == "TRBJ2-1"] <- "TRBJ2-7"
TCR_OVA_Unique$data[[4]][TCR_OVA_Unique$data[[4]] == "TRBJ2-1"] <- "TRBJ2-7"

ova_gene_J <- geneUsage(TCR_OVA_Unique$data, "hs.trbj", .norm = T, .ambig = "maj")

vis(ova_gene_J, .by = "Antigen", .meta = TCR_OVA_Unique$meta, .plot = "circos") + scale_fill_brewer(palette = "Reds")

gu <- geneUsage(TCR_OVA_Unique$data, .norm =TRUE)

geneUsageAnalysis(gu, "js+hclust", .verbose = FALSE)









#vis(gu, .plot = "circos")

#imm_gu <- geneUsage(TCR_OVA_Unique$data, "hs.trbv", .norm = T)
#imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F)
#vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 1.5)
#vis_circos(ova_gene_J, .by = "Antigen", .meta = TCR_OVA_Unique$meta)

chordDiagram(ova_gene)
Test_circ = data.frame(from = rep(rownames(ova_gene), time = ncol(ova_gene)), to = rep(colnames(ova_gene), each = nrow(ova_gene)), value = as.factor(ova_gene), stringsAsFactors =  False)

#FVIII
FVIII_gene <- geneUsage(TCR_FVIII_unique$data, "hs.trbv", .norm = T, .ambig = "maj")
TCR_FVIII_unique$data[[1]][TCR_FVIII_unique$data[[1]] == "TRBV7-4"] <- "TRBV20-1"
TCR_FVIII_unique$data[[1]][TCR_FVIII_unique$data[[1]] == "TRBV7-3"] <- "TRBV3-1"
TCR_FVIII_unique$data[[1]][TCR_FVIII_unique$data[[1]] == "TRBV4-1"] <- "TRBV3-1"
TCR_FVIII_unique$data[[2]][TCR_FVIII_unique$data[[2]] == "TRBV7-4"] <- "TRBV20-1"
TCR_FVIII_unique$data[[2]][TCR_FVIII_unique$data[[2]] == "TRBV7-3"] <- "TRBV3-1"
TCR_FVIII_unique$data[[2]][TCR_FVIII_unique$data[[2]] == "TRBV4-1"] <- "TRBV3-1"
TCR_FVIII_unique$data[[3]][TCR_FVIII_unique$data[[3]] == "TRBV7-4"] <- "TRBV3-1"
TCR_FVIII_unique$data[[3]][TCR_FVIII_unique$data[[3]] == "TRBV7-3"] <- "TRBV20-1"
TCR_FVIII_unique$data[[3]][TCR_FVIII_unique$data[[3]] == "TRBV4-1"] <- "TRBV3-1"
TCR_FVIII_unique$data[[4]][TCR_FVIII_unique$data[[4]] == "TRBV7-4"] <- "TRBV3-1"
TCR_FVIII_unique$data[[4]][TCR_FVIII_unique$data[[4]] == "TRBV7-3"] <- "TRBV3-1"
TCR_FVIII_unique$data[[4]][TCR_FVIII_unique$data[[4]] == "TRBV4-1"] <- "TRBV3-1"
TCR_FVIII_unique$data[[5]][TCR_FVIII_unique$data[[5]] == "TRBV7-4"] <- "TRBV3-1"
TCR_FVIII_unique$data[[5]][TCR_FVIII_unique$data[[5]] == "TRBV7-3"] <- "TRBV3-1"
TCR_FVIII_unique$data[[5]][TCR_FVIII_unique$data[[5]] == "TRBV4-1"] <- "TRBV20-1"
TCR_FVIII_unique$data[[6]][TCR_FVIII_unique$data[[6]] == "TRBV7-4"] <- "TRBV20-1"
TCR_FVIII_unique$data[[6]][TCR_FVIII_unique$data[[6]] == "TRBV7-3"] <- "TRBV20-1"
TCR_FVIII_unique$data[[7]][TCR_FVIII_unique$data[[7]] == "TRBV7-4"] <- "TRBV20-1"
TCR_FVIII_unique$data[[7]][TCR_FVIII_unique$data[[7]] == "TRBV7-3"] <- "TRBV20-1"
TCR_FVIII_unique$data[[8]][TCR_FVIII_unique$data[[8]] == "TRBV7-4"] <- "TRBV20-1"
TCR_FVIII_unique$data[[8]][TCR_FVIII_unique$data[[8]] == "TRBV7-3"] <- "TRBV20-1"
TCR_FVIII_unique$data[[9]][TCR_FVIII_unique$data[[9]] == "TRBV7-4"] <- "TRBV20-1"
TCR_FVIII_unique$data[[9]][TCR_FVIII_unique$data[[9]] == "TRBV7-3"] <- "TRBV20-1"
TCR_FVIII_unique$data[[10]][TCR_FVIII_unique$data[[10]] == "TRBV7-4"] <- "TRBV20-1"
TCR_FVIII_unique$data[[10]][TCR_FVIII_unique$data[[10]] == "TRBV7-3"] <- "TRBV20-1"
FVIII_gene <- geneUsage(TCR_FVIII_unique$data, "hs.trbv", .norm = T, .ambig = "maj")

vis(FVIII_gene, .by = "Antigen", .meta = TCR_FVIII_unique$meta, .plot = "hist") + scale_fill_brewer(palette = "Greens")



TCR_FVIII_unique$data[[1]][TCR_FVIII_unique$data[[1]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[2]][TCR_FVIII_unique$data[[2]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[2]][TCR_FVIII_unique$data[[2]] == "TRBJ2-7"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[3]][TCR_FVIII_unique$data[[3]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[3]][TCR_FVIII_unique$data[[3]] == "TRBJ2-7"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[4]][TCR_FVIII_unique$data[[4]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[5]][TCR_FVIII_unique$data[[5]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[5]][TCR_FVIII_unique$data[[5]] == "TRBJ2-7"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[6]][TCR_FVIII_unique$data[[6]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[7]][TCR_FVIII_unique$data[[7]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[8]][TCR_FVIII_unique$data[[8]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[9]][TCR_FVIII_unique$data[[9]] == "TRBJ2-7"] <- "TRBJ1-2"
#TCR_FVIII_unique$data[[9]][TCR_FVIII_unique$data[[9]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_FVIII_unique$data[[10]][TCR_FVIII_unique$data[[10]] == "TRBJ2-1"] <- "TRBJ1-2"
FVIII_gene_J <- geneUsage(TCR_FVIII_unique$data, "hs.trbj", .norm = T, .ambig = "maj")

vis(FVIII_gene_J, .by = "Antigen", .meta = TCR_FVIII_unique$meta, .plot = "hist") + scale_fill_brewer(palette = "Greens")



#POOL

#POOL_gene <- geneUsage(TCR_POOL_unique$data, "hs.trbv", .norm = T, .ambig = "maj")
#vis(POOL_gene)
#FVIII_gene <- geneUsage(TCR_POOL_unique$data, "hs.trbv", .norm = T, .ambig = "maj")
TCR_POOL_unique$data[[1]][TCR_POOL_unique$data[[1]] == "TRBV7-4"] <- "TRBV20-1"
TCR_POOL_unique$data[[1]][TCR_POOL_unique$data[[1]] == "TRBV7-3"] <- "TRBV3-1"
TCR_POOL_unique$data[[1]][TCR_POOL_unique$data[[1]] == "TRBV4-1"] <- "TRBV20-1"
TCR_POOL_unique$data[[2]][TCR_POOL_unique$data[[2]] == "TRBV7-4"] <- "TRBV3-1"
TCR_POOL_unique$data[[2]][TCR_POOL_unique$data[[2]] == "TRBV7-3"] <- "TRBV20-1"
TCR_POOL_unique$data[[2]][TCR_POOL_unique$data[[2]] == "TRBV4-1"] <- "TRBV3-1"
TCR_POOL_unique$data[[3]][TCR_POOL_unique$data[[3]] == "TRBV7-4"] <- "TRBV20-1"
TCR_POOL_unique$data[[3]][TCR_POOL_unique$data[[3]] == "TRBV7-3"] <- "TRBV3-1"
TCR_POOL_unique$data[[3]][TCR_POOL_unique$data[[3]] == "TRBV4-1"] <- "TRBV20-1"
TCR_POOL_unique$data[[4]][TCR_POOL_unique$data[[4]] == "TRBV7-4"] <- "TRBV3-1"
TCR_POOL_unique$data[[4]][TCR_POOL_unique$data[[4]] == "TRBV7-3"] <- "TRBV20-1"
TCR_POOL_unique$data[[4]][TCR_POOL_unique$data[[4]] == "TRBV4-1"] <- "TRBV3-1"
TCR_POOL_unique$data[[5]][TCR_POOL_unique$data[[5]] == "TRBV7-4"] <- "TRBV20-1"
TCR_POOL_unique$data[[5]][TCR_POOL_unique$data[[5]] == "TRBV5-6-1"] <- "TRBV3-1"
POOL_gene <- geneUsage(TCR_POOL_unique$data, "hs.trbv", .norm = T, .ambig = "maj")

vis(POOL_gene, .by = "Antigen", .meta = TCR_POOL_unique$meta, .plot = "hist") + scale_fill_brewer(palette = "Blues")



TCR_POOL_unique$data[[1]][TCR_POOL_unique$data[[1]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_POOL_unique$data[[2]][TCR_POOL_unique$data[[2]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_POOL_unique$data[[2]][TCR_POOL_unique$data[[2]] == "TRBJ2-7"] <- "TRBJ1-2"
TCR_POOL_unique$data[[3]][TCR_POOL_unique$data[[3]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_POOL_unique$data[[3]][TCR_POOL_unique$data[[3]] == "TRBJ2-7"] <- "TRBJ1-2"
TCR_POOL_unique$data[[4]][TCR_POOL_unique$data[[4]] == "TRBJ2-1"] <- "TRBJ1-2"
TCR_POOL_unique$data[[5]][TCR_POOL_unique$data[[5]] == "TRBJ2-1"] <- "TRBJ1-1"
TCR_POOL_unique$data[[5]][TCR_POOL_unique$data[[5]] == "TRBJ2-7"] <- "TRBJ1-1"



POOL_gene_J <- geneUsage(TCR_POOL_unique$data, "hs.trbj", .norm = T, .ambig = "maj")

vis(POOL_gene_J, .by = "Antigen", .meta = TCR_POOL_unique$meta, .plot = "hist") + scale_fill_brewer(palette = "BuPu")










#######################################################Circulize plot



#####FVIII & OVA





#VJ gene usage        ###################### FVIII

chordDiagram(Matrice_VDJ_F8_F8, annotationTrack = "grid", preAllocateTracks =  1)


# col_fun = colorRamp2(c(20, 100, 200), c("blue", "white", "red"))
# col_fun(seq(-5, 1, by = 1))
# 




# chordDiagram(Matrice_VDJ_F8_F8, grid.col = c("#99CC99"), annotationTrack = "grid", preAllocateTracks = 1)
# chordDiagram(Matrice_VDJ_F8_F8, grid.col = c("#99CC99"))
# circos.text(10, 5, "clockwise", facing = "clockwise", adj = c(0.5, 0),
#             cex = 0.8)
# chordDiagram(Matrice_VDJ_F8_F8, annotationTrack = "grid", preAllocateTracks = 1)



circos.trackPlotRegion(track.index = 2, panel.fun = function(x,y){
  xlim=get.cell.meta.data("xlim")
  ylim=get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
})




######SANS VALEURS AXES
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  #print labels 
  circos.text(mean(xlim), ylim[1] + 1.5, sector.name, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.5)
  
}, bg.border = NA)

col = c('TRBV10-3'="#99CC00",'TRBV11-2'="#CCFF00",'TRBV16'="#99FF00",
        'TRBV20-1'="#66CC66", 'TRBV25-1'="#33FF00",'TRBV3-1'="#00CC00",
        'TRBV6-3'="#33CC33",'TRBV7-2'="#006600", 'TRBV7-3'="#336633",
        'TRBV7-4'="#66CC66",'TRBV7-8'="#33CC00",'TRBV7-9'="#66FF33",
        'TRBV4-1'="#99FF66",'TRBV10-2'="#CCFF66",'TRBV11-3'="#99CC33",
        'TRBV28'="#66CC33", 'TRBV17'="#339900",'TRBV29-1'="#336600",
        'TRBV14'="#00FF00",'TRBV27'="#33FF33", 'TRBV5-4'="#669933",
        'TRBV12-3'="#669900",'TRBV21-1'="#666633",'TRBV24-1'="#666600",
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
        'TRBJ2-6'="#00FF33",'TRBJ2-7'="#003300")

chordDiagram(Matrice_VDJ_F8_F8, grid.col = col, annotationTrack = "grid", preAllocateTracks =  1)



#Saving plots
#dev.copy(jpeg,'F8_circular_VDJ.png', width=8, height=8, units="in", res=500)


#####AVEC VALEURS AXES
# circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
#   xlim = get.cell.meta.data("xlim")
#   ylim = get.cell.meta.data("ylim")
#   sector.name = get.cell.meta.data("sector.index")
#   
#   #print labels 
#   circos.text(mean(xlim), ylim[1] + 2.5, sector.name, 
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.6)
#   
#   #print axis
#   circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, 
#               sector.index = sector.name, track.index = 2)
# }, bg.border = NA)

#############################################################################
#############################################################################
#############################################################################
#############################################################################



#VJ gene usage #################### OVA

chordDiagram(Matrice_VDJ_OVA, annotationTrack = "grid", preAllocateTracks =  1)


# col_fun = colorRamp2(c(20, 100, 200), c("blue", "white", "red"))
# col_fun(seq(-5, 1, by = 1))
# 




# chordDiagram(Matrice_VDJ_OVA, grid.col = c("#99CC99"), annotationTrack = "grid", preAllocateTracks = 1)
# chordDiagram(Matrice_VDJ_OVA, grid.col = c("#99CC99"))
# circos.text(10, 5, "clockwise", facing = "clockwise", adj = c(0.5, 0),
#             cex = 0.8)
# chordDiagram(Matrice_VDJ_OVA, annotationTrack = "grid", preAllocateTracks = 1)



circos.trackPlotRegion(track.index = 2, panel.fun = function(x,y){
  xlim=get.cell.meta.data("xlim")
  ylim=get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
})




######SANS VALEURS AXES
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  
  #print labels 
  circos.text(mean(xlim), ylim[1] + 1.5, sector.name, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.5)
  
}, bg.border = NA)

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

chordDiagram(Matrice_VDJ_OVA, grid.col = col, annotationTrack = "grid", preAllocateTracks =  1)



#Saving plots
#dev.copy(jpeg,'F8_circular_VDJ.png', width=8, height=8, units="in", res=500)


#####AVEC VALEURS AXES
# circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
#   xlim = get.cell.meta.data("xlim")
#   ylim = get.cell.meta.data("ylim")
#   sector.name = get.cell.meta.data("sector.index")
#   
#   #print labels 
#   circos.text(mean(xlim), ylim[1] + 2.5, sector.name, 
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.6)
#   
#   #print axis
#   circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, 
#               sector.index = sector.name, track.index = 2)
# }, bg.border = NA)

