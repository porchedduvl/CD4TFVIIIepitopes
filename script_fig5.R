
#Graph pour FVIII&Pool_Unique
TCR_FVIII_POOL_unique <- repLoad(.path = "/mnt/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/_Analyse_fichiers_corriges_FVIII&POOL_unique")

public_rep_1 <- pubRep(
  TCR_FVIII_POOL_unique$data  ,
  .col = "aa",
  .quant = c("prop"),
  .coding = TRUE,
  .min.samples = 1,
  .max.samples = NA,
  .verbose = TRUE
)

public_rep_2F <- pubRepFilter(public_rep_1, TCR_FVIII_POOL_unique$meta, .by = c(Antigen = "FVIII"))

vis(public_rep_2F, "freq", .type = "none")


vis(public_rep_1, "freq" , .type = "none" ) +geom_point() + xlim(0.01, )

public_1 <- subset(public_rep_1, Samples > 1)
public_rep_1[is.na(public_rep_1)]=0
vis(public_1, "freq")

write.table(public_1, "Public_rep_FVIII_POOL.csv", sep=";",row.names=FALSE)


#Graph pour FVIII&Pool_FULL
TCR_FVIII_POOL_Full <- repLoad(.path = "/mnt/NGS_YeastDisplay/Pipeline_NGS_TCR/_Valeria/_Analyse_fichiers_corriges_FVIII&POOL_unique")

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

p <- ggplot(df_pub_rep_1_plot, aes(x = Donneur, y = Value, color = as.factor(Samples))) + geom_point()+xlab("Samples")+ylab("Proportion of clonotypes")
p+theme_light()


##########################################################################################################
####FIGURE 5D

df_pub_rep_2_plot <- data.frame(Value = c(unlist(public_rep_1$RUN_2_1238_FVIII),
                                          unlist(public_rep_1$RUN_2_1245_FVIII),
                                          unlist(public_rep_1$RUN_2_1255_FVIII),
                                          unlist(public_rep_1$RUN_2_1285_FVIII),
                                          unlist(public_rep_1$RUN_3_1327_FVIII)))

df_pub_rep_2_plot$Donneur <- c(rep.int("1238", 19347),
                               rep.int("1245", 19347),
                               rep.int("1255", 19347),
                               rep.int("1285", 19347),
                               rep.int("1327", 19347))

df_pub_rep_2_plot$Samples <- rep(public_rep_1$Samples,5)
df_pub_rep_2_plot$Seq <- rep(public_rep_1$CDR3.aa,5)
#on vire toutes les s?quences qui ne sont pas pr?sentes chez 2 samples
df_public_FVIII_POOL_rep <- subset(df_pub_rep_2_plot, Samples > 1)
view(df_public_FVIII_POOL_rep)



p <- ggplot(df_public_FVIII_POOL_rep, aes(x = Samples, y = Value, color = as.factor(Donneur))) + geom_jitter(size=1, width = 0.1, height = 0.1) +xlab("Number of donors sharing the clonotypes")+ylab("Proportion of public clonotypes")
p+theme_light()+ggtitle("FVIII")+theme(plot.title = element_text(hjust = 0.5))+scale_colour_discrete("Donor") 


p <- ggplot(df_public_FVIII_POOL_rep, aes(x = Samples, y = Value, color = as.factor(Donneur))) + geom_point(size=3, position=position_jitterdodge(dodge=0.5)) +xlab("Number of donors sharing the clonotypes")+ylab("Proportion of public clonotypes")
p+theme_light()+ggtitle("FVIII")+theme(plot.title = element_text(hjust = 0.5))+scale_colour_discrete("Donor") 
#####################################################################


library("immunarch")

p3 <- ggplot(df_public_FVIII_POOL_rep, aes(x = Donneur, y = Value, color = as.factor(Samples))) + geom_point()+xlab("Samples")+ylab("Proportion of public clonotypes")
p3+theme_light()+ggtitle("FVIII")+theme(plot.title = element_text(hjust = 0.5))+scale_colour_discrete("Number of donors sharing the clonotypes")

p3 <- ggplot(df_public_FVIII_POOL_rep, aes(x = Samples, y = Value, color = as.factor(Donneur))) + geom_point()+xlab("Donor")+ylab("Proportion of public clonotypes")
p3+theme_light()+ggtitle("FVIII")+theme(plot.title = element_text(hjust = 0.5))+scale_colour_discrete("Number of donors sharing the clonotypes")

LC <- public_rep_2$CDR3.aa[1:96]

p4 <- trackClonotypes(TCR_FVIII_POOL_unique$data[c(2,4,6,8,10)], public_rep_2$CDR3.aa[public_rep_2$Samples > 3])
vis(p4) + theme(legend.position = "none")

gaga <- public_rep_2
for (ligne in 1:nrow(public_rep_1)) {
  if (gaga$RUN_2_1238_FVIII_POOL[ligne] != 0 &
      gaga$RUN_2_1238_FVIII[ligne] != 0)  {
    gaga$Samples[ligne] <- gaga$Samples[ligne] - 1
  }
  if (gaga$RUN_2_1245_FVIII_POOL[ligne] != 0 &
      gaga$RUN_2_1245_FVIII[ligne] != 0)  {
    gaga$Samples[ligne] <- gaga$Samples[ligne] - 1
  }
  if (gaga$RUN_2_1255_FVIII_POOL[ligne] != 0 &
      gaga$RUN_2_1255_FVIII[ligne] != 0)  {
    gaga$Samples[ligne] <- gaga$Samples[ligne] - 1
  }
  if (gaga$RUN_2_1285_FVIII_POOL[ligne] != 0 &
      gaga$RUN_2_1285_FVIII[ligne] != 0)  {
    gaga$Samples[ligne] <- gaga$Samples[ligne] - 1
  }
  if (gaga$RUN_3_1327_FVIII_POOL[ligne] != 0 &
      gaga$RUN_3_1327_FVIII[ligne] != 0)  {
    gaga$Samples[ligne] <- gaga$Samples[ligne] - 1
    }
}


df_gaga <- data.frame(Value = c(unlist(gaga$RUN_2_1238_FVIII),
                                          unlist(gaga$RUN_2_1245_FVIII),
                                          unlist(gaga$RUN_2_1255_FVIII),
                                          unlist(gaga$RUN_2_1285_FVIII),
                                          unlist(gaga$RUN_3_1327_FVIII)))

df_gaga$Donneur <- c(rep.int("1238", 19347),
                               rep.int("1245", 19347),
                               rep.int("1255", 19347),
                               rep.int("1285", 19347),
                               rep.int("1327", 19347))

df_gaga$Samples <- rep(gaga$Samples,5)
df_gaga$Seq <- rep(gaga$CDR3.aa,5)


p3 <- ggplot(df_gaga, aes(x = Samples, y = Value, color = as.factor(Donneur))) + geom_point()+xlab("Donor")+ylab("Proportion of public clonotypes")
p3+theme_light()+ggtitle("FVIII")+theme(plot.title = element_text(hjust = 0.5))+scale_colour_discrete("Number of donors sharing the clonotypes")



df_gaga2 <- subset(df_gaga, Samples > 1)
p4 <- trackClonotypes(TCR_FVIII_POOL_unique$data[c(2,4,6,8,10)], df_gaga2$Seq)
vis(p4) + theme(legend.position = "left")
