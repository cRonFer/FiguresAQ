setwd("C:/Users/Usuario/Desktop/Sci Rev2 Figures/Materials Fig3")
setwd("C:/Users/Usuario/Desktop/GLOBI")

Seqs <- read.csv("terrentalife_seqs.csv")
Indi <- read.csv("terrentalife_BBDD.csv")

Merge <- merge(Seqs, Indi, by="SpeciesAccepted")

#Merge[Merge$Kingdom.x == "Amoebozoa",] <- "SAR/Protist"
Merge[1652,6] <- "SAR/Protist"

vectorreinos <- as.data.frame(unique(Merge$Kingdom.x))

Mat<-matrix(nrow=length(unique(Merge$Kingdom.x)),ncol=7)
rownames(Mat) <- unique(Merge$Kingdom.x)
colnames(Mat) <- c("Kingdom", "TotalSp", "SpsWithData", "BinaryFraction", "TotalSeqs", "Average", "Kingdom")

library(ggplot2)


for (i in 1:length(vectorreinos[,1])) {
  tryCatch({
a <- Merge[Merge$Kingdom.x == paste(vectorreinos[i,1]),]
b <- nrow(a)
c <- sum(a$info_numb > 0, na.rm=TRUE)
d <- c/b
e <- sum(as.numeric(a$info_numb, na.omit=TRUE))
f <- mean(as.numeric(a$info_numb))
Mat[i,1] <- paste(vectorreinos[i,1])
Mat[i,2] <- b
Mat[i,3] <- c
Mat[i,4] <- d
Mat[i,5] <- e
Mat[i,6] <- f
Mat[i,7] <- paste(unique(a$Kingdom.x))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

Mat2 <- as.data.frame(Mat)

colnames(Mat2)[1] <- "Kingdom2"

Mat2$BinaryFraction <- as.numeric(Mat2$BinaryFraction)


#Mat3 <- Mat2[Mat2$Kingdom == "Fungi",]


#############BY FUNCTIONAL GROUP


setwd("C:/Users/Usuario/Desktop/Sci Rev2 Figures/Materials Fig3")

Seqs <- read.csv("terrentalife_seqs.csv")
Indi <- read.csv("terrentalife_BBDD.csv")

Merge <- merge(Seqs, Indi, by="SpeciesAccepted")

#Merge[Merge$Kingdom.x == "Amoebozoa",] <- "SAR/Protist"
#Merge[1652,6] <- "SAR/Protist"

vectorfunc <- as.data.frame(unique(Merge$Func_group))

Mat<-matrix(nrow=length(unique(Merge$Func_group)),ncol=7)
rownames(Mat) <- unique(Merge$Func_group)
colnames(Mat) <- c("FuncGr", "TotalSp", "SpsWithData", "BinaryFraction", "TotalSeqs", "Average", "Kingdom")

library(ggplot2)



for (i in 1:length(vectorfunc[,1])) {
  tryCatch({
    Da <- Merge[Merge$Func_group == paste(vectorfunc[i,1]),]
    Db <- nrow(Da)
    Dc <- sum(Da$info_numb > 0, na.rm=TRUE)
    Dd <- Dc/Db
    De <- sum(as.numeric(Da$info_numb, na.omit=TRUE))
    Df <- mean(as.numeric(Da$info_numb))
    Mat[i,1] <- paste(vectorfunc[i,1])
    Mat[i,2] <- Db
    Mat[i,3] <- Dc
    Mat[i,4] <- Dd
    Mat[i,5] <- De
    Mat[i,6] <- Df
    Mat[i,7] <- paste(names(which.max(table(Da$Kingdom.x))))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

Mat2 <- as.data.frame(Mat)


Mat2$BinaryFraction <- as.numeric(Mat2$BinaryFraction)
Mat2$TotalSeqs <- as.numeric(Mat2$TotalSeqs)
Mat2$Average <- as.numeric(Mat2$Average)

#Mat3 <- Mat2[Mat2$Kingdom == "Fungi",]


Mat2SeqFrac <- Mat2 %>% 
  arrange(Kingdom, BinaryFraction) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = BinaryFraction, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none") +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Fraction of sp. with DNA sequences")+
  xlab("")+
  theme(axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) #+ scale_fill_hue(c=50, l=50)

Mat2SeqTot <- Mat2 %>% 
  arrange(Kingdom, TotalSeqs) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = TotalSeqs, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none") +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Total DNA sequences")+
  xlab("")+
  theme(axis.text.x = element_blank())  + scale_y_sqrt()
  

Mat2SeqAvg <- Mat2 %>% 
  arrange(Kingdom, Average) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = Average, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none") +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Average DNA sequences")+
  xlab("")+
  theme(axis.text.x = element_blank())  + scale_y_log10(breaks=c(10,100,1000,10000),labels=c(10,100,1000,10000))

#scale_color_manual(values = c ('blue', 'red', 'yellow1', 'green'))

barplot(Mat2$BinaryFraction, names.arg = Mat2$FuncGr, ylim = c( 0 , 1 ))

write.table(Mat2, "Genetic coverage.csv")


#############FOR SPATIAL OCCURRENCE DATA

setwd("C:/Users/Usuario/Desktop/Sci Rev2 Figures/Materials Fig3")

Seqs <- read.csv("terrentalife_seqs.csv")
Indi <- read.csv("terrentalife_BBDD.csv")

OccData <- read.csv("Occ_Data_incOcean.csv", sep=";")

OccData$sppName[OccData$sppName == "Pygoscelis antarctica"] <- "Pygoscelis antarcticus"

sum(OccData$sppName == "Pygoscelis antarctica")
sum(OccData$sppName == "Pygoscelis antarcticus")

#for spatially filtered records
library(dplyr)
OccDataCount <- OccData %>% 
  select(sppName, lat_round, lon_round) %>% 
  unique() %>% 
  group_by(sppName) %>% 
  summarize(n = n())

###for all raw records
library(dplyr)
OccDataCount2 <- OccData %>% 
  select(sppName, Lat) %>% 
  group_by(sppName) %>% 
  summarize(n = n())

colnames(OccDataCount)[1] <- "SpeciesAccepted"

MergeXY <- merge(Indi, OccDataCount, by="SpeciesAccepted", all.x=TRUE)

MergeXY$n[is.na(MergeXY$n)] <- 0

#Merge[Merge$Kingdom.x == "Amoebozoa",] <- "SAR/Protist"
#Merge[1652,6] <- "SAR/Protist"

vectorfuncXY <- as.data.frame(unique(MergeXY$Func_group))

MatXY<-matrix(nrow=length(unique(MergeXY$Func_group)),ncol=7)
rownames(MatXY) <- unique(MergeXY$Func_group)
colnames(MatXY) <- c("FuncGr", "TotalSp", "SpsWithData", "BinaryFraction", "TotalRecords", "Average", "Kingdom")

library(ggplot2)

for (i in 1:length(vectorfuncXY[,1])) {
  tryCatch({
    Wa <- MergeXY[MergeXY$Func_group == paste(vectorfuncXY[i,1]),]
    Wb <- nrow(Wa)
    Wc <- sum(Wa$n > 0, na.rm=TRUE)
    Wd <- Wc/Wb
    We <- sum(as.numeric(Wa$n, na.omit=TRUE))
    Wf <- mean(as.numeric(Wa$n))
    MatXY[i,1] <- paste(vectorfuncXY[i,1])
    MatXY[i,2] <- Wb
    MatXY[i,3] <- Wc
    MatXY[i,4] <- Wd
    MatXY[i,5] <- We
    MatXY[i,6] <- Wf
    MatXY[i,7] <- paste(names(which.max(table(Wa$Kingdom.x))))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

MatXY2 <- as.data.frame(MatXY)


MatXY2$BinaryFraction <- as.numeric(MatXY2$BinaryFraction)
MatXY2$TotalRecords <- as.numeric(MatXY2$TotalRecords)
MatXY2$Average <- as.numeric(MatXY2$Average)

#Mat3 <- Mat2[Mat2$Kingdom == "Fungi",]

#p<-ggplot(data=Mat2, aes(x=Mat2$FuncGr, y=Mat2$BinaryFraction)) +
 # geom_bar(stat="identity", fill) + scale_y_continuous(breaks = seq(0, 1, 100000)) + 
  #  p

MatXY2Frac <- MatXY2 %>% 
  arrange(Kingdom, BinaryFraction) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = BinaryFraction, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none") +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Fraction of sp. with spatial records")+
  xlab("")+
  theme(axis.text.x = element_blank()) +scale_y_sqrt()

MatXY2Tot <-MatXY2 %>% 
  arrange(Kingdom, TotalRecords) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = TotalRecords, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none") +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Total spatial records")+
  xlab("")+
  theme(axis.text.x = element_blank()) +scale_y_sqrt()

MatXY2Avg <- MatXY2  %>% 
  arrange(Kingdom, Average) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = Average, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none") +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Average spatial records")+
  xlab("")+   theme(axis.text.x = element_blank()) +
  scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000),labels=c(10,100,1000,10000,100000,1000000))


#scale_color_manual(values = c ('blue', 'red', 'yellow1', 'green'))

#barplot(Mat2$BinaryFraction, names.arg = Mat2$FuncGr, ylim = c( 0 , 1 ))

write.table(Mat2, "Genetic coverage.csv")



#############FOR ELTONIAN INTERACTION DATA

setwd("C:/Users/Usuario/Desktop/GLOBI")
#species <- read.delim("C:/Users/Usuario/Desktop/GLOBI/terraANTALIFE_eukariotic_v1.txt")
Interact <- read.csv("C:/Users/Usuario/Desktop/GLOBI/fullareaSPdata.csv")
Indi <- read.csv("terrentalife_BBDD.csv")

#OccData$sppName[OccData$sppName == "Pygoscelis antarctica"] <- "Pygoscelis antarcticus"

#sum(OccData$sppName == "Pygoscelis antarctica")
#sum(OccData$sppName == "Pygoscelis antarcticus")


library(dplyr)
InteractCount <- Interact %>% 
  select(source_taxon_name) %>% 
  group_by(source_taxon_name) %>% 
  summarize(n = n())

colnames(InteractCount)[1] <- "SpeciesAccepted"

MergeElt <- merge(Indi, InteractCount, by="SpeciesAccepted", all.x=TRUE)

MergeElt$n[is.na(MergeElt$n)] <- 0

#Merge[Merge$Kingdom.x == "Amoebozoa",] <- "SAR/Protist"
#Merge[1652,6] <- "SAR/Protist"

vectorfuncElt <- as.data.frame(unique(MergeElt$Func_group))

MatElt<-matrix(nrow=length(unique(MergeElt$Func_group)),ncol=7)
rownames(MatElt) <- unique(MergeElt$Func_group)
colnames(MatElt) <- c("FuncGr", "TotalSp", "SpsWithData", "BinaryFraction", "TotalRecords", "Average", "Kingdom")

library(ggplot2)

for (i in 1:length(vectorfuncElt[,1])) {
  tryCatch({
    Wa <- MergeElt[MergeElt$Func_group == paste(vectorfuncElt[i,1]),]
    Wb <- nrow(Wa)
    Wc <- sum(Wa$n > 0, na.rm=TRUE)
    Wd <- Wc/Wb
    We <- sum(as.numeric(Wa$n, na.omit=TRUE))
    Wf <- mean(as.numeric(Wa$n))
    MatElt[i,1] <- paste(vectorfuncElt[i,1])
    MatElt[i,2] <- Wb
    MatElt[i,3] <- Wc
    MatElt[i,4] <- Wd
    MatElt[i,5] <- We
    MatElt[i,6] <- Wf
    MatElt[i,7] <- paste(names(which.max(table(Wa$Kingdom.x))))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

MatElt2 <- as.data.frame(MatElt)


MatElt2$BinaryFraction <- as.numeric(MatElt2$BinaryFraction)
MatElt2$TotalRecords <- as.numeric(MatElt2$TotalRecords)
MatElt2$Average <- as.numeric(MatElt2$Average)

#Mat3 <- Mat2[Mat2$Kingdom == "Fungi",]

#p<-ggplot(data=Mat2, aes(x=Mat2$FuncGr, y=Mat2$BinaryFraction)) +
# geom_bar(stat="identity", fill) + scale_y_continuous(breaks = seq(0, 1, 100000)) + 
#  p

MatElt2Frac <-MatElt2 %>% 
  arrange(Kingdom, BinaryFraction) %>%
  mutate(FuncGr=factor(FuncGr, levels = FuncGr), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = BinaryFraction, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none")  +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Fraction of sp. with interactions data")+
  xlab("")+
  theme(axis.text.x = element_blank()) +scale_y_sqrt()


MatElt2Tot <- MatElt2 %>% 
  arrange(Kingdom, TotalRecords) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = TotalRecords, fill = Kingdom)) +
  geom_bar(stat = "identity") + theme(legend.position = "none") +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Total interactions")+
  xlab("")+ 
  theme(axis.text.x = element_blank()) +scale_y_sqrt()

MatElt2Avg <- MatElt2 %>% 
  arrange(Kingdom, Average) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = Average, fill = Kingdom)) + 
  geom_bar(stat = "identity")+ theme(legend.position = "none") + 
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Average interactions")+
  xlab("")+
  theme(axis.text.x = element_blank()) +scale_y_sqrt()
 
write.csv(MatElt2, "MatElt2.csv")

#############FOR HUTCHINSON DATA DATA

setwd("C:/Users/Usuario/Desktop/GLOBI")
#species <- read.delim("C:/Users/Usuario/Desktop/GLOBI/terraANTALIFE_eukariotic_v1.txt")
Interact <- read.csv("C:/Users/Usuario/Desktop/GLOBI/fullareaSPdata.csv")

Hutchi <- read_excel("Hutchinsonian.xlsx")
Indi <- read.csv("terrentalife_BBDD.csv")

#OccData$sppName[OccData$sppName == "Pygoscelis antarctica"] <- "Pygoscelis antarcticus"

#sum(OccData$sppName == "Pygoscelis antarctica")
#sum(OccData$sppName == "Pygoscelis antarcticus")


library(dplyr)
HutchiCount <- Hutchi %>% 
  select(SpeciesAccepted) %>% 
  group_by(SpeciesAccepted) %>% 
  summarize(n = n())

colnames(HutchiCount)[1] <- "SpeciesAccepted"

MergeHutchi <- merge(Indi, HutchiCount, by="SpeciesAccepted", all.x=TRUE)

MergeHutchi$n[is.na(MergeHutchi$n)] <- 0

#Merge[Merge$Kingdom.x == "Amoebozoa",] <- "SAR/Protist"
#Merge[1652,6] <- "SAR/Protist"

vectorfuncHutchi <- as.data.frame(unique(MergeHutchi$Func_group))

MatHutchi<-matrix(nrow=length(unique(MergeHutchi$Func_group)),ncol=7)
rownames(MatHutchi) <- unique(MergeHutchi$Func_group)
colnames(MatHutchi) <- c("FuncGr", "TotalSp", "SpsWithData", "BinaryFraction", "TotalRecords", "Average", "Kingdom")

library(ggplot2)

i=2
i=3
i=4
i=5
i=6
i=7
i=8
  

for (i in 1:length(vectorfuncHutchi)) {
  tryCatch({
    Wa <- MergeHutchi[MergeHutchi$Func_group == paste(vectorfuncHutchi[i,1]),]
    Wb <- nrow(Wa)
    Wc <- sum(Wa$n > 0, na.rm=TRUE)
    Wd <- Wc/Wb
    We <- sum(as.numeric(Wa$n, na.omit=TRUE))
    Wf <- mean(as.numeric(Wa$n))
    MatHutchi[i,1] <- paste(vectorfuncHutchi[i,1])
    MatHutchi[i,2] <- Wb
    MatHutchi[i,3] <- Wc
    MatHutchi[i,4] <- Wd
    MatHutchi[i,5] <- We
    MatHutchi[i,6] <- Wf
    MatHutchi[i,7] <- paste(names(which.max(table(Wa$Kingdom.x))))
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

MatHutchi2 <- as.data.frame(MatHutchi)


MatHutchi2$BinaryFraction <- as.numeric(MatHutchi2$BinaryFraction)
MatHutchi2$TotalRecords <- as.numeric(MatHutchi2$TotalRecords)
MatHutchi2$Average <- as.numeric(MatHutchi2$Average)

#Mat3 <- Mat2[Mat2$Kingdom == "Fungi",]

#p<-ggplot(data=Mat2, aes(x=Mat2$FuncGr, y=Mat2$BinaryFraction)) +
# geom_bar(stat="identity", fill) + scale_y_continuous(breaks = seq(0, 1, 100000)) + 
#  p

levels = c("Megafauna","Soil invertebrate","Parasitic invertebrate",
           "Lichenic symbiont", "Non-lichenized fungi",
           "Embryophytes",
           "Algae","Heterotrophic protist") # customize vector by order

MatHutchi2Frac<- MatHutchi2 %>% 
  arrange(Kingdom, BinaryFraction) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = BinaryFraction, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none") +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Fraction of sp. with thermotolerance data")+
  xlab("")+
  theme(axis.text.x = element_blank()) +scale_y_sqrt()

MatHutchi2Tot <- MatHutchi2 %>% 
  arrange(Kingdom, TotalRecords) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), 
         Kingdom = factor(Kingdom)) %>% 
  ggplot(aes(x = FuncGr, y = TotalRecords, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none") +
  scale_fill_manual(values = c ('#F57A77', '#7CA1CC', '#e7c1a6', '#e5e297'))+
  ylab("Total species thermotolerances data")+
  xlab("")+
  theme(axis.text.x = element_blank()) +scale_y_sqrt()

write.csv(MatHutchi2, "MatHutchi2.csv")


levels = c("Megafauna","Soil invertebrate","Parasitic invertebrate",
           "Lichenic symbiont", "Non-lichenized fungi",
           "Embryophytes",
           "Algae","Heterotrophic protist") # customize vector by order
MatHutchi2Avg <- MatHutchi2 %>%
  arrange(Kingdom, Average) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), # order of bars
         Kingdom = factor(Kingdom)) %>%
  ggplot(aes(x = FuncGr, y = Average, fill = Kingdom)) +
  geom_bar(stat = "identity")+ theme(legend.position = "none") + 
  scale_fill_manual(values = c('#F57A77', '#7CA1CC',
                               '#e7c1a6', '#e5e297'))+
  ylab("Average sp. thermotolerances data") +
  xlab("") +
  theme(axis.text.x = element_blank()) +scale_y_sqrt()

library(patchwork)

plot(Mat2SeqFrac) | plot(MatXY2Frac) | plot(MatElt2Frac) | plot(MatHutchi2Frac)

plot(Mat2SeqTot) | plot(MatXY2Tot) | plot(MatElt2Tot) | plot(MatHutchi2Tot)

plot(Mat2SeqAvg) | plot(MatXY2Avg) | plot(MatElt2Avg) | plot(MatHutchi2Avg)


MatHutchi2AvgL <- MatHutchi2 %>%
  arrange(Kingdom, Average) %>%
  mutate(FuncGr=factor(FuncGr, levels = levels), # order of bars
         Kingdom = factor(Kingdom)) %>%
  ggplot(aes(x = FuncGr, y = Average, fill = Kingdom)) +
  geom_bar(stat = "identity")+  
  scale_fill_manual(values = c('#F57A77', '#7CA1CC',
                               '#e7c1a6', '#e5e297'))+
  ylab("Average sp. thermotolerances data") +
  xlab("") +
  theme(axis.text.x = element_blank()) +scale_y_sqrt()