RNA <- read.csv("~/Pulpit/analiza_statystyczna/RNA/counts_rna.csv")
DE <- read.csv("~/Pulpit/analiza_statystyczna/korelacja/nowa_korelacja/analiza/Najróżniej_ekspresjonujace_geny.csv")
RRBS <- read.csv("~/Pulpit/analiza_statystyczna/korelacja/counts_rrbs.csv")

##################################### przygotowanie danych ########################################################
#geny różnie ekspresjonowane
DE_RNA <- merge(DE,RNA,by="Geneid",all.x=T)
DE_RNA <- DE_RNA[,-c(2:6)]

#normalizacja
rownames(DE_RNA) <- DE_RNA[,1]
DE_RNA <- DE_RNA[,-1]
counts <- as.matrix(DE_RNA)
counts <- apply(counts, c(1,2), function(x) log(x+10))

#grupowanie
RRBS <- RRBS[,-1]
RRBS = aggregate(RRBS,
                by = list(RRBS$Geneid),
                FUN = mean)
RRBS <- replace(RRBS, is.na(RRBS), 0)

#usuniecie kolumn gdzie jest tylko 0
RRBS <- RRBS[, colSums(RRBS != 0) > 0]

#połączenie DE_RNA z metylacją
rownames(RRBS) <- RRBS[,1]
RRBS <- RRBS[,-1]
methyl <- as.matrix(RRBS)
DE_METHYL <- merge(counts,methyl,by="row.names",all.x=TRUE)
DE_METHYL <- replace(DE_METHYL, is.na(DE_METHYL), 0)

#usuniecie wierszy, gdzie metylacja byla 0 w kazdej probie
cleaned_DE_METHYL = DE_METHYL[ rowSums(DE_METHYL[,c(59:115)])!=0, ] 

write.csv(cleaned_DE_METHYL,"~/Pulpit/metylacja_ekspresja.csv")


################################### analiza korelacji #################################################

#przygotowanie danych
rownames(cleaned_DE_METHYL) <- cleaned_DE_METHYL[,1]
cleaned_DE_METHYL <- cleaned_DE_METHYL[,-1]

#transpozycja
data <- t(cleaned_DE_METHYL)

#rozdzielenie danych na rrbs i rna
rna <- data[1:57,]
names_rna <- colnames(rna)
rrbs <- data[58:114,]
names_rrbs <- colnames(rrbs)


#odpowiednie uporządkowanie wierszy 
target_rrbs <-substr(rownames(rrbs)[],1,10)
target_rna <- substr(rownames(rna)[],1,10)

rna <- cbind(rna,target_rna)
rna <- rna[order(factor(target_rna,levels=target_rrbs)),]

rrbs <- cbind(rrbs,target_rrbs)
rrbs<-rrbs[order(factor(target_rrbs,levels=target_rrbs)),]


#obliczenie korelacji
cor <- c()

rna_numeric <- mapply(rna, FUN=as.numeric)
rna_numeric <- matrix(data=rna_numeric, ncol=dim(rna)[2], nrow=dim(rna)[1])

rrbs_numeric <- mapply(rrbs, FUN=as.numeric)
rrbs_numeric <- matrix(data=rrbs_numeric, ncol=dim(rrbs)[2], nrow=dim(rrbs)[1])

for(i in 1:202)
{
  cor[i] <- cor.test(rna_numeric[,i],rrbs_numeric[,i], methon="spearman")$estimate
}


################################### wykresy ######################################

# jakie wartości przyjmuje korelacja 
plot(cor)

result <- rbind(rrbs,rna)
result <- rbind(result,cor)

