directory <- getwd()
#załaduj dane
rawCountTable <- read.csv("counts_rna.csv",header=TRUE,sep=',',row.names=2)
#usuń niepotrzebą kolumnę
rawCountTable <- rawCountTable[,-c(0:1)]
#dodaj informacje, które próby pochodzą od osób zdrowych, a które od chorych
info <-read.csv("info.csv",header=TRUE, sep=',',row.names=1)

#Utworzenie obiektu DGEList
dgeFull <- DGEList(rawCountTable, group=info$condition)

#usunięcie genów, które dla każdej próby mają 0 odczytów
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts,1,sum) != 0,], group=dgeFull$samples$group)

################################ ANALIZA JAKOŚCI #########################################
#utworzenie wykresu pokazującego rozmiar biblioteki dla każdej próby
par(cex.axis=0.6)
barplot(dgeFull$samples$lib.size,names=colnames(dgeFull),las=2)

# utworzenie wykresu z liczbą odczytów dla każdej próby
pseudoCount = log2(rawCountTable + 1)
pseudoCount<-has.data.frame(pseudoCount)

df = melt(pseudoCount, varnames=c('Samples', 'value'))
df = data.frame(df, Condition = c(rep('HO',12), rep('NC',15),rep('STEA',15),rep('EARLY',15)))

x <- ggplot(df, aes(x = variable, y = value, fill = Condition)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1))) + scale_fill_manual(values = c("#619CFF", "#F564E3","#5cd65c",'gold'))

x + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#wykres gęstość odczytów dla każdej próby
ggplot(df, aes(x = value, colour = variable)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.2, size = 0.25)  +
  theme(legend.position = "top") + xlab(expression(log[2](count + 1)))

#wykres MDS pokazujacy jak próbki się grupują 
plotMDS(dgeFull)

############################## ANALIZA RÓŻNICOWA ##################################

directory <- getwd()
#załaduj dane
rawCountTable <- read.csv("counts_rna.csv",header=TRUE,sep=',',row.names=2)
#usuń niepotrzebą kolumnę
rawCountTable <- rawCountTable[,-c(0:1)]
#dodaj informacje, które próby pochodzą od osób zdrowych, a które od chorych
sampleInfo <-read.csv("info.csv",header=TRUE, sep=',',row.names=1)

#Utworzenie obiektu DGEList
dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition)

#usunięcie genów, które dla każdej próby mają 0 odczytów
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ], group=dgeFull$samples$group)

# obliczenie faktora normalizacji
dgeFull <- calcNormFactors(dgeFull, method="TMM")

# normalizacja danych 
eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors
normCounts <- cpm(dgeFull)
pseudoNormCounts <- log2(normCounts + 1)

#obliczenie dyspersji
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)

# wykres rozkładu dyspersji w probach
plotBCV(dgeFull)

# test różnicowej ekspresji
dgeTest <- exactTest(dgeFull)

# filtrowanie genów o niskiej ekspresji 
filtData <- HTSFilter(dgeFull)$filteredData
dgeTestFilt <- exactTest(filtData)

# wykres FC
plotMD(dgeTestFilt)

# zapisanie wyniku do pliku 
resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
wynik <- sigDownReg <- resFilt$table[resFilt$table$FDR<0.01,]
dim(wynik)
write.csv(wynik,"DEG_RNA.csv")
