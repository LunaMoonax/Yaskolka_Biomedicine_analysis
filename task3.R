# Trečia užduotis
# Viktorija Ramonaitė, Skaistė Bartkutė

#setwd("C:/Users/Viktorija Ramonaite/Desktop/UNIVERAS/3 KURSAS/BIOMEDICINOS DUOMENU ANALIZE/1 uzduotis/Yaskolka_Biomedicine_analysis")
setwd("C:/Users/skais/Desktop/Universitetas/Šeštas semestras/Biomedicina/task1/Yaskolka_Biomedicine_analysis")

library(annmatrix)
library(ggplot2)
library(matrixTests)
library(effsize)
library(effectsize)
library(gridExtra)
library(qqman)
library(pheatmap)
library(limma)

data <- readRDS("yaskolka.rds")

# Pirmiausia, pašalinami mėginiai nustatyti išskirtimis.
outlier_names = c("295_CENTRAL_T0", "295_CENTRAL_T18", "144_CENTRAL_T0", "144_CENTRAL_T18", 
                  "18_CENTRAL_T0", "18_CENTRAL_T18", "266_CENTRAL_T0", "266_CENTRAL_T18")

data <- data[,!colnames(data) %in% outlier_names]

# 1. Raskite senėjimo tendencijas

# Pasitikriname galimus kokybinius kintamuosius, kuriuos būtų galima įtraukti į modelį
# Pagrinde ieškome tokių, kuriuose būtų variabilumas
head(colanns(data))
table(data$celltype)
table(data$diagnosis)
table(data$diet)
table(data$sex)
table(data$stimulus)
table(data$timepoint)
table(data$timepoint, data$stimulus)

# Naudojame limma paketą, kuris yra optimizuotas tiesinės regresijos modeliams taikyti
# didelėms genomikos matricoms. Viduje naudoja tą patį tiesinį regresijos modelį kaip lm(),
# tačiau yra žymiai greitesnis, nes apdoroja visas CpG pozicijas vienu metu,
# o ne po vieną cikle.

# Sukuriame design matricą su visais kintamaisiais:
# age - mus dominantis kintamasis (senėjimo efektas)
# sex, diet, stimulus, timepoint - kontroliniai kintamieji,
# kuriuos kontroliuojame.
design <- model.matrix(~ age + sex + diet + stimulus + timepoint, data = colanns(data))

# Pritaikome tiesinį regresijos modelį visoms CpG pozicijoms vienu metu.
fit <- lmFit(as.matrix(data), design)

# Ištraukiame amžiaus beta koeficientus
betas_age <- fit$coefficients[,"age"]
# P-vertes skaičiuojame iš t-statistikos ir laisvės laipsnių
t_stat <- fit$coefficients[,"age"] / fit$stdev.unscaled[,"age"] / fit$sigma
pvalues <- 2 * pt(-abs(t_stat), df = fit$df.residual)

# P-verčių korekcija FDR metodu
pvalues_adj <- p.adjust(pvalues, method = "fdr")

# Kiek citozinų rodo su amžiumi susijusias tendencijas?
significant <- pvalues_adj < 0.05
sum(significant)
sum(significant) / nrow(data) * 100

# Kokios modifikavimo kryptys? (hipermetilinimas vs hipometilinimas)
sum(betas_age[significant] > 0)   # hipermetilinta
sum(betas_age[significant] < 0)   # hipometilinta

# Kokiuose genominiuose kontekstuose šie citozinai randami?
table(rowanns(data[significant,])$relation_to_island)
table(unlist(strsplit(rowanns(data[significant,])$ucsc_refgene_group, ";")))

# Volcano grafikas
# Susidarome dataframe reikalingą grafikui atvaizduoti
volcano_df <- data.frame(
  beta = betas_age,
  logp = -log10(pvalues_adj)
)
# Suteikiame pavadinimus kategorijomis, kad grafike matytusi
volcano_df$category <- "Nereikšmingas"
volcano_df$category[pvalues_adj < 0.05 & betas_age > 0] <- "Hipermetilinta"
volcano_df$category[pvalues_adj < 0.05 & betas_age < 0] <- "Hipometilinta"

# Piešiame plotą
ggplot(volcano_df, aes(x = beta, y = logp, color = category)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Hipermetilinta" = "red", "Hipometilinta" = "blue", "Nereikšmingas" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Senėjimo tendencijos (Volcano grafikas)",
    x = "Beta koeficientas (amžiaus efektas per metus)",
    y = expression(-log[10](p[adj])),
    color = "Kategorija"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Genominio konteksto stulpelinė diagrama su kryptimis
# Susidarome dataframe reikalingą grafikui atvaizduoti
sig_context <- data.frame(
  region = rowanns(data[significant,])$relation_to_island,
  direction = ifelse(betas_age[significant] > 0, "Hipermetilinta", "Hipometilinta")
)
# Piešiame plotą
ggplot(sig_context, aes(x = region, fill = direction)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("Hipermetilinta" = "red", "Hipometilinta" = "blue")) +
  labs(
    title = "Su amžiumi susijusių CpG pasiskirstymas pagal genominį kontekstą",
    x = "Regionas", y = "CpG kiekis", fill = "Kryptis"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# 2. Patikimiausių skirtumų profiliai.

# Duomenys surikiuojami pagal p-vertę didėjimo tvarka (pradedant nuo mažiausių p verčių).
sorted_pvalues_adj <- sort(pvalues_adj)

# Atrenkame 10 mažiausią p-vertę turinčių citozinų.
top10_CPGs <- names(sorted_pvalues_adj[1:10])
top10_CPGs

# Tai citozinai: cg16867657, cg17268658, cg21572722, cg22454769, cg06639320,
# cg24724428, cg08097417, cg07553761, cg26947034, cg17403084.

# Top 10 CpG scatter plot'ai — stipriausi amžiaus ir metilinimo ryšiai pagal pvalue.

p1 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[1],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[1]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[2],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[2]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p1, p2)

p3 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[3],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[3]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[4],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[4]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p3, p4)

p5 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[5],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[5]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p6 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[6],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[6]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p5, p6)

p7 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[7],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[7]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p8 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[8],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[8]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p7, p8)

p9 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[9],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[9]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

p10 <- ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[10],])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, aes(color=factor(colanns(data)$timepoint))) +
  scale_color_manual(values = c("0" = "blue", "18" = "red")) +
  geom_smooth(method = "lm", color = "green") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", top10_CPGs[10]),
    x = "Amžius", y = "Metilinimas", color = "Laiko momentas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(p9, p10)

# 3. P-verčių histograma.

# P verčių histograma.
ggplot(data.frame(pvalue = pvalues_adj), aes(x=pvalue)) +
  labs(title = "P verčių pasiskirstymas", x = "P vertė", y = "Kiekis") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_histogram(binwidth=0.01, fill = "darkolivegreen")

# Statistiškai patikimų citozinų dalis apskaičiuota anksčiau.

# 5. Manhattan grafikas

# Sudarom duomenų rinkinį manhattann plotui su reikalingais duomenimis
manhattan_df <- data.frame(
  SNP = rownames(data),
  CHR = as.numeric(gsub("chr", "", data@chr)),
  BP = data@pos,
  P = pvalues_adj
)

# Pašaliname tuščias reikšmes, kad plotas būtų atvaizduotas tiksliai
manhattan_df <- manhattan_df[!is.na(manhattan_df$CHR) & !is.na(manhattan_df$P) & !is.na(manhattan_df$BP), ]

# Braižome manhattan plotą
manhattan(manhattan_df, main = "Manhattan grafikas (visos chromosomos)",
          suggestiveline = FALSE, genomewideline = -log10(0.05),
          col = c("steelblue", "darkorange"))
# Dauguma taškų yra patikimi (virš raudonos linijos, kuri žymi p=0.05). Kuo aukštesnis taškas
# tuo patikimesnis skirtumas (p rekšmė mažiausia). 
# Grafike dauguma chromosomų atrodo panašios, ir itin didelių skirtumų sunku rasti.
# Pasirenkame išsiskiriančia 17 chromosoma - palyginus su kitomis chromosomomis
# ji turi daugiausiai taškų, kurie yra aukštai (ne tik pavienį tašką).

# Kadangi norime atvaizduoti tik vieną chromosomą ir suprantamai, čia pasinaudosime ggplot
# Paruošiame duomenis grafikui, logoritmuojame ir konvertuojame pozicijas į megabazes (Mb)
chr17 <- manhattan_df[manhattan_df$CHR == 17, ]
chr17$logp <- -log10(chr17$P)
chr17$pos_mb <- chr17$BP / 1e6

# Braižome manhattan plotą 17 chromosomai
ggplot(chr17, aes(x = pos_mb, y = logp)) +
  geom_point(size = 1, alpha = 0.7, color = "steelblue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(0, 80, by = 5)) +
  labs(
    title = "Manhattan grafikas (17 chromosoma)",
    x = "Genominė pozicija (Mb)",
    y = expression(-log[10](p))
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# 5. Horvath epigenetinis laikrodis
DNAm <- read.csv("C:/Users/skais/Desktop/Universitetas/Šeštas semestras/Biomedicina/task1/13059_2013_3156_MOESM3_ESM.csv", skip = 2)

length(DNAm$CpGmarker)
# Kartu su Intercept yra 353 CpG markeriai epigenetiniam amžiui nustatyti.

age_data <- data[rownames(data) %in% DNAm$CpGmarker,]
dim(age_data)
# Duomenyse atrasti 334 CpG markeriai iš 353 markerių.

# Išsisaugome intercept reikšmę.
intercept <- DNAm$CoefficientTraining[1]

# Pasiliekame tik tuos epigenetinio laikrodžio CpGs, kurie buvo rasti duomenyse.
# Perrikiuojame age data pagal CpGmarker, kad modifikacijos įverčiai būtų padauginti iš teisingų koeficientų.
DNAm_CPGs <- DNAm[DNAm$CpGmarker %in% rownames(age_data),]
age_data <- age_data[DNAm_CPGs$CpGmarker,]
rownames(age_data) == DNAm_CPGs$CpGmarker
length(DNAm_CPGs$CpGmarker)

# Padauginame koeficientus atitinkamo CpG iš jo modifikacijos įverčio ir sudedame į galutinę mėginio sumą.
epigenetic_age <- colSums(age_data * DNAm_CPGs$CoefficientTraining) + intercept

# Norime pridėti efektus ir tų CpG, kurių nebuvo duomenyse.
# Todėl padauginsime nepanaudotų CpG medianas iš jų koeficiento.
DNAm_CPGs_not_found <- DNAm[!DNAm$CpGmarker %in% rownames(age_data),]

# Pašaliname intercept, kad nesikartotų.
DNAm_CPGs_not_found <- DNAm_CPGs_not_found[-1, ]

# Dauginame CpG medianas iš atitinkamo jų koeficiento ir gauname galutinę sumą jų efekto.
coeff <- sum(DNAm_CPGs_not_found$CoefficientTraining * DNAm_CPGs_not_found$medianByCpG)
epigenetic_age_w_coeff <- epigenetic_age + coeff

# Horvath apibrėžta gautų įverčių transformavimo funkcija (supplementary failuose).
anti.trafo <- function(x,adult.age=20) { 
   if(x<0) {
     (1+adult.age)*exp(x)-1
   } else {
     (1+adult.age)*x+adult.age 
   }
}

# Kiekvienam tiriamajam pritaikome amžiaus transformaciją.
epigenetic_age_w_coeff <- sapply(epigenetic_age_w_coeff, anti.trafo)
epigenetic_age <- sapply(epigenetic_age, anti.trafo)

# Sudarome data frames grafikams: realus amžius ir epigentinis amžius su papildomais CpG
age_df_1 <- data.frame(
  ages = c(age_data$age, epigenetic_age_w_coeff),
  groups = c(rep("real", length(age_data$age)),
            rep("epi", length(epigenetic_age_w_coeff))),
  cols = colnames(age_data),
  timepoint = age_data$timepoint,
  donor = age_data$donor
)

# Realus amžius ir epigentinis amžius be papildomų CpG
age_df_2 <- data.frame(
  ages = c(age_data$age, epigenetic_age),
  groups = c(rep("real", length(age_data$age)),
             rep("epi", length(epigenetic_age))),
  cols = colnames(age_data),
  timepoint = age_data$timepoint,
  donor = age_data$donor
)

# Dėl tvarkingumo pašalinami mėginių vardai.
ggplot(age_df_1, aes(x=cols,y=ages,color=groups,group=cols)) +
  labs(title = "Tikras amžius vs. epigenetinis amžius (su papildomais CpG)", x = "Tiriamieji", y = "Amžius") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_point(alpha = 0.5) +
  geom_line(color = "grey70") +
  scale_color_manual(values = c("real" = "red", "epi" = "blue"),
                     labels = c("real" = "realus amžius", "epi" = "epigenetinis amžius"),
                     name = "Amžiaus kategorija") +
  theme(axis.text.x = element_blank())

ggplot(age_df_2, aes(x=cols,y=ages,color=groups,group=cols)) +
  labs(title = "Tikras amžius vs. epigenetinis amžius (be papildomų CpG)", x = "Tiriamieji", y = "Amžius") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_point(alpha = 0.5) +
  geom_line(color = "grey70") +
  scale_color_manual(values = c("real" = "red", "epi" = "blue"),
                     labels = c("real" = "realus amžius", "epi" = "epigenetinis amžius"),
                     name = "Amžiaus kategorija") +
  theme(axis.text.x = element_blank())

# Vidutiniai skirtumai tarp epigenetinio amžiaus ir realaus amžiaus.
mean(epigenetic_age_w_coeff - age_data$age)
mean(epigenetic_age - age_data$age)

# Standartinės paklaidos tarp epigenetinio amžiaus ir realaus amžiaus skirtumų.
sd(epigenetic_age_w_coeff - age_data$age) /sqrt(length(epigenetic_age_w_coeff))
sd(epigenetic_age - age_data$age) /sqrt(length(epigenetic_age))

# Gauname vienodas standartines paklaidas tiek naudojant papildomus CpG, tiek nenaudojant.
# Tačiau skaičiavimas nenaudojant papildomų CpG atrodo tikslesnis, 
# kadangi vidutinis amžiaus skirtumas tarp chronologinio ir epigenetinio mažesnis.

# Apžvalginė analizė

# Žinome, kad chronologinis amžiaus skirtumas prieš ir po analizės yra 1,5 metų,
# tačiau norime patyrinėti kaip pasikeitė epigenetinis amžius prieš ir po tyrimo.

ggplot(age_df_2[age_df_2$groups == "epi",], aes(x=donor,y=ages, color=factor(timepoint), group=donor)) +
  labs(title = "Epigenetinio amžiaus amplitudė", x = "Tiriamieji", y = "Amžius", color = "Mėginio paėmimo momentas") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_point(alpha = 0.5) +
  geom_line(color = "grey70") +
  theme(axis.text.x = element_blank())

# Iš grafiko matyti, kad kai kurių tiriamųjų epigenetinis amžius tyrimo eigoje sumažėjo.
# Atsirenkame epigenetinius amžius, apskaičiuojame jų skirtumą tyrimo pradžioje ir pabaigoje,
# šiuos skirtumus priskiriame tiriamiesiems.
epi_age_df <- age_df_2[age_df_2$groups == "epi",]
epi_age_diff <- epi_age_df[epi_age_df$timepoint == 18,]$ages - epi_age_df[epi_age_df$timepoint == 0,]$ages
unique_donor <- epi_age_df[!duplicated(epi_age_df$donor),]
unique_donor$age_diff <- epi_age_diff

# Atsirenkame tuos tiriamuosius, kurių epigenetinis amžius sumažėjo tyrimo eigoje.
donors_neg_diff <- unique_donor[unique_donor$age_diff < 0,]$donor
length(donors_neg_diff)

# Atvaizduosime epigenetiškai atjaunėjusių tiriamųjų charakteristikas.
# Kad dėl poruotų mėginių nesidubliuotų stulpelių rezultatai diet, stimulus, pasiliekame tik T0 mėginius
# (diet ir stimulus tam pačiam tiriamajam  T0 ir T18 nurodomi vienodi)
age_factors <- colanns(age_data)[colanns(age_data)$timepoint == 0 & colanns(age_data)$donor %in% donors_neg_diff,]

g1 <- ggplot(age_factors, aes(x=diet)) +
  labs(title = "Epigenetiškai atjaunėjusių tiriamųjų dieta", x = "Dietos tipas", y = "Tiriamųjų kiekis") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_bar(fill = "slateblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

g2 <- ggplot(age_factors, aes(x=stimulus)) +
  labs(title = "Epigenetiškai atjaunėjusių tiriamųjų fizinis aktyvumas", x = "Fizinis aktyvumo (ne)buvimas", y = "Tiriamųjų kiekis") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_bar(fill = "slateblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

g3 <- ggplot(age_factors, aes(x=sex)) +
  labs(title = "Epigenetiškai atjaunėjusių tiriamųjų lytis", x = "Lytis", y = "Tiriamųjų kiekis") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_bar(fill = "slateblue") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

grid.arrange(g1, g2, g3)