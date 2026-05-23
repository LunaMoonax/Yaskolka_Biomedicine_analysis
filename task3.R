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

# Top CpG scatter plotas — stipriausias amžiaus ir metilinimo ryšys
top1 <- which.min(pvalues_adj)
ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top1,])),
       aes(x = age, y = methylation)) +
  geom_point(alpha = 0.5, color = "steelblue") +
  geom_smooth(method = "lm", color = "red") +
  labs(
    title = paste("Metilinimo priklausomybė nuo amžiaus:", rownames(data)[top1]),
    x = "Amžius", y = "Metilinimas"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# 2. Patikimiausių skirtumų profiliai.

# Duomenys surikiuojami pagal p-vertę didėjimo tvarka (pradedant nuo mažiausių p verčių).
pvalues_adj <- sort(pvalues_adj)

# Atrenkame 10 mažiausią p-vertę turinčių citozinų.
top10_CPGs <- names(pvalues_adj[1:10])
top10_CPGs

# Tai citozinai: cg16867657, cg17268658, cg21572722, cg22454769, cg06639320,
# cg24724428, cg08097417, cg07553761, cg26947034, cg17403084.

# Top 10 CpG scatter plot'ai — stipriausi amžiaus ir metilinimo ryšiai pagal pvalue.

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[1],])),
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

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[2],])),
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

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[3],])),
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

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[4],])),
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

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[5],])),
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

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[6],])),
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

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[7],])),
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

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[8],])),
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

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[9],])),
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

ggplot(data.frame(age = colanns(data)$age, methylation = as.numeric(data[top10_CPGs[10],])),
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
  P = p_values_ttest$pvalue
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
# Pasirenkame išsiskiriančia 19 chromosoma - palyginus su kitomis chromosomomis
# ji turi daugiausiai taškų, kurie yra aukštai (ne tik pavienį tašką).

# Kadangi norime atvaizduoti tik vieną chromosomą ir suprantamai, čia pasinaudosime ggplot
# Paruošiame duomenis grafikui, logoritmuojame ir konvertuojame pozicijas į megabazes (Mb)
chr19 <- manhattan_df[manhattan_df$CHR == 19, ]
chr19$logp <- -log10(chr19$P)
chr19$pos_mb <- chr19$BP / 1e6

# Braižome manhattan plotą 19 chromosomai
ggplot(chr19, aes(x = pos_mb, y = logp)) +
  geom_point(size = 1, alpha = 0.7, color = "steelblue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(0, 60, by = 5)) +
  labs(
    title = "Manhattan grafikas (19 chromosoma)",
    x = "Genominė pozicija (Mb)",
    y = expression(-log[10](p))
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))