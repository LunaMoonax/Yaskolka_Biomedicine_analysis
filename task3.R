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

data <- readRDS("yaskolka.rds")

# 1. Grupių palyginimas naudojant statistinį testą.
# Apskaičiuosime p-vertes DNR modifikacijoms lygindami grupes pradžioje tyrimo (T1)
# ir tyrimo pabaigoje po intervencijos (T18).

# Pirmiausia, pašalinami mėginiai nustatyti išskirtimis praeitoje užduotyje.
outlier_names = c("295_CENTRAL_T0", "295_CENTRAL_T18", "144_CENTRAL_T0", "144_CENTRAL_T18", 
                  "18_CENTRAL_T0", "18_CENTRAL_T18", "266_CENTRAL_T0", "266_CENTRAL_T18")

data <- data[,!colnames(data) %in% outlier_names]

# Duomenų matrica atskiriama į dvi dalis pagal tai, kokiam paėmimo laikotarpiui priklauso mėginiai.
# Statistiniams testams bus naudojamas 'paired' variantas, kadangi mėginiai poromis susiję:
# tas pats donoras duoda mėginį tyrimo pradžioje ir jo pabaigoje.
# Todėl reiktų patikrinti ar mėginiai bus teisingai suporuojami 'paired' testuose.

# Ištraukiami tiriamųjų numeriai iš stulpelių pavadinimų.
colnames(data[,data$timepoint == 0])
colnames(data[,data$timepoint == 18])

t0 <- sub("_.*", "", colnames(data[,data$timepoint == 0]))
t18 <- sub("_.*", "", colnames(data[,data$timepoint == 18]))

# Patikrinimas ar sutampa abiejų grupių tiriamųjų numeriai bei jų eilės tvarka.
identical(t0, t18)

# Apskaičiuojami p_values tiek t testo, tiek Vilkoksono testo metodais.
p_values_ttest <- row_t_paired(data[,data$timepoint == 18], data[,data$timepoint == 0])
p_values_wilcoxon <- row_wilcoxon_paired(data[,data$timepoint == 18], data[,data$timepoint == 0])

# P reikšmių korekcija "fdr" metodu.
p_values_ttest$pvalue <- p.adjust(p_values_ttest$pvalue, method = "fdr")
p_values_wilcoxon$pvalue <- p.adjust(p_values_wilcoxon$pvalue, method = "fdr")

# Efekto dydis prieinamas p_values_ttest$mean.diff

# 2. Patikimiausių skirtumų profiliai.

# Duomenys surikiuojami pagal p-vertę didėjimo tvarka (pradedant nuo mažiausių p verčių)
p_values_ttest_ordered <- p_values_ttest[order(p_values_ttest$pvalue),]
p_values_wilcoxon_ordered <- p_values_wilcoxon[order(p_values_wilcoxon$pvalue),]

# Pagal abu testus atrenkami 10 mažiausią p-vertę turinčių cg pozicijų.
top10_ttest <- rownames(p_values_ttest_ordered[1:10,])
top10_wilcoxon <- rownames(p_values_wilcoxon_ordered[1:10,])

top10_ttest
top10_wilcoxon

# Pagal t testą gaunamos cg pozicijos: "cg26210521", "cg10992198", "cg16872172", "cg13315471", "cg07769732",
# "cg27481720", "cg01522525", "cg10509965", "cg20903764", "cg04272309"

# Pagal Vilkoksono testą gaunamos cg pozicijos: "cg26210521", "cg16872172", "cg10992198", "cg13315471", 
# "cg01522525", "cg17376730", "cg27481720", "cg09145071", "cg20903764", "cg02364038"

intersect(top10_ttest, top10_wilcoxon)
# Pagal abu testus sutampa 7 pozicijų pavadinimai (nebūtinai iš eilės).

# Susikuriame df objektus iš atrinktų cg pozicijų.
top10_ttest_df <- stack(data[rownames(data) %in% top10_ttest,])
top10_wilcoxon_df <- stack(data[rownames(data) %in% top10_wilcoxon,])

# 10 patikimiausių cg pozicijų grafiškai pagal t.testą
p1 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg26210521",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg26210521", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

p2 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg10992198",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg10992198", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

grid.arrange(p1, p2, nrow = 1)

p3 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg16872172",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg16872172", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

p4 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg13315471",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg13315471", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

grid.arrange(p3, p4,nrow = 1)

p5 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg07769732",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg07769732", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

p6 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg27481720",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg27481720", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

grid.arrange(p5, p6,nrow = 1)

p7 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg01522525",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg01522525", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

p8 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg10509965",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg10509965", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

grid.arrange(p7, p8, nrow = 1)

p9 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg20903764",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg20903764", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

p10 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg04272309",], 
              mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg04272309", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

grid.arrange(p9, p10, nrow = 1)

# 3 nesutapusios CpG pozicijos Vilksono teste lyginant su t testu.
setdiff(top10_wilcoxon, top10_ttest)

p11 <- ggplot(top10_wilcoxon_df[top10_wilcoxon_df$name == "cg17376730",], 
              mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg17376730", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

p12 <- ggplot(top10_wilcoxon_df[top10_wilcoxon_df$name == "cg09145071",], 
              mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg09145071", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) + 
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

p13 <- ggplot(top10_wilcoxon_df[top10_wilcoxon_df$name == "cg02364038",], 
              mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint), fill = as.factor(timepoint))) +
  labs(title = "cg02364038", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(width = 0.3, alpha = 0.8) +
  scale_color_manual(values = c("darkolivegreen", "slateblue")) +
  scale_fill_manual(values = c("darkolivegreen", "slateblue"))

grid.arrange(p11, p12, p13, nrow = 1)

# 3. P-verčių histograma.

# P verčių histogramos pagal t testą ir Vilkoksono testą.
ggplot(p_values_ttest, aes(x=pvalue)) +
  labs(title = "T testo p verčių pasiskirstymas", x = "p reikšmė", y = "Kiekis") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_histogram(binwidth=0.01, fill = "darkolivegreen")

# Statiškai patikimų (<0.05) p verčių kiekis pagal t testą.
sum(p_values_ttest$pvalue < 0.05)

# Kokią dalį visų CpG pozicijų sudaro statistiškai patikimi CpG (procentais).
sum(p_values_ttest$pvalue < 0.05) / length(rownames(data)) * 100

ggplot(p_values_wilcoxon, aes(x=pvalue)) +
  labs(title = "Vilkoksono testo p verčių pasiskirstymas", x = "p reikšmė", y = "Kiekis") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_histogram(binwidth=0.01, fill = "darkolivegreen")

# Statiškai patikimų (<0.05) p verčių kiekis pagal Vilkoksono testą.
sum(p_values_wilcoxon$pvalue < 0.05)

# Kokią dalį visų CpG pozicijų sudaro statistiškai patikimi CpG.
sum(p_values_wilcoxon$pvalue < 0.05) / length(rownames(data)) * 100

# 4. Volcano grafikas

# Sukuriame duomenų lentelę, skirtą volcano grafikui pavaizduoti
volcano_df <- data.frame(
  mean_diff = p_values_ttest$mean.diff,
  logp = -log10(p_values_ttest$pvalue)
)

# Tam, kad grafikas būtų aiškesnis ir duotų daugiau informacijos, duomenis paskirstome kategorijomis
# Visus duomenis kol kas priskiriame nereikšmingus.
# Reikšmingi tampa tie, kurių p_value dydis yra mažesnis negu 0.05 ir atitinkamai efekto dydis
# priskiriamas hipermetilintai, arba hipometilintam.
volcano_df$category <- "Nereikšmingas"
volcano_df$category[p_values_ttest$pvalue < 0.05 & volcano_df$mean_diff > 0.01] <- "Hipermetilinta"
volcano_df$category[p_values_ttest$pvalue < 0.05 & volcano_df$mean_diff < -0.01] <- "Hipometilinta"

# Pavaizduojame volcano grafiką
ggplot(volcano_df, aes(x = mean_diff, y = logp, color = category)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("Hipermetilinta" = "red", "Hipometilinta" = "blue", "Nereikšmingas" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.01, 0.01), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano grafikas",
    x = "Vidurkių skirtumas (T18 - T0)",
    y = expression(-log[10](p[adj])),
    color = "Kategorija"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

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