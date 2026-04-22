# Antra užduotis
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

data <- readRDS("yaskolka.rds")

# 1. Grupių palyginimas naudojant statistinį testą.
# Apskaičiuosime p-vertes DNR modifikacijoms lygindami grupes pradžioje tyrimo (T1)
# ir tyrimo pabaigoje po intervencijos (T18).

# Pirmiausia, pašalinami mėginiai nustatyti išskirtimis praeitoje užduotyje.
outlier_names = c("295_CENTRAL_T0", "295_CENTRAL_T18", "144_CENTRAL_T0", "144_CENTRAL_T18", 
                  "18_CENTRAL_T0", "18_CENTRAL_T18", "266_CENTRAL_T0", "266_CENTRAL_T18")

data_without_outliers <- data[,!colnames(data) %in% outlier_names]

# Duomenų matrica padalinama į dvi dalis pagal tai, kokiam paėmimo laikotarpiui priklauso mėginiai.
data_t0 <- data_without_outliers[,data_without_outliers$timepoint == 0]
data_t18 <- data_without_outliers[,data_without_outliers$timepoint == 18]

# Statistiniams testams naudojamas 'paired' variantas, kadangi mėginiai poromis susiję:
# tas pats donoras duoda mėginį tyrimo pradžioje ir jo pabaigoje.

# Apskaičiuojami p_values tiek t testo, tiek Vilkoksono testo metodais.
p_values_ttest <- row_t_paired(data_t0, data_t18)
p_values_wilcoxon <- row_wilcoxon_paired(data_t0, data_t18)

# Vienas iš būdų apibūdinti efekto dydį - Cohen's d.
# Jis apibrėžia efekto dydį kaip vidurkių skirtumą standartinio nuokrypio vienetais.

# Funkcija pritaikoma kiekvienai matricų atitinkamų eilučių porai ir išsaugoma [[3]] - Cohen's d reikšmė.
eff_val_cohensd <- mapply(function(i) {
  cohen.d(data_t18[i, ], data_t0[i, ], paired = TRUE)[[3]]
}, 1:nrow(data_t0))

# Kitas efekto dydis: r - koreliacija apskaičiuojama rangais panašiai kaip Wilcoxon suporuotų duomenų testu.

# Funkcija pritaikoma kiekvienai matricų atitinkamų eilučių porai.
# Išsaugoma [[1]] r reikšmė iš gauto rank_biserial objekto.
eff_val_r <- mapply(function(i) {
  rank_biserial(data_t18[i, ], data_t0[i, ], paired = TRUE)[[1]]
}, 1:nrow(data_t0))

# 2. Patikimiausių skirtumų profiliai.

# Duomenys surikiuojami pagal p-vertę didėjimo tvarka (pradedant nuo mažiausių p verčių)
cg_ordered_by_ttest <- data_without_outliers[order(p_values_ttest$pvalue),]
cg_ordered_by_wilcoxon <- data_without_outliers[order(p_values_wilcoxon$pvalue),]

# Pagal abu testus atrenkami 10 mažiausią p-vertę turinčių cg pozicijų.
top10_ttest <- cg_ordered_by_ttest[1:10,]
top10_wilcoxon <- cg_ordered_by_wilcoxon[1:10,]

top10_ttest@id
top10_wilcoxon@id
# Pagal t testą gaunamos cg pozicijos: "cg26210521", "cg10992198", "cg16872172", "cg13315471", "cg07769732",
# "cg27481720", "cg01522525", "cg20903764", "cg10509965", "cg04272309"

# Pagal Vilksono testą gaunamos cg pozicijos: "cg26210521", "cg10992198", "cg16872172", "cg13315471", "cg01522525",
# "cg17376730", "cg27481720", "cg09145071", "cg20903764", "cg02364038"

intersect(top10_ttest@id, top10_wilcoxon@id)
# Pagal abu testus sutampa 7 pozicijų pavadinimai (nebūtinai iš eilės).

# Susikuriame df objektus iš atrinktų cg pozicijų.
top10_ttest_df <- stack(top10_ttest)
top10_wilcoxon_df <- stack(top10_wilcoxon)

# 10 patikimiausių citozinų grafiškai pagal t.testą
p1 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg26210521",], 
  mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg26210521", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

p2 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg10992198",], 
  mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg10992198", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

p3 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg16872172",], 
  mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg10992198", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

p4 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg13315471",], 
  mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg13315471", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

grid.arrange(p1, p2, p3, p4)

p5 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg07769732",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg07769732", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

p6 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg27481720",], 
  mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg27481720", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

p7 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg01522525",], 
  mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg01522525", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

p8 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg20903764",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg20903764", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

grid.arrange(p5, p6, p7, p8)

p9 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg10509965",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg10509965", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

p10 <- ggplot(top10_ttest_df[top10_ttest_df$name == "cg04272309",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg04272309", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

grid.arrange(p9, p10, nrow = 1)

# 3 nesutapę citozinai Vilksono teste lyginant su t testu.
setdiff(top10_wilcoxon@id, top10_ttest@id)

p11 <- ggplot(top10_wilcoxon_df[top10_wilcoxon_df$name == "cg17376730",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg17376730", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

p12 <- ggplot(top10_wilcoxon_df[top10_wilcoxon_df$name == "cg09145071",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg09145071", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

p13 <- ggplot(top10_wilcoxon_df[top10_wilcoxon_df$name == "cg02364038",], 
             mapping=aes(x=as.factor(timepoint), y=value, color=as.factor(timepoint))) +
  labs(title = "cg02364038", x = "Laiko momentas tyrime", y = "Modifikacijos įvertis") +
  theme(plot.title = element_text(hjust=0.5, size=12), legend.position = "none") +
  geom_boxplot()

grid.arrange(p11, p12, p13, nrow = 1)

# 3. P-verčių histograma.

# P verčių histogramos pagal t testą ir Vilkoksono testą.
ggplot(p_values_ttest, aes(x=pvalue)) +
  labs(title = "T testo p reikšmių pasiskirstymas", x = "p reikšmė", y = "Kiekis") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_histogram(binwidth=0.01, fill = "darkseagreen")

# Statiškai patikimų (<=0.05) p verčių kiekis pagal t testą.
sum(p_values_ttest$pvalue <= 0.05)

ggplot(p_values_wilcoxon, aes(x=pvalue)) +
  labs(title = "Vilkoksono testo p reikšmių pasiskirstymas", x = "p reikšmė", y = "Kiekis") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_histogram(binwidth=0.01, fill = "darkseagreen")

# Statiškai patikimų (<=0.05) p verčių kiekis pagal Vilkoksono testą.
sum(p_values_wilcoxon$pvalue <= 0.05)

# 5. Manhattan grafikas

# Sudarom duomenų rinkinį manhattann plotui su reikalingais duomenimis
manhattan_df <- data.frame(
  SNP = rownames(data_without_outliers),
  CHR = as.numeric(gsub("chr", "", data_without_outliers@chr)),
  BP = data_without_outliers@pos,
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
# Grafike labiausiai išsiskiria 19 chromosoma - palyginus su kitomis chromosomomis
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

# 6. GO analizė

# Naudosime GOrilla internetinį įrankį

# Apaskaičiuojame vidurkių skirtumą (T18 - T0)
# Vidurkius panaudosime skirstant į grupes.
mean_diff <- rowMeans(data_t18) - rowMeans(data_t0)

# Pasirenkame tik patikimus skirtumus
significant_cpg <- p_values_ttest$pvalue < 0.05

# Padalijame į grupes pagal metilinimą - kur viena grupė labiau metilinta už kitą ir atvirkščiai.
hyper <- significant_cpg & mean_diff > 0   # T18 > T0
hypo  <- significant_cpg & mean_diff < 0   # T18 < T0

# Funkcija, kuri panaikins kablaitaškius ir paruoš pavadinimus pagal reikiamą formatą analizei
split_genes <- function(gene_vector) {
  genes <- unlist(strsplit(gene_vector, ";"))
  genes <- genes[genes != ""]
  unique(genes)
}

# Išrenkame visų genų pavadinimus
genes_background <- split_genes(data_without_outliers@ucsc_refgene_name)

# Padaliname genų pavadinimus pagal sukurtas grupes
genes_hyper <- split_genes(data_without_outliers@ucsc_refgene_name[hyper])
genes_hypo  <- split_genes(data_without_outliers@ucsc_refgene_name[hypo])

# Išsaugome į failus (GO analize atliksime puslapyje)
writeLines(genes_hyper, "genes_hyper.txt")
writeLines(genes_hypo, "genes_hypo.txt")
writeLines(genes_background, "genes_background.txt")
