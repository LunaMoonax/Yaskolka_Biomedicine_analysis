# Antra užduotis
# Viktorija Ramonaitė, Skaistė Bartkutė

setwd("C:/Users/Viktorija Ramonaite/Desktop/UNIVERAS/3 KURSAS/BIOMEDICINOS DUOMENU ANALIZE/1 uzduotis/Yaskolka_Biomedicine_analysis")
#setwd("C:/Users/skais/Desktop/Universitetas/Šeštas semestras/Biomedicina/task1/Yaskolka_Biomedicine_analysis")

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
p_values_ttest <- row_t_paired(data[,data$timepoint == 0], data[,data$timepoint == 18])
p_values_wilcoxon <- row_wilcoxon_paired(data[,data$timepoint == 0], data[,data$timepoint == 18])

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

# 10 patikimiausių citozinų grafiškai pagal t.testą
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

# Kokią dalį visų CpG pozicijų sudaro statistiškai patikimi CpG.
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

# 6. GO analizė

# Naudosime GOrilla internetinį įrankį

# Apaskaičiuojame vidurkių skirtumą (T18 - T0)
# Vidurkius panaudosime skirstant į grupes.
mean_diff <- rowMeans(data[,data$timepoint == 18]) - rowMeans(data[,data$timepoint == 0])

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

# Išrenkame visų unikalių genų pavadinimus
genes_background <- split_genes(data@ucsc_refgene_name)

# Padaliname unikalius genų pavadinimus pagal sukurtas grupes
genes_hyper <- split_genes(data@ucsc_refgene_name[hyper])
genes_hypo  <- split_genes(data@ucsc_refgene_name[hypo])

# Išsaugome į failus (GO analize atliksime puslapyje)
writeLines(genes_hyper, "genes_hyper.txt")
writeLines(genes_hypo, "genes_hypo.txt")
writeLines(genes_background, "genes_background.txt")

# Sipriausių efekto dydžių apžvalginė analizė
# Vienas iš būdų apibūdinti efekto dydį - Cohen's d.
# Jis apibrėžia efekto dydį kaip vidurkių skirtumą standartinio nuokrypio vienetais.

# Funkcija pritaikoma kiekvienai matricų atitinkamų eilučių porai ir išsaugoma [[3]] - Cohen's d reikšmė.
p_values_ttest$eff_val_cohensd <- mapply(function(i) {
  cohen.d(data[i, data$timepoint == 18], data[i, data$timepoint == 0], paired = TRUE)[[3]]
}, 1:nrow(data))

# Kitas efekto dydis: r - koreliacija apskaičiuojama rangais panašiai kaip Wilcoxon suporuotų duomenų testu.

# Funkcija pritaikoma kiekvienai matricų atitinkamų eilučių porai.
# Išsaugoma [[1]] r reikšmė iš gauto rank_biserial objekto.
p_values_wilcoxon$eff_val_r <- mapply(function(i) {
  rank_biserial(data[i, data$timepoint == 18], data[i, data$timepoint == 0], paired = TRUE)[[1]]
}, 1:nrow(data))

# Atrenkamos tik CpG pozicijos patikimos pagal p vertę. 
significant_ttest <- p_values_ttest[p_values_ttest$pvalue < 0.05,]
# Peržiūrimos didžiausias vidurkių skirtumas ir pagal tai nustatomas slenkstis.
max(abs(significant_ttest$mean.diff))
# Didžiausias skirtumas 0.06465848, todėl bus naudojama apie pusė šio skirtumo (0.03).
eff_val_diff <- significant_ttest[abs(significant_ttest$mean.diff) > 0.03,]

# Cohen's d ir r slenktis bus 0.5 (vidutinis efektas).
eff_val_cohensd_strong <- significant_ttest[abs(significant_ttest$eff_val_cohensd) > 0.5,]
significant_wilcoxon <- p_values_wilcoxon[p_values_wilcoxon$pvalue < 0.05,]
eff_val_r_strong <- significant_wilcoxon[abs(significant_wilcoxon$eff_val_r) > 0.5,]

# Į abiejus sąrašus pagal Cohen's d ir r efekto dydžius patenkančios CpG pozicijos.
intersect(rownames(eff_val_cohensd_strong), rownames(eff_val_r_strong))
length(intersect(rownames(eff_val_cohensd_strong), rownames(eff_val_r_strong)))

# Į abiejus sąrašus pagal Cohen's d ir vidurkių skirtumo efekto dydžius patenkančios CpG pozicijos.
intersect(rownames(eff_val_cohensd_strong), rownames(eff_val_diff))
length(intersect(rownames(eff_val_cohensd_strong), rownames(eff_val_diff)))

# Į abiejus sąrašus pagal r ir vidurkių skirtumo efekto dydžius patenkančios CpG pozicijos.
intersect(rownames(eff_val_r_strong), rownames(eff_val_diff))
length(intersect(rownames(eff_val_r_strong), rownames(eff_val_diff)))

# Kadangi pagal vidurkių skirtumo sąrašą su kitais sąrašais sutampa tik 7 CpG pozicijos,
# vidurkių skirtumo slenkstis sumažinamas.
eff_val_diff <- significant_ttest[abs(significant_ttest$mean.diff) > 0.01,]

# Sutapimų skaičius padidėjo iki 124.
intersect(rownames(eff_val_cohensd_strong), rownames(eff_val_diff))
length(intersect(rownames(eff_val_cohensd_strong), rownames(eff_val_diff)))

# Sutapimų skaičius padidėjo iki 361.
intersect(rownames(eff_val_r_strong), rownames(eff_val_diff))
length(intersect(rownames(eff_val_r_strong), rownames(eff_val_diff)))

# Toliau bus pavaizduojama kiek CpG pozicijų pagal vidutinius ir stiprius efekto dydžius patenka į skirtingus regionus.
# Susidarome dataframe grafikams.
cohensd_df <- rowanns(data[rownames(data) %in% rownames(eff_val_cohensd_strong),])
r_df <- rowanns(data[rownames(data) %in% rownames(eff_val_r_strong),])
diff_df <- rowanns(data[rownames(data) %in% rownames(eff_val_diff),])

ggplot(cohensd_df, aes(x=relation_to_island)) +
  labs(title = "CpG pasiskirstymas tarp regionų (atrinktos pagal vidutinį Cohen's d efekto dydį)", x = "Regionas", y = "CpG pozicijų kiekis regione") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_bar(fill = "slateblue")

ggplot(r_df, aes(x=relation_to_island)) +
  labs(title = "CpG pasiskirstymas tarp regionų (atrinktos pagal vidutinį r (biserial rank) efekto dydį)", x = "Regionas", y = "CpG pozicijų kiekis regione") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_bar(fill = "slateblue")

ggplot(diff_df, aes(x=relation_to_island)) +
  labs(title = "CpG pasiskirstymas tarp regionų (atrinktos pagal vidurkių skirtumo (>0.01) efekto dydį)", x = "Regionas", y = "CpG pozicijų kiekis regione") +
  theme(plot.title = element_text(hjust=0.5, size=12)) +
  geom_bar(fill = "slateblue")

# Apskaičiavimas kiek CpG pozicijų patenka į kiekvieną sąrašą.

length(rownames(eff_val_cohensd_strong))
length(rownames(eff_val_r_strong))
length(rownames(eff_val_diff))

# Rikiavimas pagal efekto dydžių stiprumą (nuo didžiausio efekto dydžio).
eff_val_cohensd_strong <- eff_val_cohensd_strong[order(abs(eff_val_cohensd_strong$eff_val_cohensd), decreasing = TRUE),]
rownames(eff_val_cohensd_strong[1:10,])

eff_val_r_strong <- eff_val_r_strong[order(abs(eff_val_r_strong$eff_val_r), decreasing = TRUE),]
rownames(eff_val_r_strong[1:10,])

eff_val_diff <- eff_val_diff[order(abs(eff_val_diff$mean.diff), decreasing = TRUE),]
rownames(eff_val_diff[1:10,])

# Sutapimai tarp top 10 didžiausių efekto dydžių (pagal skirtingus efekto dydžių tipus).
intersect(rownames(eff_val_cohensd_strong[1:10,]), rownames(eff_val_r_strong[1:10,]))
# cg01522525

intersect(rownames(eff_val_cohensd_strong[1:10,]), rownames(eff_val_diff[1:10,]))
# cg05001044 ir cg22519184

intersect(rownames(eff_val_r_strong[1:10,]), rownames(eff_val_diff[1:10,]))
# nėra bendrų CpG pozicijų.

# Heatmaps pagal didžiausius efekto dydžius.
heatmap(data[rownames(data) %in% rownames(eff_val_cohensd_strong[1:10,]),])
heatmap(data[rownames(data) %in% rownames(eff_val_r_strong[1:10,]),])
heatmap(data[rownames(data) %in% rownames(eff_val_diff[1:10,]),])