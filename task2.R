# Antra užduotis
# Viktorija Ramonaitė, Skaistė Bartkutė

#setwd("C:/Users/Viktorija Ramonaite/Desktop/UNIVERAS/3 KURSAS/BIOMEDICINOS DUOMENU ANALIZE/1 uzduotis/Yaskolka_Biomedicine_analysis")
setwd("C:/Users/skais/Desktop/Universitetas/Šeštas semestras/Biomedicina/task1/Yaskolka_Biomedicine_analysis")

library(annmatrix)
library(ggplot2)
library(matrixTests)
library(effectsize)
library(gridExtra)

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