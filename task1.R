# Pirma užduotis
setwd("C:/Users/Viktorija Ramonaite/Desktop/UNIVERAS/3 KURSAS/BIOMEDICINOS DUOMENU ANALIZE/1 uzduotis/Yaskolka_Biomedicine_analysis")
setwd("C:/Users/skais/Desktop/Universitetas/Šeštas semestras/Biomedicina/task1/Yaskolka_Biomedicine_analysis")

library(annmatrix)
library(ggplot2)
library(dplyr)
library(patchwork)

data <- readRDS("yaskolka.rds")

# Kokybės kontrolė #1

# Annmatrix reikšmių peržiūra
rowanns(data)
colanns(data)

# Peržiūrime CpG salų reikšmes
data@relation_to_island

# Suskaičiuojame metilinimo lygių (beta) vidutines reikšmes eilutėms
mean_beta <- rowMeans(data, na.rm = TRUE)

# Sukuriame dataframe, kuri bus naudojama grafike
df <- data.frame(
  mean_beta = mean_beta,
  region = data@relation_to_island
)

# Suskirstome relation_to_island reikšmes 4 grupes, kurios mums yra reikalingos kokybės kontrolei patvirtinti
df$region_grouped <- case_when(
  df$region == "Island" ~ "Island",
  df$region == "OpenSea" ~ "OpenSea",
  df$region %in% c("N_Shelf", "S_Shelf") ~ "Shelf",
  df$region %in% c("N_Shore", "S_Shore") ~ "Shore"
)
# Sugrupuojame lygiais, kad grafike būtų atvaizduota mums reikalinga tvarka
df$region_grouped <- factor(df$region_grouped, levels = c("Island", "Shore", "Shelf", "OpenSea"))

# Piešiame grafiką
ggplot(df, aes(x = mean_beta, color = region_grouped)) +
  geom_density(linewidth = 1) +
  labs(
    x = "Vidutinis metilinimo lygis (beta)",
    y = "Tankis",
    color = "Regionas",
    title = "Metilinimo (beta) reikšmių pasiskirstymas pagal CpG regioną"
  ) +
  theme_minimal()
# Grafike matoma, kad salų regionai yra mažiausiai metilinti, ir didžioji jų dalis yra nemetilinta.
# Salų krantai yra pusiau nemetilinti ir metilinti.
# Salų šleifai yra daugiau metilinti nei krantai ar pačios salos.
# Atviros jūros regionai yra panašiai metilinti kaip šleifai - didžiausia dalis yra metilinta.
# Grafikas patvirtina, kad CpG salos yra nemetilintos, daugiau yra krantai, o daugiausiai
# šleifai ir atviros jūros regionai. Vadinasi mūsų mėginiai atitinka kokybės kontrolę.

# Kokybės kontrolė #2

# Apskaičiuojame koreliacijas
cor_matrix <- cor(data, use = "pairwise.complete.obs")

# Paimame visus stulpelius
sample_ann <- colanns(data)

# Paimame tik apatinį trikampį (be diagonalės)
lower <- which(lower.tri(cor_matrix), arr.ind = TRUE)
pair_cor <- cor_matrix[lower]

# Sukuriame data.frame, kur sužymime, kurie mėginiai yra tarp, o kurie grupės viduje.
# Bus naudojama grafike.
df <- data.frame(
  cor = pair_cor,
  sex_same = ifelse(sample_ann$sex[lower[,1]] == sample_ann$sex[lower[,2]], "Grupės viduje", "Tarp grupių"),
  donor_same = ifelse(sample_ann$donor[lower[,1]] == sample_ann$donor[lower[,2]], "Grupės viduje", "Tarp grupių"),
  diet_same = ifelse(sample_ann$diet[lower[,1]] == sample_ann$diet[lower[,2]], "Grupės viduje", "Tarp grupių"),
  stimulus_same = ifelse(sample_ann$stimulus[lower[,1]] == sample_ann$stimulus[lower[,2]], "Grupės viduje", "Tarp grupių")
)

# Sukuriame visų analizuojamų grupių grafikai
p1 <- ggplot(df, aes(x = cor, color = sex_same)) + geom_density(linewidth = 1) +
  labs(title = "Lytis", x = "Koreliacija", color = "Ta pati grupė") + theme_minimal()

p2 <- ggplot(df, aes(x = cor, color = donor_same)) + geom_density(linewidth = 1) +
  labs(title = "Donoras", x = "Koreliacija", color = "Ta pati grupė") + theme_minimal()

p3 <- ggplot(df, aes(x = cor, color = diet_same)) + geom_density(linewidth = 1) +
  labs(title = "Dieta", x = "Koreliacija", color = "Ta pati grupė") + theme_minimal()

p4 <- ggplot(df, aes(x = cor, color = stimulus_same)) + geom_density(linewidth = 1) +
  labs(title = "Aktyvumas", x = "Koreliacija", color = "Ta pati grupė") + theme_minimal()

p1 / p2 / p3 / p4

# Lyties ir Donoro grafikuose matoma, kad panašumas grupės viduje yra didesnis nei tarp skirtingų grupių.
# Tai yra rezultatas, kurio ir tikėjomės, vadinasi kokybės kontrolę atitinka.
# Papildomai patikrinta Dietos ir Aktyvumo grupės. Kadangi šios grupės nebuvo tarp tų, kurios būtų žinomos,
# kad DNR modifikacija skiriasi, matome kad tiek grupės viduje, tiek tarp grupių yra vienodas panašumas.

# Išskirčių paieška

# Išskirčių paieškai bus naudojamas Inter-Array Correlation metodas
# Pirmas žingsnis: koreliacijų matrica
# Koreliacijų matrica jau apskaičiuota: cor_matrix

# Kiekvienam mėginiui apskaičiuojame vidutinę koreliaciją su kitais mėginiais.
# Norint teisingai apskaičiuoti vidutinę koreliaciją su kitais mėginiais,
# koreliacijų matricoje turime ignoruoti mėginių koreliacijas su savimi, t.y. matricos diagonales reikšmes 1.
# Matricos diagonalę pakeičiame nuliais, o skaičiuojant vidurkį iš eilučių kiekio atimame vieną.
cor_matrix_adj <- cor_matrix
diag(cor_matrix_adj) <- 0
col_mean_cor <- colSums(cor_matrix_adj) / (nrow(cor_matrix_adj) - 1)
col_mean_cor

# Trečias žingsnis: atrinkti mėginius, kurių vidutinė koreliacija daugiau negu trimis standartiniais
# nuokrypiais mažesnė už vidutinę.
mean_cor <- mean(col_mean_cor)
sd_cor <- sd(col_mean_cor)
outliers <- col_mean_cor
outliers[outliers < mean_cor - 3 * sd_cor]

# Gautos išskirtys: 182_CENTRALT0, 198_CENTRAL_T0, 23_CENTRAL_T18, 245_CENTRAL_T0
# Susikuriame laikiną duomenų lentelę su pašalintomis išimtimis.
outlier_names = c("182_CENTRAL_T0","198_CENTRAL_T0","23_CENTRAL_T18","245_CENTRAL_T0")
data_iteration1 <- data[, !colnames(data) %in% outlier_names]

# Atliekame antrą metodo iteraciją.
cor_matrix <- cor(data_iteration1, use = "pairwise.complete.obs")

cor_matrix_adj <- cor_matrix
diag(cor_matrix_adj) <- 0
col_mean_cor <- colSums(cor_matrix_adj) / (nrow(cor_matrix_adj) - 1)
col_mean_cor

mean_cor <- mean(col_mean_cor)
sd_cor <- sd(col_mean_cor)
outliers <- col_mean_cor
outliers[outliers < mean_cor - 3 * sd_cor]

# Gautos išskirtys: 144_CENTRAL_T0, 175_CENTRAL_T0, 18_CENTRAL_T0, 18_CENTRAL_T18, 266_CENTRAL_T18
outlier_names = c("144_CENTRAL_T0","175_CENTRAL_T0","18_CENTRAL_T0","18_CENTRAL_T18","266_CENTRAL_T18")
data_iteration2 <- data_iteration1[, !colnames(data_iteration1) %in% outlier_names]

# Atliekame trečią metodo iteraciją.
cor_matrix <- cor(data_iteration2, use = "pairwise.complete.obs")

cor_matrix_adj <- cor_matrix
diag(cor_matrix_adj) <- 0
col_mean_cor <- colSums(cor_matrix_adj) / (nrow(cor_matrix_adj) - 1)
col_mean_cor

mean_cor <- mean(col_mean_cor)
sd_cor <- sd(col_mean_cor)
outliers <- col_mean_cor
outliers[outliers < mean_cor - 3 * sd_cor]

# Gautos išskirtys: 126_CENTRAL_T0, 175_CENTRAL_T18, 233_CENTRAL_T18, 245_CENTRAL_T18, 295_CENTRAL_T0
outlier_names = c("126_CENTRAL_T0","175_CENTRAL_T18","233_CENTRAL_T18","245_CENTRAL_T18","295_CENTRAL_T0")

# Visų outlierių sąrašas.
all_outliers <- c("182_CENTRAL_T0","198_CENTRAL_T0","23_CENTRAL_T18",
                  "245_CENTRAL_T0","144_CENTRAL_T0","175_CENTRAL_T0",
                  "18_CENTRAL_T0","18_CENTRAL_T18","266_CENTRAL_T18",
                  "126_CENTRAL_T0","175_CENTRAL_T18","233_CENTRAL_T18",
                  "245_CENTRAL_T18","295_CENTRAL_T0")

# Padaliname duomenis: vien tik išskirčių duomenų rinkinys ir duomenų rinkinys be išskirčių.
outlier_data <- data[,colnames(data) %in% all_outliers]
data_without_outliers <- data[,!colnames(data) %in% all_outliers]

# Norint vizualiai atvaizduoti išskirčių mėginius, sudarome data frame iš išskirčių duomenų rinkinio.
# Pridedame papildomą stulpelį, kur pažymėtas mėginio sample vardas (jį gauname paėmę dalį eilutės pavadinimo).
outlier_df <- stack(outlier_data)
outlier_df$samples <- sub(".*:", "", rownames(outlier_df))

# Apskaičiuojame vidurkius eilučių reikšmėms be išskirčių.
mean_without_outliers <- rowMeans(data_without_outliers)

# Šiems vidurkiams sudarome data.frame (id paimtas kaip atsitiktinis stulpelis, kad ggplot priimtų šį df).
df <- data.frame(
  mean = mean_without_outliers,
  id = data_without_outliers@id
)

# Brėžiame bendrą grafiką išskirtims ir lyginame su vidutiniu reikšmių pasiskirstymu be išskirčių.
p5 <- ggplot() +
  geom_density(data = outlier_df, aes(x=value, group=samples, fill=samples, color=samples), alpha = 0.4) +
  geom_density(data = df, aes(x=mean), fill="white", color="black", alpha = 0.1) +
  labs(title = "Beta reikšmių pasiskirstymas", x = "Beta reikšmė", y = "Tankis", fill = "Mėginys") +
  guides(color = "none") +
  theme_minimal() 
  
p5

# Iš beta reikšmių tankio grafikų matome, jog visi išskirčių tankio grafikai yra vienodai pasiskirstę,
# nestipriai skiriasi tik išskirčių pikų aukščiai nuo vidutinio pasiskirstymo be išskirčių (juoda linija).
# Kadangi nėra stiprių skirtumų tankio grafikų formose, nėra pagrindo šalinti biologinių išskirčių.