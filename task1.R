# Pirma užduotis
setwd("C:/Users/Viktorija Ramonaite/Desktop/UNIVERAS/3 KURSAS/BIOMEDICINOS DUOMENU ANALIZE/1 uzduotis/Yaskolka_Biomedicine_analysis")

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
  theme_minimal
# Grafike matoma, kad salų regionai yra mažiausiai metilinti, ir didžioji jų dalis yra nemetilinta.
# Salų krantia yra pusiau nemetilinti ir metilinti.
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
# Papildomai patikrinta Dietos ir Aktyvumo grupės. Kadangi šios grupės nebuvo tapr tų, kurios būtų žinomos,
# kad DNR modifikacija skiriasi, matome kad tiek grupės viduje, tiek tarp grupių yra vienodas panašumas.

