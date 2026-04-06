# Pirma užduotis
# Viktorija Ramonaitė, Skaistė Bartkutė

setwd("C:/Users/Viktorija Ramonaite/Desktop/UNIVERAS/3 KURSAS/BIOMEDICINOS DUOMENU ANALIZE/1 uzduotis/Yaskolka_Biomedicine_analysis")
#setwd("C:/Users/skais/Desktop/Universitetas/Šeštas semestras/Biomedicina/task1/Yaskolka_Biomedicine_analysis")

library(annmatrix)
library(ggplot2)
library(dplyr)
library(patchwork)

data <- readRDS("yaskolka.rds")

# Duomenų aprašymas

# Pasinaudojant lentelėmis ir grafikais, nagrinėsime duotus duomenis siekiant juos aprašyti.

# naudodami metaduomenis (colanns) ir CpG anotacijas (rowanns)
meta <- colanns(data)
probe_ann <- rowanns(data)

# Bendras mėginių skaičius
cat("Mėginių skaičius:", nrow(meta), "\n")

# Unikalių donorų skaičius (pašaliname pasikartojimus dėl dviejų laiko taškų)
donors <- meta[!duplicated(meta$donor), ]
cat("Unikalių donorų skaičius:", nrow(donors), "\n")

# Lyčių pasiskirstymas tarp donorų
cat("Lyčių pasiskirstymas (M - vyrai, F - moterys):")
table(donors$sex)

# Amžiaus vidurkis
meta_t0 <- meta[meta$timepoint == 0, ]
cat("Amžiaus vidurkis:", mean(meta_t0$age), "\n")

# Amžiaus pasiskirstymas pagal lytį
ggplot(meta_t0, aes(x = age, fill = sex)) +
  geom_histogram(bins = 12, alpha = 0.7) +
  facet_wrap(~sex, scales = "free_y") +
  scale_fill_manual(values = c("F" = "#E8A0BF", "M" = "#7FC8D8")) +
  labs(title = "Donorų amžiaus pasiskirstymas pagal lytį",
       x = "Amžius (metai)", y = "Skaičius") +
  theme_minimal()

# Mėginių skaičius pagal intervencijos grupę
ggplot(meta_t0, aes(x = diet, fill = stimulus)) +
  geom_bar(position = "dodge") +
  labs(title = "Donorų pasiskirstymas pagal dietą ir fizinį aktyvumą",
       x = "Dieta", y = "Donorų skaičius", fill = "Fizinis aktyvumas") +
  theme_minimal()

# Bendras matuotų CpG pozicijų skaičius 
cat("CpG pozicijų skaičius:", nrow(probe_ann), "\n")

# Iš duomenų lentelės išrinksime kiek CpG pozicijų yra kiekvienoje chromosomoje
chr_counts <- table(probe_ann$chr)
chr_df <- data.frame(chr = names(chr_counts), count = as.numeric(chr_counts))
# Surikiuojame chromosomas
chr_df$chr <- factor(chr_df$chr, levels = paste0("chr", c(1:22, "X", "Y")))
chr_df <- chr_df[!is.na(chr_df$chr), ]

# CpG pasiskirstymas tarp chromosomų
ggplot(chr_df, aes(x = chr, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "CpG pozicijų pasiskirstymas tarp chromosomų",
       x = "Chromosoma", y = "CpG skaičius") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Patikriname, ar rinkinyje yra daugiau nei vienas audinio tipas ar diagnozė
# Ląstelės tipai
cat("Ląstelės tipai:")
table(meta$celltype)

# Diagnozės tipai
cat("Diagnozės tipai:")
table(meta$diagnosis)

# Kokybės kontrolė #1

# Tikriname, ar beta reikšmių pasiskirstymas atitinka biologiškai tikėtiną modelį:
# - CpG salose (Island) metilinimo lygis turėtų būti žemas,
# - atvirose jūros srityse (OpenSea) — aukštas,
# - pakrančių (Shore) ir šelfų (Shelf) regionuose — tarpinis.

# Suskaičiuojame metilinimo lygių (beta) vidutines reikšmes eilutėms
mean_beta <- rowMeans(data, na.rm = TRUE)

# Sukuriame dataframe, kuris bus naudojama grafike
df <- data.frame(
  mean_beta = mean_beta,
  region = data@relation_to_island
)

# Suskirstome relation_to_island reikšmes į 4 grupes, kurios mums yra reikalingos kokybės kontrolei patvirtinti
df$region_grouped <- case_when(
  df$region == "Island" ~ "Island",
  df$region == "OpenSea" ~ "OpenSea",
  df$region %in% c("N_Shelf", "S_Shelf") ~ "Shelf",
  df$region %in% c("N_Shore", "S_Shore") ~ "Shore"
)
# Sugrupuojame, kad grafike būtų atvaizduota mums reikalinga tvarka
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

# Tikriname, ar mėginiai, priklausantys tai pačiai biologinei grupei, yra panašesni
# tarpusavyje nei mėginiai iš skirtingų grupių.

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

# Sukuriame visų analizuojamų grupių grafikus
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
# Reference reikšmės reikalingos vėlesniems palyginimams.
diag(cor_matrix) <- 0
col_mean_cor <- colSums(cor_matrix) / (nrow(cor_matrix) - 1)
reference <- col_mean_cor

# Trečias žingsnis: atrinkti mėginius, kurių vidutinė koreliacija daugiau negu trimis standartiniais
# nuokrypiais mažesnė už vidutinę.
mean_cor <- mean(col_mean_cor)
sd_cor <- sd(col_mean_cor)
outlier_names <- col_mean_cor[col_mean_cor < mean_cor - 3 * sd_cor]
outlier_names

reference_mean <- mean_cor

# Gautos išskirtys: 182_CENTRALT0, 198_CENTRAL_T0, 23_CENTRAL_T18, 245_CENTRAL_T0
# Susikuriame laikiną duomenų lentelę su pašalintomis išskirtimis pridėdami 
# atitinkamus papildomus mėginius pagal mėginio paėmimo laikotarpį.
outlier_names = c("182_CENTRAL_T0","182_CENTRAL_T18","198_CENTRAL_T0","198_CENTRAL_T18",
                  "23_CENTRAL_T0","23_CENTRAL_T18","245_CENTRAL_T0", "245_CENTRAL_T18")
data_without_outliers <- data[, !colnames(data) %in% outlier_names]

# Atliekame antrą metodo iteraciją.
cor_matrix <- cor(data_without_outliers, use = "pairwise.complete.obs")

diag(cor_matrix) <- 0
col_mean_cor <- colSums(cor_matrix) / (nrow(cor_matrix) - 1)
col_mean_cor

mean_cor <- mean(col_mean_cor)
sd_cor <- sd(col_mean_cor)
outlier_names <- col_mean_cor[col_mean_cor < mean_cor - 3 * sd_cor]
outlier_names

# Gautos išskirtys: 144_CENTRAL_T0, 175_CENTRAL_T0, 18_CENTRAL_T0, 18_CENTRAL_T18, 266_CENTRAL_T18, 233_CENTRAL_T18,
# Pridedame jas prie egzistuojančių išskirčių sąrašo.
outlier_names = c("144_CENTRAL_T0","144_CENTRAL_T18","175_CENTRAL_T0","175_CENTRAL_T18","18_CENTRAL_T0",
                  "18_CENTRAL_T18","266_CENTRAL_T0","266_CENTRAL_T18","233_CENTRAL_T0","233_CENTRAL_T18",
                  "182_CENTRAL_T0","182_CENTRAL_T18","198_CENTRAL_T0","198_CENTRAL_T18",
                  "23_CENTRAL_T0","23_CENTRAL_T18","245_CENTRAL_T0", "245_CENTRAL_T18")
data_without_outliers <- data[, !colnames(data) %in% outlier_names]

# Atliekame trečią metodo iteraciją.
cor_matrix <- cor(data_without_outliers, use = "pairwise.complete.obs")

diag(cor_matrix) <- 0
col_mean_cor <- colSums(cor_matrix) / (nrow(cor_matrix) - 1)
col_mean_cor

mean_cor <- mean(col_mean_cor)
sd_cor <- sd(col_mean_cor)
outlier_names <- col_mean_cor[col_mean_cor < mean_cor - 3 * sd_cor]
outlier_names

# Gautos išskirtys: 126_CENTRAL_T0, 295_CENTRAL_T0, 305_CENTRAL_T0, 305_CENTRAL_T18, 309_CENTRAL_T0
# Sudarom bendrą visų outlierių sąrašą
outlier_names = c("126_CENTRAL_T0", "126_CENTRAL_T18", "295_CENTRAL_T0", "295_CENTRAL_T18", 
                  "305_CENTRAL_T0", "305_CENTRAL_T18", "309_CENTRAL_T0", "309_CENTRAL_T18",
                  "144_CENTRAL_T0","144_CENTRAL_T18","175_CENTRAL_T0","175_CENTRAL_T18",
                  "18_CENTRAL_T0", "18_CENTRAL_T18","266_CENTRAL_T0","266_CENTRAL_T18",
                  "233_CENTRAL_T0","233_CENTRAL_T18", "182_CENTRAL_T0","182_CENTRAL_T18",
                  "198_CENTRAL_T0","198_CENTRAL_T18","23_CENTRAL_T0","23_CENTRAL_T18",
                  "245_CENTRAL_T0", "245_CENTRAL_T18")

# Padaliname duomenis: vien tik išskirčių duomenų rinkinys ir duomenų rinkinys be išskirčių.
outlier_data <- data[,colnames(data) %in% outlier_names]
data_without_outliers <- data[,!colnames(data) %in% outlier_names]

# Norint vizualiai atvaizduoti išskirčių mėginius, sudarome data frame iš išskirčių duomenų rinkinio.
# Pridedame papildomą stulpelį, kur pažymėtas mėginio sample vardas (jį gauname paėmę dalį eilutės pavadinimo).
outlier_df <- stack(outlier_data)
outlier_df$samples <- sub(".*:", "", rownames(outlier_df))

# Apskaičiuojame vidurkius eilučių reikšmėms be išskirčių.
mean_beta <- rowMeans(data_without_outliers)

# Šiems vidurkiams sudarome data.frame (id paimtas kaip atsitiktinis stulpelis, kad ggplot priimtų šį df).
df <- data.frame(
  mean = mean_beta,
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
# Dėl to atsiranda triukšmingumas duomenyse, bet pagal išskirčių pasiskirstymo formas, neatrodo, kad yra techninių klaidų.

# Išskirčių klinikinių duomenų lentelė.
outlier_df <- data.frame(
  id = outlier_data$id,
  diet = outlier_data$diet,
  stimulus = outlier_data$stimulus,
  sex = outlier_data$sex,
  age = outlier_data$age
)
outlier_df

female_samples <- data$id[as.character(data$sex) == "F"]
female_samples

# Moterų mėginiai sudaro nedidelę dalį visų mėginių.
# Tai mėginiai: 126, 175, 182, 23, 233, 235, 245, 255, 305, 309.
# Su išskirtimis sutampa: 126, 175, 182, 23, 233, 245, 305, 309.
# Tai biologinė priežastis, paaiškinanti šių išskirčių buvimą, todėl šie mėginiai bus paliekami duomenų rinkinyje.

# Likę išskirčių mėginiai: 295, 144, 18, 266, 198.

summary(data$age)
outlier_df

# Peržvelgus likusių išskirčių mėginių donorų amžius matyti, kad didžioji dalis asmenų nėra arti tyrimo amžiaus ribų,
# išskyrus mėginio 198, kur amžius yra 30.12 - 31.62, arti mažiausio amžiaus ribos 28.75, todėl šis mėginys paliekamas duomenų rinkinyje.

# Palyginamos mėginių vidutinės koreliacijos su vidutine koreliacija.
outlier_names = c("295_CENTRAL_T0", "295_CENTRAL_T18", "144_CENTRAL_T0", "144_CENTRAL_T18", 
                  "18_CENTRAL_T0", "18_CENTRAL_T18", "266_CENTRAL_T0", "266_CENTRAL_T18")
reference[outlier_names]
reference_mean
# Išskirtys, tokios kaip: 295_T0, 295_T18, 144_T0, 18_T0, 18_T18, 266_T18, turi vidutinę koreliaciją su 
# kitais mėginiais apie 0.982 - 0.983, o vidutinė koreliacijos reikšmė yra apie 0.988, todėl tai sukelia triukšmingumą duomenyse.

# Dėl triukšmingumo pašalinami mėginiai: 295, 144, 18, 266, kurie neturėjo daugiau biologinių indikacijų.
# Mėginiai šalinami poromis dėl simetriškumo.
data_without_outliers <- data[,!colnames(data) %in% outlier_names]

# Klasterizavimas

# Pirmiausia, sudaroma mėginių atstumų matrica, atstumą skaičiuojant kaip 1 - cor.
# Ši matrica paverčiama į objektą dist, reikalingą klasterizavimo funkcijai hclust.
library(WGCNA)

cor_dist_matrix <- as.dist(1 - stats::cor(data_without_outliers, use = "pairwise.complete.obs"))
clusters <- stats::hclust(cor_dist_matrix, method = "complete")

# Nubraižomas klasterizavimo grafikas (be mėginių pavadinimų norint išlaikyti tvarkingumą
# pamatyti vizualias grupes)
plotDendroAndColors(clusters, NULL, dendroLabels = FALSE)

# Dendrogramoje brėžiant liniją apytiklsiai ties viduriniu aukščiu, matomos 5 ryškiausios grupės.
groups <- cutree(clusters, k = 5)

# Gaunamas sąrašas su mėginiais ir priskirtuoju grupės numeriu, iš naujo braižomas klasterizavimo grafikas
# su spalvomis priskirtai grupei.
colors <- c("purple","red","green","blue","yellow")
dendogram_colors <- colors[groups]
plotDendroAndColors(clusters, dendogram_colors, dendroLabels = FALSE)

# Gauname vieną labai didelę grupę, tris smulkesnes grupes, ir vieną labai mažą grupę, kuri gali rodyti išskirtis.
# Padaliname medį į daugiau dalių, norint gauti grupes, kurios mažiau besiskiria savo dydžiu.
# Kartojame procesą.

groups <- cutree(clusters, k = 6)
colors <- c("purple","red","green","blue","yellow","orange")
dendogram_colors <- colors[groups]
plotDendroAndColors(clusters, dendogram_colors, dendroLabels = FALSE)

# Padalinus medį į šešias dalis, šiek tiek sumažėjo skirtumas tarp grupių dydžių, todėl toks suskirstymas
# atrodo geriau, bet vis dar matoma ta pati smulki, tikriausiai išskirčių grupė kaip prieš tai.

# Dėl įdomumo suskirstykime medį į septynias grupes.
groups <- cutree(clusters, k = 7)
colors <- c("purple","red","green","blue","yellow","orange","pink")
dendogram_colors <- colors[groups]
plotDendroAndColors(clusters, dendogram_colors, dendroLabels = FALSE)

# Toks suskirstymas nepagerina situacijos, susidaro dar viena smulki išskirčių grupė (rožinė).
# Sudarome data frame, grupių ir mėginių peržiūrai.

group_df <- data.frame(groups, dendogram_colors)
outlier_group <- group_df[group_df$dendogram_colors == "pink",]
outlier_samples <- rownames(outlier_group)
outlier_samples

# Rožinei grupei priklauso mėginių pora: 198 (arti mažiausios amžiaus ribos), tačiau šis mėginys 
# sėkmingai priskiriamas žaliai grupei, kai medis dalinamas į šešias dalis.
# Geriausias variantas: 6 grupės.

groups <- cutree(clusters, k = 6)
colors <- c("purple","red","green","blue","yellow","orange")
dendogram_colors <- colors[groups]

# Smulkiausia grupė (oranžinės spalvos), taip pat yra galimos išskirtys.
group_df <- data.frame(groups, dendogram_colors)
outlier_group <- group_df[group_df$dendogram_colors == "orange",]
outlier_samples <- rownames(outlier_group)
outlier_samples

# Tai mėginių poros: 175, 182, 309
# Šie mėginiai pateko į išskirčių sąrašą jau anksčiau, nes šių mėginių donorės buvo moterys.
# Todėl šie mėginiai biologiškai reikšmingi tyrimui.

# Atrenkame mėginius iš kiekvienos grupės.
group_one <- group_df[group_df$group == 1,]
group_one_samples <- rownames(group_one)

group_two <- group_df[group_df$group == 2,]
group_two_samples <- rownames(group_two)

group_three <- group_df[group_df$group == 3,]
group_three_samples <- rownames(group_three)

group_four <- group_df[group_df$group == 4,]
group_four_samples <- rownames(group_four)

group_five <- group_df[group_df$group == 5,]
group_five_samples <- rownames(group_five)

group_six <- group_df[group_df$group == 6,]
group_six_samples <- rownames(group_six)

# Kiekvienai grupei sudarome data frame su tos grupės klinikiniais duomenimis.
group_one_data <- data_without_outliers[,colnames(data_without_outliers) %in% group_one_samples]
group_one_df <- data.frame(
  sample_name = group_one_samples,
  diet = group_one_data$diet,
  stimulus = group_one_data$stimulus,
  timepoint = group_one_data$timepoint,
  sex = group_one_data$sex,
  age = group_one_data$age
)

group_two_data <- data_without_outliers[,colnames(data_without_outliers) %in% group_two_samples]
group_two_df <- data.frame(
  sample_name = group_two_samples,
  diet = group_two_data$diet,
  stimulus = group_two_data$stimulus,
  timepoint = group_two_data$timepoint,
  sex = group_two_data$sex,
  age = group_two_data$age
)

group_three_data <- data_without_outliers[,colnames(data_without_outliers) %in% group_three_samples]
group_three_df <- data.frame(
  sample_name = group_three_samples,
  diet = group_three_data$diet,
  stimulus = group_three_data$stimulus,
  timepoint = group_three_data$timepoint,
  sex = group_three_data$sex,
  age = group_three_data$age
)

group_four_data <- data_without_outliers[,colnames(data_without_outliers) %in% group_four_samples]
group_four_df <- data.frame(
  sample_name = group_four_samples,
  diet = group_four_data$diet,
  stimulus = group_four_data$stimulus,
  timepoint = group_four_data$timepoint,
  sex = group_four_data$sex,
  age = group_four_data$age
)

group_five_data <- data_without_outliers[,colnames(data_without_outliers) %in% group_five_samples]
group_five_df <- data.frame(
  sample_name = group_five_samples,
  diet = group_five_data$diet,
  stimulus = group_five_data$stimulus,
  timepoint = group_five_data$timepoint,
  sex = group_five_data$sex,
  age = group_five_data$age
)

group_six_data <- data_without_outliers[,colnames(data_without_outliers) %in% group_six_samples]
group_six_df <- data.frame(
  sample_name = group_six_samples,
  diet = group_six_data$diet,
  stimulus = group_six_data$stimulus,
  timepoint = group_six_data$timepoint,
  sex = group_six_data$sex,
  age = group_six_data$age
)

# Kiekvienai grupei priskiriame stulpelį su grupės pavadinimu (reikės grafikui).
group_one_df$group <- "one"
group_two_df$group <- "two"
group_three_df$group <- "three"
group_four_df$group <- "four"
group_five_df$group <- "five"
group_six_df$group <- "six"

# Visas grupes apjungiame į vieną data frame.
all_groups <- rbind(group_one_df, group_two_df, group_three_df, group_four_df, group_five_df, group_six_df)

# Braižome grafikus klinikiniams duomenims atsižvelgdami į grupes.
ggplot(all_groups, aes(x=group, fill=as.factor(diet))) +
  geom_bar(position = "dodge") +
  labs(
    x = "Grupė",
    y = "Mėginių kiekis",
    fill = "Dietos tipas",
    title = "Mėginių grupių pasiskirstymas pagal dietos tipą"
  ) +
  theme_minimal()

# Pasiskirstymas grupėse tarp dietų labai panašus, ryškesni skirtumai matomi trečioje, ketvirtoje, penktoje grupėje.
# Ketvirtoje grupėje daugiau besimaitinusių žemo angliavandenių kiekio dieta, o trečioje ir penktoje - žemo riebalų kiekio dieta.

ggplot(all_groups, aes(x=group, fill=as.factor(stimulus))) +
  geom_bar(position = "dodge") +
  labs(
    x = "Grupė",
    y = "Mėginių kiekis",
    fill = "Fizinio aktyvumo (ne)buvimas",
    title = "Mėginių grupių pasiskirstymas pagal fizinio aktyvumo (ne)buvimą"
  ) +
  theme_minimal()

# Pagal fizinio aktyvumo buvimą/nebuvimą vienodai pasiskirsčiusios pirma ir antra (beveik vienodai) grupės.
# Nedideli skirtumai pastebimi šeštoje grupėje, nes ten tik trys tiriamieji.
# Trečioje, ketvirotje grupėse daugiau asmenų, kurie užsiėmė fiziniu aktyvumu.
# Penktoje grupėje vien tik asmenys, kurie nebuvo fiziškai aktyvūs.

ggplot(all_groups, aes(x=group, fill=as.factor(timepoint))) +
  geom_bar(position = "dodge") +
  labs(
    x = "Grupė",
    y = "Mėginių kiekis",
    fill = "Mėginių paėmimo laikotarpis",
    title = "Mėginių grupių pasiskirstymas pagal mėginių paėmimo laikotarpį"
  ) +
  theme_minimal()

# Visose grupėse mėginiai pasiskirstę vienodai pagal mėginių paėmimo laikotarpius.

ggplot(all_groups, aes(x=group, fill=as.factor(sex))) +
  geom_bar(position = "dodge") +
  labs(
    x = "Grupė",
    y = "Mėginių kiekis",
    fill = "Lytis",
    title = "Mėginių grupių pasiskirstymas pagal lytį"
  ) +
  theme_minimal()

# Pirmoje, antroje, trečioje, ketvirtoje grupėse vien tik vyrai.
# Penktoje ir šeštoje grupėje vien tik moterys.

ggplot(all_groups, aes(x=group, y=age)) +
  geom_boxplot(position = "dodge") +
  labs(
    x = "Grupė",
    y = "Amžius",
    title = "Mėginių grupių pasiskirstymas pagal amžių"
  ) +
  theme_minimal()

# Penktoje grupėje kvartilinis amžiaus intervalas siauriausias: apie 50-56 metai,
# Plačiausias kvartilinis amžiaus intervalas trečioje grupėje: apie 35-58 metai
# Kitose grupėse kvartiliniai amžiaus intervalai apima maždaug nuo 40(45)-60 metų.

# Heatmap

# Apskaičiuota variacija kiekvienai duomenų matricos eilutei.
cg_variance <- apply(data_without_outliers, 1, var)

# Matrica išrikiuota mažėjimo tvarka pagal eilučių variaciją.
cg_matrix <- data_without_outliers[order(cg_variance, decreasing = TRUE),]

# Heatmap su 1000 eilučių.
heatmap(cg_matrix[1:1000,])

# Heatmap su 100 eilučių.
heatmap(cg_matrix[1:100,])

# Heatmap su 10 eilučių.
heatmap(cg_matrix[1:10,])

# Apžvalginė analizė

# PCA — Pagrindinių komponenčių analizė
# Atliekame PCA, kad nustatytume, kokie faktoriai labiausiai lemia mėginių metilinimo variaciją.

# Ttransponuojame matricą (t), nes PCA tikisi mėginius eilutėse, o CpG stulpeliuose
pca <- prcomp(t(data_without_outliers), center = TRUE, scale. = FALSE)

# Sukuriame data.frame su pirmosiomis trimis komponentėmis ir klinikiniais duomenimis
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3],
  colanns(data_without_outliers)
)

# Apskaičiuojame, kiek variacijos (%) paaiškina kiekviena komponentė
variance_explained <- summary(pca)$importance[2, 1:5] * 100
cat("Variacijos dalis (%), kurią paaiškina pirmosios 5 komponentės:\n")
print(round(variance_explained, 2))

# Grafikas pagal lytį
p_sex <- ggplot(pca_df, aes(x = PC1, y = PC2, color = sex)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("F" = "#E8A0BF", "M" = "#7FC8D8")) +
  labs(title = "Lytis",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal()

# Grafikas pagal laiko tašką (T0 vs T18)
p_time <- ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(timepoint))) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Laiko taškas", color = "Timepoint",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal()

# Grafikas pagal dietos tipą
p_diet <- ggplot(pca_df, aes(x = PC1, y = PC2, color = diet)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Dieta",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal()

# Grafikas pagal fizinį aktyvumą
p_stim <- ggplot(pca_df, aes(x = PC1, y = PC2, color = stimulus)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "Fizinis aktyvumas",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 1), "%)")) +
  theme_minimal()

# Sudedame visus keturis grafikus į vieną paveikslą
(p_sex | p_time) / (p_diet | p_stim)

# Žiūrime PC3, ar joje matosi kiti faktoriai
ggplot(pca_df, aes(x = PC1, y = PC3, color = sex)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("F" = "#E8A0BF", "M" = "#7FC8D8")) +
  labs(title = "PCA: PC1 vs PC3, spalvinta pagal lytį",
       x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PC3 (", round(variance_explained[3], 1), "%)")) +
  theme_minimal()

# Metilinimo pokytis tarp T0 ir T18
# Kiekvienam donorui apskaičiuojame metilinimo pokytį (T18 - T0)
# kiekvienoje CpG pozicijoje. Tada žiūrime, kurios pozicijos vidutiniškai
# pakito labiausiai per visus donorus.

# Pasiimame mėginių metaduomenis
meta_clean <- colanns(data_without_outliers)

# Išskiriame T0 ir T18 donorus
t0_donors <- meta_clean$donor[meta_clean$timepoint == 0]
t18_donors <- meta_clean$donor[meta_clean$timepoint == 18]

# Surandame donorus, kurie turi abu laiko taškus
common_donors <- intersect(t0_donors, t18_donors)

# Surandame stulpelių indeksus atitinkamiems T0 ir T18 mėginiams
t0_idx <- match(paste0(common_donors, "_T0"), colnames(data_without_outliers))
t18_idx <- match(paste0(common_donors, "_T18"), colnames(data_without_outliers))

# Skaičiuojame vidutinį delta beta (T18 - T0) kiekvienai CpG pozicijai per visus donorus
delta_beta <- rowMeans(data_without_outliers[, t18_idx] - data_without_outliers[, t0_idx], na.rm = TRUE)

# Vizualizuojame delta beta pasiskirstymą — tikimės, kad dauguma pokyčių bus arti 0
df_delta <- data.frame(delta = delta_beta)

ggplot(df_delta, aes(x = delta)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +  # raudona linija ties 0
  labs(title = "Metilinimo pokytis (T18 - T0) per visas CpG pozicijas",
       x = "Vidutinis beta pokytis (delta)", y = "CpG pozicijų skaičius") +
  theme_minimal()

# Amžiaus koreliacija su metilinimu (epigenetinis laikrodis)

# ikriname, ar yra CpG pozicijų, kurių bazinis metilinimo lygis (T0)
# stipriai koreliuoja su donoro amžiumi.

# Atskiriame tik T0 mėginius (baziniai duomenys prieš intervenciją)
data_t0 <- data_without_outliers[, meta_clean$timepoint == 0]

# Paimame donorų amžius
ages <- colanns(data_t0)$age

# Skaičiuojame Pearson koreliaciją tarp kiekvienos CpG pozicijos beta reikšmės ir amžiaus
cor_with_age <- apply(data_t0, 1, function(x) cor(x, ages, use = "pairwise.complete.obs"))

# Vizualizuojame koreliacijos pasiskirstymą — tikimės normalaus pasiskirstymo apie 0
# su keliais stipriau koreliuojančiais CpG
df_age_cor <- data.frame(correlation = cor_with_age)

ggplot(df_age_cor, aes(x = correlation)) +
  geom_histogram(bins = 100, fill = "darkorange", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "CpG metilinimo koreliacija su amžiumi (T0)",
       x = "Pearson koreliacija", y = "CpG pozicijų skaičius") +
  theme_minimal()
