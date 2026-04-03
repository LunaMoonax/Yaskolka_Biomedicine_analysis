# Pirma užduotis
# Viktorija Ramonaitė, Skaistė Bartkutė

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
# Susikuriame laikiną duomenų lentelę su pašalintomis išimtimis pridėdami 
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
reference[outlier_names]
reference_mean
# Išskirtys, tokios kaip: 295_T0, 295_T18, 144_T0, 18_T0, 18_T18, 266_T18, turi vidutinę koreliaciją su 
# kitais mėginiais apie 0.982 - 0.983, o vidutinė koreliacijos reikšmė yra apie 0.988, todėl tai sukelia triukšmingumą duomenyse.

# Dėl triukšmingumo pašalinami mėginiai: 295, 144, 18, 266, kurie neturėjo daugiau biologinių indikacijų.
# Mėginiai šalinami poromis dėl simetriškumo.
outlier_names = c("295_CENTRAL_T0", "295_CENTRAL_T18", "144_CENTRAL_T0", "144_CENTRAL_T18", 
                  "18_CENTRAL_T0", "18_CENTRAL_T18", "266_CENTRAL_T0", "266_CENTRAL_T18")
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

# Penktoje grupėje amžiaus intervalas siauriausias: apie 50-56 metai,
# Plačiausias amžiaus intervalas trečioje grupėje: apie 35-58 metai
# Kitose grupėse amžiaus intervalai apima maždaug nuo 40(45)-60 metų.

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