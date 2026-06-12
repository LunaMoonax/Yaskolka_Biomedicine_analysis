# Ketvirta užduotis
# Viktorija Ramonaitė, Skaistė Bartkutė

setwd("C:/Users/Viktorija Ramonaite/Desktop/UNIVERAS/3 KURSAS/BIOMEDICINOS DUOMENU ANALIZE/1 uzduotis/Yaskolka_Biomedicine_analysis")
#setwd("C:/Users/skais/Desktop/Universitetas/Šeštas semestras/Biomedicina/task1/Yaskolka_Biomedicine_analysis")

library(annmatrix)
library(dplyr) 
library(glmnet)


# Nuskaitome duomenų rinkinius
data_train_1 <- readRDS("AgeData/fraga.rds")

#colanns(data_train_1)
#rowanns(data_train_1)

data_train_2 <- readRDS("AgeData/johansson.rds")
data_train_3 <- readRDS("AgeData/kalyakulina.rds")
data_train_4 <- readRDS("AgeData/kurushima.rds")
data_train_5 <- readRDS("AgeData/mahdi.rds")
data_train_6 <- readRDS("AgeData/xu.rds")

# quinn.rds pasiliekame kaip testuojamą duomenų rinkinį.
data_test <- readRDS("AgeData/quinn.rds")

# Testuojamų rinkinių sąrašas
train_sets <- list(
  fraga       = data_train_1,
  johansson   = data_train_2,
  kalyakulina = data_train_3,
  kurushima   = data_train_4,
  mahdi       = data_train_5,
  xu          = data_train_6
)

lapply(train_sets, function(x) colnames(colanns(x)))

# Surenkame tik bendrus CpG visiems duomenų rinkiniams
common_cpg <- Reduce(intersect, lapply(train_sets, rownames))

# Funkcija paruošti duomenų rinkinį išflitravus išskirtis
# ir patogiam duomenų formate
prep_set <- function(x, cohort_name) {
  ca <- colanns(x)
  
  # QC išskirčių taisyklės
  drop <- coalesce(ca$qc_iac < -3,        FALSE) |
    coalesce(ca$qc_detection < 0.95, FALSE) |
    coalesce(ca$qc_badsex == TRUE,   FALSE) |
    coalesce(ca$qc_badsnp == TRUE,   FALSE)
  
  keep <- !drop
  
  # beta matrica, apkarpyta iki bendrų CpG, transponuota -> mėginiai × CpG
  beta <- t(x[common_cpg, keep])
  
  data.frame(
    cohort = cohort_name,
    age    = ca$age[keep],
    beta,
    check.names = FALSE
  )
}

# Sujungiam visus 6 duomenų runkinius į vieną aibę
# Duomenų rinkiniai tuo pačiu yra paruošiami
full <- do.call(rbind, Map(prep_set, train_sets, names(train_sets)))

# Patikriname
#dim(full)
#table(full$cohort)
#sum(is.na(full$age))
# buvo 7 reikšmės NA

# Išmetame mėginius, kur age yra NA
full <- full[!is.na(full$age), ]

# Patikriname ar yra trūkstamų rekšmių kitur, nes glmnet metodui netinka jei yra NA
#beta_cols <- full[, -(1:2)]
#sum(is.na(beta_cols))
# 0 reiksmiu

# Beta reikšmes padarome matrix, kad būtų lengviau saičiuoti
full_beta <- as.matrix(full[, -(1:2)])
# Išsaugome netransformuotus age reikšmes
age_raw <- full$age
# Logoritmuojame amžių
age_log <- log(full$age)

cohort  <- full$cohort

# full nebereikalingas kaip didelė matrica – atlaisvinam atmintį
rm(full)
# priverstinai atlaisviname atmintį
gc()

# Grąžina N CpG indeksų, labiausiai susijusių su amžiumi
top_cpg <- function(rows, n_keep = 10000) {
  cors <- as.vector(cor(age_log[rows], full_beta[rows, ]))
  cors[is.na(cors)] <- 0
  order(abs(cors), decreasing = TRUE)[seq_len(n_keep)]
}

# LOCO ciklas - leave-one-cohort-out
# Testuojame ir mokome modelį - tai ne galutinis modelis.
# Naudosime elastic-net modeli is bilbiotekos glmnet, pagal literatura tai turetu but tiksliausias.
loco <- do.call(rbind, lapply(unique(cohort), function(test_c) {
  train <- which(cohort != test_c)
  test <- which(cohort == test_c)
  
  idx <- top_cpg(train, 10000)
  
  fit <- cv.glmnet(full_beta[train, idx], age_log[train],
                   alpha = 0.5, type.measure = "mae", nfolds = 5)
  
  pred_log <- predict(fit, full_beta[test, idx], s = "lambda.min")
  
  out <- data.frame(
    cohort = test_c,
    true   = age_raw[test],
    pred   = exp(as.vector(pred_log))
  )
  
  rm(fit, idx); gc()
  out
}))

# klaida metais
loco$err <- loco$pred - loco$true

# Bendras MAE (vidutinė absoliuti klaida metais)
mae <- mean(abs(loco$err))
cat("LOCO MAE (metais):", round(mae, 2), "\n")
# LOCO MAE (metais): 3.51 

# Klaida pagal rinkinį - ar koks nors iškrenta
aggregate(abs(err) ~ cohort, data = loco, FUN = mean)
#       cohort abs(err)
#1       fraga 3.515834
#2   johansson 3.700486
#3 kalyakulina 4.409006
#4   kurushima 3.392192
#5       mahdi 3.443585
#6          xu 2.598103

# Pearson koreliacija tikras vs spėtas - kaip gerai seka
cat("Pearson r:", round(cor(loco$true, loco$pred), 3), "\n")
# Pearson r: 0.971 


# Išvada kol kas:
# Mokomas modelis jau labai tikslus pagal literatūra:
# epigenetiniai laikrodziai kaip Horvath duoda MAE tarp 3-5 metu
# koreliacija jei virs 0.95 - labai gera

# Ka dar reikia padaryt:
# Galima pabandyt dar pakoreguot modelio apmokyma su alpha
# Siuo metu alfa parinktas random 0.5, bet derinant (mazinant ar aukstinant) butu galima
# gal dar tikslesni rezultata gaut ir parodyt kaip buvo derinama
# Pabaigti ir sukurti galutini modeli ant 5 duomenu rinkiniu, padaryti is jo funkcija ir istestuoti ant
# palikto duomenu rinkinio quinn