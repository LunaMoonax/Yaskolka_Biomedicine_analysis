# Ketvirta užduotis
# Viktorija Ramonaitė, Skaistė Bartkutė

#setwd("C:/Users/Viktorija Ramonaite/Desktop/UNIVERAS/3 KURSAS/BIOMEDICINOS DUOMENU ANALIZE/1 uzduotis/Yaskolka_Biomedicine_analysis")
setwd("C:/Users/skais/Desktop/Universitetas/Šeštas semestras/Biomedicina/task1/Yaskolka_Biomedicine_analysis")

library(annmatrix)
library(dplyr) 
library(glmnet)
library(ggplot2)


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

# Ištestuojame alpha reikšmes loco ciklui intervale nuo 0.0 iki 1.0 (inkrementuojant kas 0.1)
alphas <- seq(0, 1, by = 0.1)

# Ciklas išbandyti skirtingiems alpha loco cikluose.
diff_alphas <- do.call(rbind, lapply(alphas, function(alpha_val) {

  loco <- do.call(rbind, lapply(unique(cohort), function(test_c) {
    train <- which(cohort != test_c)
    test <- which(cohort == test_c)
  
    idx <- top_cpg(train, 10000)
  
    fit <- cv.glmnet(full_beta[train, idx], age_log[train],
                     alpha = alpha_val, type.measure = "mae", nfolds = 5)
  
    pred_log <- predict(fit, full_beta[test, idx], s = "lambda.min")
  
    out <- data.frame(
      cohort = test_c,
      true   = age_raw[test],
      pred   = exp(as.vector(pred_log))
    )
  
    rm(fit, idx); gc()
    out
  }))
  
  alpha_out <- data.frame(
    alpha = alpha_val,
    mae = mae <- mean(abs(loco$pred - loco$true)),
    corrs = cor(loco$true, loco$pred)
  )
  alpha_out
}))

# Skirtingų alpha reikšmių rezultatai.
print(diff_alphas)

# alpha      mae     corrs
# 1    0.0 4.059213 0.9619025
# 2    0.1 3.533190 0.9709141
# 3    0.2 3.514222 0.9712781
# 4    0.3 3.505293 0.9714638
# 5    0.4 3.504516 0.9715078
# 6    0.5 3.505156 0.9714920
# 7    0.6 3.505789 0.9714905
# 8    0.7 3.507992 0.9714620
# 9    0.8 3.502974 0.9715427
# 10   0.9 3.504274 0.9715227
# 11   1.0 3.506883 0.9714905

# Išsirenkame geriausią alpha - mažiausią pagal MAE ir didžiausią pagal koreliaciją.
best_alpha <- diff_alphas[order(diff_alphas$mae, -diff_alphas$corrs), ][1, ]
best_alpha

# alpha      mae     corrs
# 9   0.8 3.502974 0.9715427

# Sudarome pilną modelį

# Sunumeruojame cohort eilutes 'idx' radimui.
cohort_num <- seq_along(cohort)

# Atsirenkam top CpG remdamiesi visais duomenų rinkiniais.
idx <- top_cpg(cohort_num, 10000)

# Ieškom geriausio "lambda".
fit <- cv.glmnet(full_beta[cohort_num, idx], age_log[cohort_num],
                 alpha = 0.8)

# Sudarome galutinį modelį:
best_fit <- glmnet(full_beta[cohort_num, idx], age_log[cohort_num], alpha = 0.8, lambda = fit$lambda.min)

# Kuriame funkciją, kuri paruoštų duomenis modelio pritaikymui ir pritaikytų modelį duomenims.
predict_age <- function(data_test, idx, best_fit) {
  
  # Pritaikome duomenims tas pačias transformacijas kaip su training duomenimis.
  data_preped <- prep_set(data_test, "data")
  
  beta <- as.matrix(data_preped[, -(1:2)])
  
  # Prognozuojame amžius.
  pred <- predict(best_fit, beta[,idx], s = "lambda.min")
  
  age_df <- data.frame(
    real = colanns(data_test)$age,
    id = colanns(data_test)$id
  )
  
  # Kadangi atliekant prep_test (qc), galėjo būti atmesta dalis mėginių, atmetame dalį mėginių
  # iš originalaus data_test (taip pat surikiuojame pagal tiriamuosius)
  age_df <- age_df[age_df$id %in% rownames(pred), ]
  
  # Prijungiame prognozuotą amžių
  age_df$pred <- (exp(as.vector(pred)))
  
  age_df
}

# Iškviečiame funkciją ir išsaugom rezultatus
age <- predict_age(data_test, idx, best_fit)

# Apskaičiuojame MAE
mean(abs(age$pred - age$real))
# Gauname 5.352225

# Apskaičiuojame koreliaciją
cat("Pearson r:", round(cor(age$real, age$pred), 3), "\n")
# Gauname 0.939

# Pavaizduosime realų vs. prognozuotą amžių

dim(age)
# Turime 157 tiriamuosius, kad būtų tvarkingesnis pavaizdavimas, darysime tris grafikus

ggplot(age[1:50,]) +
  geom_col(aes(x = seq_along(real), y = real, fill = "Realus")) +
  geom_col(aes(x = seq_along(pred), y = pred, fill = "Prognozuojamas"), alpha = 0.5) +
  scale_fill_manual(values = c("Realus" = "thistle4","Prognozuojamas" = "thistle")) +
  labs(title = "Realus vs. prognozuojamas amžius (rinkinys 1)",
       x = "Tiriamieji",
       y = "Metai",
       fill = "Amžius") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_blank())

ggplot(age[51:100,]) +
  geom_col(aes(x = seq_along(real), y = real, fill = "Realus")) +
  geom_col(aes(x = seq_along(pred), y = pred, fill = "Prognozuojamas"), alpha = 0.5) +
  scale_fill_manual(values = c("Realus" = "thistle4","Prognozuojamas" = "thistle")) +
  labs(title = "Realus vs. prognozuojamas amžius (rinkinys 2)",
       x = "Tiriamieji",
       y = "Metai",
       fill = "Amžius") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_blank())

ggplot(age[101:157,]) +
  geom_col(aes(x = seq_along(real), y = real, fill = "Realus")) +
  geom_col(aes(x = seq_along(pred), y = pred, fill = "Prognozuojamas"), alpha = 0.5) +
  scale_fill_manual(values = c("Realus" = "thistle4","Prognozuojamas" = "thistle")) +
  labs(title = "Realus vs. prognozuojamas amžius (rinkinys 3)",
       x = "Tiriamieji",
       y = "Metai",
       fill = "Amžius") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_blank())

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