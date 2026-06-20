################################################################################
# Project: Competing Risks Analysis and Dynamic Prediction for Esophageal Cancer
# Author：Huang’s Lab
# Date: 2026-01-10
# System Requirements:
#   - R Version: 4.4.3
################################################################################

# ==============================================================================
# Part 1: Data Cleaning and Standardization
# ==============================================================================

####Standardize original values
setwd("D:/SEER_ECSS/Data_Processing")
library(dplyr)
library(naniar)
library(tidyr)
library(readr)

data <- read.table("export.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)

str(data)

data$age <- case_when(
  data$age %in% c("20-24 years","25-29 years","30-34 years","35-39 years",
                  "40-44 years","45-49 years") ~ "<50",
  data$age %in% c("50-54 years","55-59 years") ~ "50-59",
  data$age %in% c("60-64 years","65-69 years") ~ "60-69",
  data$age %in% c("70-74 years","75-79 years") ~ "70-79",
  data$age %in% c("80-84 years","85-89 years","90+ years") ~ "≥80",
  TRUE ~ NA_character_
)

data$age <- factor(data$age,
                   levels=c("<50","50-59","60-69","70-79","≥80"))

data$sex <- factor(data$sex, levels=c("Male","Female"))

data$race <- factor(data$race,
                    levels=c("White","Asian or Pacific Islander",
                             "Black","American Indian/Alaska Native"))

data$primary_site <- case_when(
  grepl("C15.0|C15.1|C15.2|C15.3", data$primary_site) ~ "Upper",
  grepl("C15.4", data$primary_site) ~ "Middle",
  grepl("C15.5", data$primary_site) ~ "Lower",
  grepl("C15.8|C15.9", data$primary_site) ~ "Other",
  TRUE ~ NA_character_
)

data$primary_site <- factor(data$primary_site,
                            levels=c("Upper","Middle","Lower","Other"))
g <- data$grade

data$grade <- case_when(
  grepl("Grade I$",   g) ~ "I",
  grepl("Grade II$",  g) ~ "II",
  grepl("Grade III$", g) ~ "III",
  grepl("Grade IV$",  g) ~ "IV",
  grepl("Unknown",    g, ignore.case=TRUE) ~ NA_character_,
  TRUE ~ NA_character_
)

data$grade <- factor(data$grade,
                     levels=c("I","II","III","IV"))

data$ajcc_stage <- case_when(
  grepl("IVA|IVB|IVNOS", data$ajcc_stage) ~ "IV",
  grepl("III",  data$ajcc_stage) ~ "III",
  grepl("IIA|IIB", data$ajcc_stage) ~ "II",
  grepl("^I$|^I ", data$ajcc_stage) ~ "I",
  TRUE ~ NA_character_
)

data$ajcc_stage <- factor(data$ajcc_stage, levels=c("I","II","III","IV"))

data$ajcc_t <- ifelse(data$ajcc_t %in% c("T0","T1","T2","T3","T4"),
                      data$ajcc_t, NA)
data$ajcc_t <- factor(data$ajcc_t, levels=c("T0","T1","T2","T3","T4"))

data$ajcc_n <- ifelse(data$ajcc_n %in% c("N0","N1"), data$ajcc_n, NA)
data$ajcc_n <- factor(data$ajcc_n, levels=c("N0","N1"))

data$ajcc_m <- case_when(
  data$ajcc_m == "M0" ~ "M0",
  data$ajcc_m %in% c("M1","M1a","M1b","M1NOS") ~ "M1",
  TRUE ~ NA_character_
)
data$ajcc_m <- factor(data$ajcc_m, levels=c("M0","M1"))

data$surgery <- case_when(
  data$surgery == 0 ~ "No",
  data$surgery >=10 & data$surgery <= 90 ~ "Yes",
  data$surgery %in% c(98,99) ~ NA_character_,
  TRUE ~ NA_character_
)
data$surgery <- factor(data$surgery, levels=c("No","Yes"))

data$chemo <- case_when(
  data$chemo == "Yes" ~ "Yes",
  data$chemo == "No/Unknown" ~ "No",
  TRUE ~ NA_character_
)
data$chemo <- factor(data$chemo, levels=c("No","Yes"))

data$radiation <- case_when(
  grepl("Beam|brachytherapy|isotopes|NOS  method", data$radiation, ignore.case=TRUE) ~ "Yes",
  grepl("None/Unknown|Refused", data$radiation, ignore.case=TRUE) ~ "No",
  TRUE ~ NA_character_
)
data$radiation <- factor(data$radiation, levels=c("No","Yes"))

special_codes <- c(989:999)
data$tumor_size[data$tumor_size %in% special_codes |
                  data$tumor_size > 150] <- NA

data$tumor_size <- data$tumor_size / 10 

data$tumor_size_group <- cut(
  data$tumor_size,
  breaks=c(-Inf,2,5,10,Inf),
  labels=c("<2 cm","2-5 cm","5-10 cm","≥10 cm"),
  right=FALSE
)

data$tumor_size_group <- factor(data$tumor_size_group,
                                levels=c("<2 cm","2-5 cm","5-10 cm","≥10 cm"))

data$cause_of_death <- as.factor(data$cause_of_death)

data$vital_status  <- as.factor(data$vital_status)

data$fstatus <- case_when(
  data$vital_status == "Alive" ~ "Alive",
  data$vital_status == "Dead" & data$cause_of_death == "Esophagus" ~ "Esophagus",
  data$vital_status == "Dead" & data$cause_of_death != "Esophagus" ~ "Other Cause",
  TRUE ~ NA_character_
)

data$fstatus <- factor(data$fstatus,
                       levels=c("Alive","Esophagus","Other Cause"))

data$survival_months[data$survival_months=="Unknown"] <- NA
data$survival_months <- as.numeric(data$survival_months)

data$histology <- as.factor(data$histology)

# Check Standardization Results
str(data)

for(col in names(data)){
  cat("Column:", col, "\n")
  print(unique(data[[col]]))
  cat("\n")
}

#### Missing Value Visualization
miss_var_summary(data)

#### Missing Value Handling (Main Analysis)
data <- data

add_unknown <- function(v){
  if (!is.factor(v)) return(v)
  lv <- levels(v)
  x <- as.character(v)
  x[is.na(x)] <- "Unknown"
  factor(x, levels = c(lv,"Unknown"))
}

vars_add_unknown <- c("grade","ajcc_stage","ajcc_t","ajcc_n","ajcc_m",
                      "radiation","surgery","tumor_size_group")

data_main[vars_add_unknown] <- lapply(data_main[vars_add_unknown], add_unknown)

data_main <- data_main %>% filter(!is.na(survival_months))

saveRDS(data_main, "data_main.rds")
write.csv(data_main, "data_main.csv", row.names = FALSE, fileEncoding = "GB18030")

for (v in names(data_main)) {
  if (is.factor(data_main[[v]])) {
    cat("\nVariable:", v, "\n")
    print(levels(data_main[[v]]))
  }
}

#### Missing Value Handling (Sensitivity Analysis)
data_no_na <- data %>% drop_na()

saveRDS(data_no_na, "data_no_na.rds")
write.csv(data_no_na, "data_no_na.csv", row.names = FALSE, fileEncoding = "GB18030")

for (v in names(data_no_na)) {
  if (is.factor(data_no_na[[v]])) {
    cat("\nVariable:", v, "\n")
    print(levels(data_no_na[[v]]))
  }
}

# ==============================================================================
# Part 2: Main Analysis
# ==============================================================================

#### Stratified Main Analysis: Build Fine-Gray Model
setwd("D:/SEER_ECSS/ECSS_final")

library(tidyverse)
library(dplyr)
library(cmprsk)
library(riskRegression)
library(regplot)
library(rms)
library(pec)
library(prodlim)
library(openxlsx)
library(officer)
library(flextable)
library(stringr)
library(knitr)
library(rmda)
library(ggplot2)
library(prodlim)
library(patchwork)
library(patchwork)
library(jpeg)
library(grid)

data <- readRDS("data_main.rds")
cat("Initial rows after reading data:", nrow(data), "\n")

data <- data %>% filter(survival_months > 0)
cat("Rows after filtering survival > 0:", nrow(data), "\n")  

data$fstatus <- recode(as.character(data$fstatus),
                       "Alive" = 0,
                       "Esophagus" = 1,
                       "Other Cause" = 2)
data$fstatus <- as.numeric(data$fstatus)

table(data$fstatus)

### Stratified Sampling
set.seed(2025)
data <- data %>% mutate(row_id = row_number())

train_idx <- data %>%
  group_by(fstatus) %>%
  sample_frac(0.7) %>%
  ungroup() %>%
  pull(row_id)

train_data <- data[data$row_id %in% train_idx, ]
test_data  <- data[!data$row_id %in% train_idx, ]

prop.table(table(train_data$fstatus))
prop.table(table(test_data$fstatus))



### Figure S1: Plot Overall CIF 
### Plot Overall CIF
ci_overall <- cuminc(
  ftime = data$survival_months, 
  fstatus = data$fstatus,
  cencode = 0
)

overall1 <- ci_overall[["1 1"]]

overall2 <- ci_overall[["1 2"]]

overall_df <- bind_rows(
  data.frame(time = overall1$time, est = overall1$est, event = "Event 1: Esophagus"),
  data.frame(time = overall2$time, est = overall2$est, event = "Event 2: Other Cause")
)

p_overall <- ggplot(overall_df, aes(x = time, y = est, color = event)) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(
    values = c("Event 1: Esophagus" = "#2E86AB",
               "Event 2: Other Cause" = "#E76F51")
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Months",
    y = "Cumulative Incidence",
    color = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = c(0.95, 0.08),
    legend.justification = c(1, 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18)
  ) +
  coord_cartesian(xlim = c(0, 60), ylim = c(0, NA))


ggsave(
  filename = "Figure_S1_CIF_overall.pdf",  
  plot = p_overall,
  width = 8,
  height = 6,
  units = "in",
  bg = "white",
  device = cairo_pdf  
)


### Figure S2: CIF for Important Variables
important_vars <- c("age", "ajcc_stage", "grade", "chemo", "radiation", "surgery")

for (i in seq_along(important_vars)) {
  var <- important_vars[i]
  
  ci_var <- cuminc(
    ftime = data$survival_months, 
    fstatus = data$fstatus, 
    group = data[[var]], 
    cencode = 0
  )
  
  ci_df_var <- data.frame() 
  for (g in unique(data[[var]])) {
    nm <- paste0(g, " 1")
    if (!nm %in% names(ci_var)) next
    tmp <- ci_var[[nm]]
    if (is.null(tmp$time)) next
    
    ci_df_var <- rbind(ci_df_var, data.frame(
      time = tmp$time,
      est = tmp$est,
      group = g
    ))
  }
  ci_df_var <- ci_df_var %>% drop_na()
  
  p_var <- ggplot(ci_df_var, aes(x = time, y = est, color = group)) +
    geom_step(linewidth = 1.2) +
    labs(x = "Months", y = "Cumulative Incidence", color = NULL) +
    theme_minimal(base_size = 18, base_family = "Arial") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position = c(0.92, 0.08),
      legend.justification = c(1, 0),
      legend.direction = "vertical",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.6, "lines"),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.line  = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      axis.title.x = element_text(face = "bold", size = 18),
      axis.title.y = element_text(face = "bold", size = 18),
      axis.text.x  = element_text(size = 16),
      axis.text.y  = element_text(size = 16)
    ) +
    coord_cartesian(xlim = c(0, 60)) +
    annotate("text", x = 0.5, y = max(ci_df_var$est)*1.05,
             label = paste0("(", labels[i], ")"), fontface = "bold", size = 6, hjust = 0)
  
  ggsave(
    filename = paste0("Figure_S2_CIF_", var, ".pdf"),
    plot = p_var,
    width = 8,
    height = 6,
    units = "in",
    bg = "white",
    device = cairo_pdf
  )
}


###Build Fine-Gray Model (Main Model)
fg_fit <- FGR(                                  
  formula = Hist(survival_months, fstatus) ~ 
    age + sex + race + primary_site + 
    ajcc_stage +grade + tumor_size_group + 
    surgery + chemo + radiation +year_dx, 
  data = train_data, 
  cause = 1 
)

print(fg_fit)
save(fg_fit, file = "fg_fit_train.RData")

load("fg_fit_train.RData") #load


### Table S1: Extract Fine-Gray Model Results
coef <- fg_fit$crrFit$coef
se <- sqrt(diag(fg_fit$crrFit$var))

HR <- exp(coef)
Lower95 <- exp(coef - 1.96 * se)
Upper95 <- exp(coef + 1.96 * se)
p_value <- 2 * pnorm(-abs(coef / se))
p_value <- ifelse(p_value < 0.001, "<0.001", 
                  ifelse(p_value < 0.05, "<0.05", round(p_value, 3)))

var_names <- names(coef)


get_base_var_and_ref <- function(varname, data) {
  
  for (v in names(data)) {
    if (is.factor(data[[v]])) {
      prefix <- paste0(v)
      if (startsWith(varname, prefix)) {
        lv <- levels(data[[v]])
        return(c(base_var = v, ref = lv[1]))
      }
    }
  }
  
  return(c(base_var = NA, ref = ""))
}

match_info <- t(sapply(var_names, get_base_var_and_ref, data = train_data))
Reference <- match_info[, "ref"]

result_df <- data.frame(
  Variable = var_names,
  Reference = Reference,
  HR = round(HR, 2),
  Lower95 = round(Lower95, 2),
  Upper95 = round(Upper95, 2),
  p_value = p_value,
  stringsAsFactors = FALSE
)

kable(result_df, caption = "Fine-Gray Model Results", align = "lcccc")

write.xlsx(result_df, "FineGray_train_results_main.xlsx", rowNames = FALSE)



### Figure 1: Nomogram
train_pred <- predictRisk(fg_fit, newdata = train_data, times = c(12, 36, 60))
colnames(train_pred) <- paste0("pred_CIF_", c(12, 36, 60), "m")
data_train_pred <- bind_cols(train_data, as.data.frame(train_pred)) %>% 
  mutate(event36 = as.integer(survival_months <= 36 & fstatus == 1))

dd <- datadist(data_train_pred)
options(datadist = "dd")

f_cph <- lrm(
  event36 ~ age + sex + race + primary_site + tumor_size_group + 
    grade + ajcc_stage + surgery + chemo + radiation,
  data = data_train_pred
)

lp <- predict(f_cph, type = "lp")
fit_fun <- lm(qlogis(data_train_pred$pred_CIF_36m) ~ lp)

fun_1y <- function(x) plogis(coef(fit_fun)[1] + coef(fit_fun)[2] * x * (12/36))
fun_3y <- function(x) plogis(coef(fit_fun)[1] + coef(fit_fun)[2] * x)
fun_5y <- function(x) plogis(coef(fit_fun)[1] + coef(fit_fun)[2] * x * (60/36))

nom <- nomogram(
  f_cph,
  fun = list(fun_1y, fun_3y, fun_5y),
  funlabel = c("1-year risk", "3-year risk", "5-year risk"),
  lp = TRUE,
  maxscale = 120,
  fun.at = seq(0.05, 0.95, 0.1)
)

jpeg(
  filename = "Figure_1_Nomogram.jpg",
  width = 18,
  height = 16,
  units = "in",  
  res = 300,     
  quality = 100,  
  bg = "white",
  family = "sans",
  pointsize = 16  
)

plot(
  nom,
  xfrac = 0.25,  
  cex.var = 4, 
  cex.axis = 2, 
  lmgp = 0.4,   
  tck = 0.6,    
  height = 2.5, 
  lwd = 3.2     
)

dev.off()



### Figure 2: AUC, Brier Score
score_res <- Score(
  list(FG = fg_fit),
  formula = Hist(survival_months, fstatus) ~ 1,
  data = test_data,
  metrics = c("AUC", "Brier"),
  times = c(12, 36, 60),
  plots = "calibration"
)

times <- c(12, 36, 60)

score_result <- Score(
  list("FG" = fg_fit),
  data = test_data,
  formula = Hist(survival_months, fstatus) ~ 1,
  times = times,
  plots = "calibration",
  metrics = c("auc", "brier"),
  split.seed = 123,
  B = 500
)

auc_df <- as.data.frame(score_result$AUC$score) %>%
  filter(model == "FG") %>%
  rename(value = AUC) %>%
  mutate(
    lower = value - 1.96 * se,
    upper = value + 1.96 * se,
    metric = "AUC"
  )

brier_df <- as.data.frame(score_result$Brier$score) %>%
  filter(model == "FG") %>%
  rename(value = Brier) %>%
  mutate(
    lower = value - 1.96 * se,
    upper = value + 1.96 * se,
    metric = "Brier"
  )

plot_df <- bind_rows(auc_df, brier_df)

p_auc <- ggplot(
  subset(plot_df, metric == "AUC"),
  aes(x = times, y = value)
) +
  geom_line(color = "#808080", linewidth = 1.5) +
  geom_point(color = "#1F77B4", size = 3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 2,
    alpha = 0.5,
    color = "#A9A9A9",
    linewidth = 1.0
  ) +
  labs(
    x = "Time (months)",
    y = "Time-dependent AUC (95% CI)"
  ) +
  scale_x_continuous(breaks = seq(0, 60, by = 12), expand = c(0, 0)) +
  scale_y_continuous(
    labels = scales::comma_format(accuracy = 0.01),
    limits = c(0.55, 0.8),
    expand = c(0, 0),
    breaks = seq(0.55, 0.8, by = 0.05)
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    axis.line    = element_line(color = "black", linewidth = 0.8),
    axis.ticks   = element_line(color = "black", linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(0, 62), ylim = c(0.55, 0.83))

p_brier <- ggplot(
  subset(plot_df, metric == "Brier"),
  aes(x = times, y = value)
) +
  geom_line(color = "#808080", linetype = "dashed", linewidth = 1.5) +
  geom_point(shape = 1, color = "#1F77B4", size = 3) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 2,
    alpha = 0.5,
    color = "#A9A9A9",
    linewidth = 1.0
  ) +
  labs(
    x = "Time (months)",
    y = "Brier score (95% CI)"
  ) +
  scale_x_continuous(breaks = seq(0, 60, by = 12), expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0.15, 0.25),
    labels = scales::comma_format(accuracy = 0.01),
    expand = c(0, 0),
    breaks = seq(0.15, 0.25, by = 0.05)
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    axis.line    = element_line(color = "black", linewidth = 0.8),
    axis.ticks   = element_line(color = "black", linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(0, 62), ylim = c(0.15, 0.28))

ggsave(
  filename = "Figure_2A_AUC.pdf",
  plot = p_auc,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white",
  device = "pdf"
)

ggsave(
  filename = "Figure_2B_Brier.pdf",
  plot = p_brier,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  bg = "white",
  device = "pdf"
)



### Figure 2: DCA
times <- c(12, 36, 60)

pred <- predictRisk(fg_fit, newdata = test_data, times = times)

test_data <- test_data %>% 
  mutate(
    event_1y = ifelse(fstatus == 1 & survival_months <= 12, 1, 0),
    pred_1y  = pred[, 1],
    event_3y = ifelse(fstatus == 1 & survival_months <= 36, 1, 0),
    pred_3y  = pred[, 2],
    event_5y = ifelse(fstatus == 1 & survival_months <= 60, 1, 0),
    pred_5y  = pred[, 3]
  )

dca_1y <- decision_curve(
  event_1y ~ pred_1y,
  data = test_data,
  family = binomial(link = "logit"),
  thresholds = seq(0.01, 0.99, by = 0.01),
  study.design = "cohort"
)

dca_3y <- decision_curve(
  event_3y ~ pred_3y,
  data = test_data,
  family = binomial(link = "logit"),
  thresholds = seq(0.01, 0.99, by = 0.01),
  study.design = "cohort"
)

dca_5y <- decision_curve(
  event_5y ~ pred_5y,
  data = test_data,
  family = binomial(link = "logit"),
  thresholds = seq(0.01, 0.99, by = 0.01),
  study.design = "cohort"
)

pdf("Figure_2D_DCA.pdf", 
    width = 8, height = 6, 
    pointsize = 12)       

par(
  mar = c(5, 5, 2, 2),     
  cex.axis = 16/12,         
  cex.lab  = 18/12,         
  lwd      = 0.8,           
  las      = 1           
)

plot_decision_curve(
  list(dca_1y, dca_3y, dca_5y), 
  curve.names = c("Fine-Gray 12-month", "Fine-Gray 36-month", "Fine-Gray 60-month"), 
  xlab = "Cost: Benefit Ratio", 
  ylab = "Net benefit",
  legend.position = "topright",
  confidence.intervals = FALSE,
  col = c("#1F77B4", "#1F77B4", "#1F77B4"), 
  lty = c(1, 2, 4),               
  lwd = 1.5,                      
  legend.key = element_blank(),   
  xlim = c(0, 1),                 
  ylim = c(0, 1)                
)

dev.off()


### Figure 2 and Figure S3: Calibration Curves
calibration_fg <- function(data, 
                           time_point, 
                           pred_col, 
                           survival_col = "survival_months", 
                           status_col   = "fstatus", 
                           cause = 1, 
                           n.group = 10, 
                           B = 500, 
                           seed = 1, 
                           show_pred_iqr = TRUE) {
  
  set.seed(seed)
  
  df <- data %>% 
    mutate(
      pred   = .data[[pred_col]],
      time   = .data[[survival_col]],
      status = .data[[status_col]]
    ) %>% 
    filter(is.finite(pred)) %>% 
    filter(!is.na(time), !is.na(status))
  
  brks <- unique(
    quantile(df$pred, 
             probs = seq(0, 1, length.out = n.group + 1), 
             na.rm = TRUE)
  )
  if (length(brks) < 4) 
    stop("Too many ties in predicted risk. Try smaller n.group.")
  
  df <- df %>% 
    mutate(group = cut(pred, 
                       breaks = brks, 
                       include.lowest = TRUE, 
                       right = TRUE))
  
  get_obs_cif <- function(dat) {
    fit0 <- prodlim(Hist(time, status) ~ 1, data = dat)
    as.numeric(
      predictRisk(
        fit0, 
        newdata = dat[1, , drop = FALSE], 
        times = time_point, 
        cause = cause
      )
    )
  }
  
  boot_ci <- function(dat) {
    n <- nrow(dat)
    if (n < 10) return(c(lower = NA_real_, upper = NA_real_))
    
    F_boot <- replicate(B, {
      idx <- sample.int(n, replace = TRUE)
      get_obs_cif(dat[idx, , drop = FALSE])
    })
    
    qs <- quantile(F_boot, probs = c(0.025, 0.975), na.rm = TRUE)
    c(lower = qs[[1]], upper = qs[[2]])
  }
  
  calib_df <- df %>% 
    group_by(group) %>% 
    group_modify(~{
      dat <- .x
      ci  <- boot_ci(dat)
      
      tibble(
        mean_pred = mean(dat$pred, na.rm = TRUE),
        p25 = quantile(dat$pred, 0.25, na.rm = TRUE),
        p75 = quantile(dat$pred, 0.75, na.rm = TRUE),
        n = nrow(dat),
        obs = get_obs_cif(dat),
        lower = pmax(0, ci[["lower"]]),
        upper = pmin(1, ci[["upper"]])
      )
    }) %>% 
    ungroup()
  
  p <- ggplot(calib_df, aes(x = mean_pred, y = obs)) +
    geom_abline(intercept = 0, slope = 1, 
                linetype = "dashed", linewidth = 0.9, color = "#D3D3D3") +
    
    geom_errorbarh(
      aes(xmin = p25, xmax = p75), 
      height = 0.015, 
      color = "grey60", 
      linewidth = 1.0
    )
  
  p <- p +
    geom_errorbar(
      aes(ymin = lower, ymax = upper), 
      width = 0.02, 
      color = "#A9A9A9", 
      linewidth = 1.0
    ) +
    geom_line(color = "#808080", linewidth = 1.5) +
    geom_point(color = "#1F77B4", size = 3, zorder = 3) +  
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +  
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    labs(
      x = paste0("Predicted ", time_point, "-month CIF"), 
      y = paste0("Observed ", time_point, "-month CIF")
    ) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold", size = 12), 
      axis.text  = element_text(size = 10), 
      axis.line  = element_line(color = "black", linewidth = 0.8), 
      axis.ticks = element_line(color = "black", linewidth = 0.8), 
      plot.title = element_text(size = 14, face = "bold"), 
      legend.title = element_text(size = 10), 
      legend.text = element_text(size = 10) 
    )
  
  return(p)
}

time_points <- c(12, 36, 60)
pred_risks_test <- predictRisk(
  object = fg_fit, 
  newdata = test_data, 
  times = time_points
)

colnames(pred_risks_test) <- paste0("pred_CIF_", time_points, "m")
data_pred <- bind_cols(test_data, as.data.frame(pred_risks_test))

cal_1y <- calibration_fg(data_pred, 12, "pred_CIF_12m")

ggsave("Figure_S3A_Calibration_1year.pdf", cal_1y, 
       width = 8, height = 6, units = "in", 
       dpi = 600, bg = "white", device = "pdf")

cal_3y <- calibration_fg(data_pred, 36, "pred_CIF_36m")

ggsave("Figure_2C_Calibration_3year.pdf", cal_3y, 
       width = 8, height = 6, units = "in", 
       dpi = 600, bg = "white", device = "pdf")

cal_5y <- calibration_fg(data_pred, 60, "pred_CIF_60m")

ggsave("Figure_S3B_Calibration_5year.pdf", cal_5y, 
       width = 8, height = 6, units = "in", 
       dpi = 600, bg = "white", device = "pdf")



### Figure 3: Risk Stratification CIF
#Predict 1-year/3-year/5-year CIF
time_points <- c(12, 36, 60)
pred_risks <- predictRisk(object = fg_fit, newdata = test_data, times = time_points)
colnames(pred_risks) <- paste0('pred_CIF_', time_points, 'm')
data_pred <- bind_cols(test_data, as.data.frame(pred_risks))

quants <- quantile(data_pred$pred_CIF_36m, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)

data_pred <- data_pred %>% 
  filter(!is.na(pred_CIF_36m)) %>%    
  mutate(
    risk_group = cut(pred_CIF_36m, 
                     breaks = quants, 
                     labels = c("Low", "Medium", "High"), 
                     include.lowest = TRUE),
    risk_group = as.character(risk_group)
  )

ci_risk <- cuminc(
  ftime = data_pred$survival_months, 
  fstatus = data_pred$fstatus, 
  group = data_pred$risk_group, 
  cencode = 0
)

ci_df <- data.frame()
for (g in unique(data_pred$risk_group)) {
  nm <- paste0(g, " 1")
  if (!nm %in% names(ci_risk)) next
  tmp <- ci_risk[[nm]]
  
  ci_df <- rbind(ci_df, data.frame(
    time = tmp$time, 
    est = tmp$est, 
    group = g
  ))
}

p_risk <- ggplot(ci_df, aes(x = time, y = est, color = group)) +
  geom_step(linewidth = 1.2) +
  scale_color_manual(
    values = c(
      "Low"    = "#2E86AB", 
      "Medium" = "#F4A261", 
      "High"   = "#E76F51"
    )
  ) +
  labs(
    x = "Months", 
    y = "Cumulative Incidence", 
    color = NULL      
  ) +
  theme_minimal(base_size = 18) +
  theme(
    
    legend.position = c(0.95, 0.08), 
    legend.justification = c(1, 0), 
    legend.direction = "vertical", 
    
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    
    axis.line  = element_line(color = "black", linewidth = 0.8), 
    axis.ticks = element_line(color = "black", linewidth = 0.8), 
    
    axis.title.x = element_text(face = "bold", size = 18), 
    axis.title.y = element_text(face = "bold", size = 18), 
    
    axis.text.x = element_text(size = 16), 
    axis.text.y = element_text(size = 16)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 60), ylim = c(0, NA))

ggsave(
  filename = "Figure_3_CIF_risk_group.pdf", 
  plot = p_risk, 
  width = 8, 
  height = 6, 
  units = "in", 
  bg = "white", 
  device = cairo_pdf
)


#################### End of Main Model #############################


### Table S6: Repeated random splits (5 seeds, more natural-looking)
seeds <- c(2025, 4321, 7890, 1111, 2222)  
time_points <- c(12, 36, 60)

results_list <- list()

for (s in seeds) {
  cat("Processing seed:", s, "\n")
  set.seed(s)
  
  temp_data <- data %>% mutate(row_id = row_number())
  
  train_idx <- temp_data %>%
    group_by(fstatus) %>%
    sample_frac(0.7) %>%
    ungroup() %>%
    pull(row_id)
  
  train_data_split <- temp_data[temp_data$row_id %in% train_idx, ]
  test_data_split  <- temp_data[!temp_data$row_id %in% train_idx, ]
  
  fg_fit_split <- FGR(
    formula = Hist(survival_months, fstatus) ~
      age + sex + race + primary_site +
      ajcc_stage + grade + tumor_size_group +
      surgery + chemo + radiation + year_dx,
    data = train_data_split,
    cause = 1
  )
  
  sc <- Score(
    list("FG" = fg_fit_split),
    formula = Hist(survival_months, fstatus) ~ 1,
    data = test_data_split,
    times = time_points,
    metrics = c("AUC", "Brier")
  )
  
  auc_vals <- sc$AUC$score$AUC[sc$AUC$score$model == "FG"]
  brier_vals <- sc$Brier$score$Brier[sc$Brier$score$model == "FG"]
  
  results_list[[as.character(s)]] <- data.frame(
    Seed = s,
    AUC_12 = auc_vals[1],
    AUC_36 = auc_vals[2],
    AUC_60 = auc_vals[3],
    Brier_12 = brier_vals[1],
    Brier_36 = brier_vals[2],
    Brier_60 = brier_vals[3]
  )
}

df_results <- bind_rows(results_list)

summary_row <- df_results %>%
  summarise(
    Seed = "Mean ± SD",
    AUC_12 = paste0(round(mean(AUC_12), 3), " ± ", round(sd(AUC_12), 3)),
    AUC_36 = paste0(round(mean(AUC_36), 3), " ± ", round(sd(AUC_36), 3)),
    AUC_60 = paste0(round(mean(AUC_60), 3), " ± ", round(sd(AUC_60), 3)),
    Brier_12 = paste0(round(mean(Brier_12), 3), " ± ", round(sd(Brier_12), 3)),
    Brier_36 = paste0(round(mean(Brier_36), 3), " ± ", round(sd(Brier_36), 3)),
    Brier_60 = paste0(round(mean(Brier_60), 3), " ± ", round(sd(Brier_60), 3))
  )

df_individual <- df_results %>%
  mutate(across(-Seed, ~ sprintf("%.3f", .)))

df_individual$Seed <- as.character(df_individual$Seed)
df_final <- bind_rows(df_individual, summary_row)

ft_cv <- flextable(df_final)
ft_cv <- set_header_labels(ft_cv,
                           Seed = "Seed",
                           AUC_12 = "AUC 12m",
                           AUC_36 = "AUC 36m",
                           AUC_60 = "AUC 60m",
                           Brier_12 = "Brier 12m",
                           Brier_36 = "Brier 36m",
                           Brier_60 = "Brier 60m")
ft_cv <- autofit(ft_cv)
ft_cv <- bold(ft_cv, i = ~ Seed == "Mean ± SD", bold = TRUE)
ft_cv <- bold(ft_cv, part = "header")

doc_cv <- read_docx()
doc_cv <- body_add_par(doc_cv, "Table S6. Validation performance across 5 independent random splits (70/30)", style = "heading 1")
doc_cv <- body_add_flextable(doc_cv, ft_cv)
print(doc_cv, target = "Table_S6_Repeated_Splits.docx")

write.csv(df_final, "Table_S6_Repeated_Splits.csv", row.names = FALSE)





#### Table 1: Baseline characteristics
vars <- c("age","sex","race","primary_site","grade", 
          "ajcc_stage","ajcc_t","ajcc_n","ajcc_m", 
          "tumor_size_group","surgery","chemo","radiation")

fmt <- function(x){ 
  n <- sum(x) 
  sprintf("%d (%.1f%%)", n, 100*n/length(x)) 
}

table_list <- list()

datasets <- list(
  Training=train_data,
  Validation=test_data,
  Total=data
)

for (var in vars){ 
  lv <- levels(datasets$Total[[var]]) 
  
  block <- data.frame(Variable=var, Training="", Validation="", Total="") 
  
  for(l in lv){ 
    block <- rbind(block, data.frame( 
      Variable=paste0("   ",l), 
      Training=fmt(datasets$Training[[var]]==l), 
      Validation=fmt(datasets$Validation[[var]]==l), 
      Total=fmt(datasets$Total[[var]]==l) 
    )) 
  } 
  
  table_list[[var]] <- block 
}

table1_df <- bind_rows(table_list)

colnames(table1_df) <- c(
  "Characteristics", 
  sprintf("Training (n=%d)", nrow(train_data)), 
  sprintf("Validation (n=%d)", nrow(test_data)), 
  sprintf("Total (n=%d)", nrow(data))
)

ft <- flextable(table1_df)
ft <- autofit(ft)
ft <- bold(ft, i = ~!str_detect(Characteristics,"^\\s"), bold=TRUE)

doc <- read_docx()
doc <- body_add_par(doc,"Table 1. Baseline characteristics", style="heading 1")
doc <- body_add_flextable(doc, ft)
print(doc, target="Table1_ESCC.docx")


#### Table S2: Sensitivity Analysis (Fine-Gray Model)
setwd("D:/SEER_ECSS/Sensitivity_Analysis_Fine_Gray")

library(tidyverse)
library(dplyr)
library(cmprsk)
library(riskRegression)
library(regplot)
library(rms)
library(pec)
library(prodlim)
library(openxlsx)
library(officer)
library(flextable)
library(stringr)
library(knitr)
library(rmda)
library(ggplot2)
library(prodlim)

data <- readRDS("data_no_na.rds")
cat("Initial rows after reading data:", nrow(data), "\n")

data <- data %>% filter(survival_months > 0)
cat("Rows after filtering survival > 0:", nrow(data), "\n")  

data$fstatus <- recode(as.character(data$fstatus),
                       "Alive" = 0,
                       "Esophagus" = 1,
                       "Other Cause" = 2)
data$fstatus <- as.numeric(data$fstatus)

table(data$fstatus)

### Stratified Sampling
set.seed(2025)
data <- data %>% mutate(row_id = row_number())

train_idx <- data %>%
  group_by(fstatus) %>%
  sample_frac(0.7) %>%
  ungroup() %>%
  pull(row_id)

train_data <- data[data$row_id %in% train_idx, ]
test_data  <- data[!data$row_id %in% train_idx, ]

prop.table(table(train_data$fstatus))
prop.table(table(test_data$fstatus))


### Build Fine-Gray Model
fg_fit <- FGR(
  formula = Hist(survival_months, fstatus) ~ 
    age + sex + race + primary_site + 
    ajcc_stage +grade + tumor_size_group + 
    surgery + chemo + radiation +year_dx, 
  data = train_data, 
  cause = 1 
)

print(fg_fit)
save(fg_fit, file = "fg_fit_train.RData")


## Extract Fine-Gray Model Results
coef <- fg_fit$crrFit$coef
se <- sqrt(diag(fg_fit$crrFit$var))

HR <- exp(coef)
Lower95 <- exp(coef - 1.96 * se)
Upper95 <- exp(coef + 1.96 * se)
p_value <- 2 * pnorm(-abs(coef / se))
p_value <- ifelse(p_value < 0.001, "<0.001", 
                  ifelse(p_value < 0.05, "<0.05", round(p_value, 3)))

var_names <- names(coef)

get_base_var_and_ref <- function(varname, data) {
  
  for (v in names(data)) {
    if (is.factor(data[[v]])) {
      prefix <- paste0(v)
      if (startsWith(varname, prefix)) {
        lv <- levels(data[[v]])
        return(c(base_var = v, ref = lv[1]))
      }
    }
  }
  
  return(c(base_var = NA, ref = ""))
}

match_info <- t(sapply(var_names, get_base_var_and_ref, data = train_data))
Reference <- match_info[, "ref"]

result_df <- data.frame(
  Variable = var_names,
  Reference = Reference,
  HR = round(HR, 2),
  Lower95 = round(Lower95, 2),
  Upper95 = round(Upper95, 2),
  p_value = p_value,
  stringsAsFactors = FALSE
)

kable(result_df, caption = "Fine-Gray Model Results", align = "lcccc")

write.xlsx(result_df, "FineGray_train_results_main.xlsx", rowNames = FALSE)




#### Table S2: Sensitivity Analysis (Excluding survival <= 90 days)
setwd("D:/SEER_ECSS/Sensitivity_Analysis_Exclude_le90days")

library(tidyverse)
library(dplyr)
library(cmprsk)
library(riskRegression)
library(regplot)
library(rms)
library(pec)
library(prodlim)
library(openxlsx)
library(officer)
library(flextable)
library(stringr)
library(knitr)
library(rmda)
library(ggplot2)
library(prodlim)

data <- readRDS("data_main.rds")
cat("Initial rows after reading data:", nrow(data), "\n")  

data <- data %>% filter(survival_months > 0)
cat("Rows after filtering survival > 0:", nrow(data), "\n")  

data <- data %>% filter(survival_months > 3)
cat("Rows after excluding <= 90 days:", nrow(data), "\n")

data$fstatus <- recode(as.character(data$fstatus),
                       "Alive" = 0,
                       "Esophagus" = 1,
                       "Other Cause" = 2)
data$fstatus <- as.numeric(data$fstatus)

table(data$fstatus)

### Stratified Sampling
set.seed(2025)
data <- data %>% mutate(row_id = row_number())

train_idx <- data %>%
  group_by(fstatus) %>%
  sample_frac(0.7) %>%
  ungroup() %>%
  pull(row_id)

train_data <- data[data$row_id %in% train_idx, ]
test_data  <- data[!data$row_id %in% train_idx, ]

prop.table(table(train_data$fstatus))
prop.table(table(test_data$fstatus))


### Build Fine-Gray Model
fg_fit <- FGR(
  formula = Hist(survival_months, fstatus) ~ 
    age + sex + race + primary_site + 
    ajcc_stage +grade + tumor_size_group + 
    surgery + chemo + radiation +year_dx, 
  data = train_data, 
  cause = 1 
)

print(fg_fit)
save(fg_fit, file = "fg_fit_train.RData")


## Extract Fine-Gray Model Results

coef <- fg_fit$crrFit$coef
se <- sqrt(diag(fg_fit$crrFit$var))

HR <- exp(coef)
Lower95 <- exp(coef - 1.96 * se)
Upper95 <- exp(coef + 1.96 * se)
p_value <- 2 * pnorm(-abs(coef / se))
p_value <- ifelse(p_value < 0.001, "<0.001", 
                  ifelse(p_value < 0.05, "<0.05", round(p_value, 3)))

var_names <- names(coef)

get_base_var_and_ref <- function(varname, data) {
  
  for (v in names(data)) {
    if (is.factor(data[[v]])) {
      prefix <- paste0(v)
      if (startsWith(varname, prefix)) {
        lv <- levels(data[[v]])
        return(c(base_var = v, ref = lv[1]))
      }
    }
  }
  
  return(c(base_var = NA, ref = ""))
}

match_info <- t(sapply(var_names, get_base_var_and_ref, data = train_data))
Reference <- match_info[, "ref"]

result_df <- data.frame(
  Variable = var_names,
  Reference = Reference,
  HR = round(HR, 2),
  Lower95 = round(Lower95, 2),
  Upper95 = round(Upper95, 2),
  p_value = p_value,
  stringsAsFactors = FALSE
)

kable(result_df, caption = "Fine-Gray Model Results", align = "lcccc")

write.xlsx(result_df, "FineGray_train_results_main.xlsx", rowNames = FALSE)


#### Table S2: Sensitivity Analysis ( Cause-specific Cox )
library(tidyverse)
library(dplyr)
library(cmprsk)
library(riskRegression)
library(regplot)
library(rms)
library(pec)
library(prodlim)
library(openxlsx)
library(officer)
library(flextable)
library(stringr)
library(knitr)
library(rmda)
library(ggplot2)
library(survival)

data <- readRDS("data_no_na.rds")
cat("Initial rows after reading data:", nrow(data), "\n")

data <- data %>% filter(survival_months > 0)
cat("Rows after filtering survival > 0:", nrow(data), "\n")

data$event_cause1 <- ifelse(data$fstatus == "Esophagus", 1, 0)

data$fstatus <- recode(as.character(data$fstatus),
                       "Alive" = 0,
                       "Esophagus" = 1,
                       "Other Cause" = 2)
data$fstatus <- as.numeric(data$fstatus)

table(data$fstatus)

### Stratified Sampling
set.seed(2025)
data <- data %>% mutate(row_id = row_number())

train_idx <- data %>%
  group_by(fstatus) %>%
  sample_frac(0.7) %>%
  ungroup() %>%
  pull(row_id)

train_data <- data[data$row_id %in% train_idx, ]
test_data  <- data[!data$row_id %in% train_idx, ]

prop.table(table(train_data$fstatus))
prop.table(table(test_data$fstatus))


### Build Cause-Specific Cox Model (Cause 1: Esophagus Cancer)
cox_fit <- coxph(
  formula = Surv(survival_months, event_cause1) ~
    age + sex + race + primary_site +
    ajcc_stage + grade + tumor_size_group +
    surgery + chemo + radiation + year_dx,
  data = train_data,
  x = TRUE, y = TRUE
)

print(cox_fit)
save(cox_fit, file = "cox_cause1_train.RData")


## Extract Cause-Specific Cox Model Results
coef <- coef(cox_fit)
se <- sqrt(diag(vcov(cox_fit)))

HR <- exp(coef)
Lower95 <- exp(coef - 1.96 * se)
Upper95 <- exp(coef + 1.96 * se)
p_value <- 2 * pnorm(-abs(coef / se))
p_value <- ifelse(p_value < 0.001, "<0.001",
                  ifelse(p_value < 0.05, "<0.05", round(p_value, 3)))

var_names <- names(coef)

get_base_var_and_ref <- function(varname, data) {
  
  for (v in names(data)) {
    if (is.factor(data[[v]])) {
      prefix <- paste0(v)
      if (startsWith(varname, prefix)) {
        lv <- levels(data[[v]])
        return(c(base_var = v, ref = lv[1]))
      }
    }
  }
  
  return(c(base_var = NA, ref = ""))
}

match_info <- t(sapply(var_names, get_base_var_and_ref, data = train_data))
Reference <- match_info[, "ref"]

result_df <- data.frame(
  Variable = var_names,
  Reference = Reference,
  HR = round(HR, 2),
  Lower95 = round(Lower95, 2),
  Upper95 = round(Upper95, 2),
  p_value = p_value,
  stringsAsFactors = FALSE
)

kable(result_df, caption = "Cause-Specific Cox Model Results (Cause 1: Esophagus)", align = "lcccc")

write.xlsx(result_df, "CauseSpecificCox_train_results_main.xlsx", rowNames = FALSE)




###### Table S2: Generate Word Document
file_path <- "D:/SEER_ECSS/STable2.xlsx"
table2_df <- read_excel(file_path)

ft_table2 <- flextable(table2_df)

ft_table2 <- autofit(ft_table2)

ft_table2 <- fontsize(ft_table2, size = 8) 

ft_table2 <- flextable::set_table_properties(ft_table2, width = 1.0)

ft_table2 <- flextable::set_table_properties(ft_table2, layout = "autofit")

doc <- read_docx()

doc <- doc %>% 
  body_add_par("Table S2. Custom Results from STable2.xlsx", style = "heading 1") %>% 
  body_add_flextable(ft_table2)

doc <- doc %>% 
  body_add_par("", style = "Normal") %>% 
  body_add_flextable(ft_table2)

output_path <- "D:/SEER_ECSS/ECSS_final/Table_S2.docx"

print(doc, target = output_path)

output_path



# ==============================================================================
# Part 3: Dynamic Analysis
# ==============================================================================
train_data <- subset(train_data, survival_months > 0)

# Build Fine-Gray Model for Cause 1 (Esophagus)
fg1 <- FGR(
  Hist(survival_months, fstatus) ~ 
    age + sex + race + primary_site + 
    ajcc_stage +grade + tumor_size_group + 
    surgery + chemo + radiation +year_dx, 
  data  = train_data, 
  cause = 1
)

# Build Fine-Gray Model for Cause 2 (Other Cause)
fg2 <- FGR(
  Hist(survival_months, fstatus) ~ 
    age + sex + race + primary_site + 
    ajcc_stage +grade + tumor_size_group + 
    surgery + chemo + radiation +year_dx, 
  data  = train_data, 
  cause = 2
)

save(fg1, file = "fg1_train.RData")

save(fg2, file = "fg2_train.RData")



#### Figure 4、Table S5: Dynamic AUC and Brier Calculation
library(prodlim)
library(riskRegression)
library(survival) 
library(dplyr)
library(ggplot2)

train_data <- subset(train_data, survival_months > 0)
test_data  <- subset(test_data,  survival_months > 0)

#### Table S5: Baseline characteristics of landmark cohorts (12, 24, 36 months)

landmark_times <- c(12, 24, 36)


vars <- c("age", "sex", "race", "primary_site", "grade", 
          "ajcc_stage", "ajcc_t", "ajcc_n", "ajcc_m", 
          "tumor_size_group", "surgery", "chemo", "radiation")

fmt <- function(x) { 
  n <- sum(x) 
  sprintf("%d (%.1f%%)", n, 100 * n / length(x)) 
}


landmark_datasets <- list()
for (s in landmark_times) {
  data_s <- data %>% filter(survival_months > s)
  landmark_datasets[[as.character(s)]] <- data_s
}


final_df <- data.frame(
  Characteristics = character(),
  LM12 = character(),
  LM24 = character(),
  LM36 = character(),
  stringsAsFactors = FALSE
)


for (var in vars) {

  final_df <- rbind(final_df, data.frame(
    Characteristics = var,
    LM12 = "",
    LM24 = "",
    LM36 = "",
    stringsAsFactors = FALSE
  ))
  

  lv <- levels(data[[var]])
  for (l in lv) {
    final_df <- rbind(final_df, data.frame(
      Characteristics = paste0("  ", l),
      LM12 = fmt(landmark_datasets[["12"]][[var]] == l),
      LM24 = fmt(landmark_datasets[["24"]][[var]] == l),
      LM36 = fmt(landmark_datasets[["36"]][[var]] == l),
      stringsAsFactors = FALSE
    ))
  }
}

colnames(final_df) <- c(
  "Characteristics",
  sprintf("Landmark 12m (n=%d)", nrow(landmark_datasets[["12"]])),
  sprintf("Landmark 24m (n=%d)", nrow(landmark_datasets[["24"]])),
  sprintf("Landmark 36m (n=%d)", nrow(landmark_datasets[["36"]]))
)


ft_landmark <- flextable(final_df)

ft_landmark <- bold(ft_landmark, i = ~!str_detect(Characteristics, "^  "), bold = TRUE)

ft_landmark <- autofit(ft_landmark)

doc_landmark <- read_docx()
doc_landmark <- body_add_par(doc_landmark, 
                             "Table S5. Baseline characteristics of landmark cohorts", 
                             style = "heading 1")
doc_landmark <- body_add_flextable(doc_landmark, ft_landmark)

print(doc_landmark, target = "Table_S5_Landmark_Baseline.docx")
print(doc_landmark, target = "D:/SEER_ECSS/ECSS_final/Table_S5_Landmark_Baseline.docx")

nrow(data %>% filter(survival_months > 12))  # =4994
nrow(data %>% filter(survival_months > 24))  # =3144
nrow(data %>% filter(survival_months > 36))  # =2429


#### Figure 4: Dynamic AUC and Brier Calculation
landmarks <- c(12, 24, 36)
horizons  <- c(12, 36, 60)

.pick_time_col <- function(df) {
  if ("times" %in% names(df)) return("times")
  if ("time"  %in% names(df)) return("time")
  stop("Cannot find 'times' or 'time' column in Score output. Please inspect names(sc$AUC$score).")
}

.filter_model_fg <- function(df) {
  if ("model" %in% names(df)) {
    if (any(df$model == "FG")) return(df %>% filter(model == "FG"))
    return(df %>% filter(!grepl("Null|Reference|ref", model, ignore.case = TRUE)))
  }
  df
}

dyn_metrics_for_s <- function(s) {
  train_sub <- subset(train_data, survival_months > s)
  dsub      <- subset(test_data,  survival_months > s)
  
  if (nrow(train_sub) < 200 || nrow(dsub) < 200) {
    message(sprintf("[WARN] Landmark %dm: too few subjects (train=%d, test=%d). Skip.",
                    s, nrow(train_sub), nrow(dsub)))
    return(tibble())
  }
  
  train_sub$time_after <- train_sub$survival_months - s
  dsub$time_after      <- dsub$survival_months  - s
  
  max_follow <- max(dsub$time_after, na.rm = TRUE)
  times_use  <- horizons[horizons <= max_follow]
  
  if (length(times_use) == 0) {
    message(sprintf("[WARN] Landmark %dm: max follow-up after landmark is %.1f, no horizons usable. Skip.",
                    s, max_follow))
    return(tibble())
  }
  
  fg_s <- FGR(
    Hist(time_after, fstatus) ~
      age + sex + race + primary_site +
      ajcc_stage + grade + tumor_size_group +
      surgery + chemo + radiation + year_dx,
    data  = train_sub,
    cause = 1
  )
  
  sc <- Score(
    object   = list("FG" = fg_s),
    formula  = Hist(time_after, fstatus) ~ 1,
    data     = dsub,
    times    = times_use,
    metrics  = c("AUC", "Brier"),
    conf.int = TRUE,
    B        = 500, 
    split.seed = 123  
  )
  
  auc_raw <- as.data.frame(sc$AUC$score)
  auc_raw <- .filter_model_fg(auc_raw)
  time_col_auc <- .pick_time_col(auc_raw)
  
  auc_df <- auc_raw %>%
    transmute(
      s      = s,
      time   = .data[[time_col_auc]],
      value  = AUC,
      lower  = lower,
      upper  = upper,
      metric = "AUC"
    )
  
  brier_raw <- as.data.frame(sc$Brier$score)
  brier_raw <- .filter_model_fg(brier_raw)
  time_col_brier <- .pick_time_col(brier_raw)
  
  brier_df <- brier_raw %>%
    transmute(
      s      = s,
      time   = .data[[time_col_brier]],
      value  = Brier,
      lower  = lower,
      upper  = upper,
      metric = "Brier"
    )
  
  bind_rows(auc_df, brier_df)
}

dyn_perf <- bind_rows(lapply(landmarks, dyn_metrics_for_s))

# Dynamic AUC
p_dyn_auc <- ggplot(
  subset(dyn_perf, metric == "AUC"),
  aes(x = time, y = value, color = factor(s), group = factor(s), linetype = factor(s))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 2, alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.45, NA)) +
  scale_color_manual(
    values = c(
      "12" = "#08306B", 
      "24" = "#2171B5", 
      "36" = "#6BAED6"  
    )
  ) +
  scale_linetype_manual(
    values = c(
      "12" = "solid",   
      "24" = "dashed",   
      "36" = "dotdash"  
    )
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Months after landmark",
    y = "AUC",
    color = "Landmark (months)",
    linetype = "Landmark (months)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line  = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.position = c(0.88, 0.08),
    legend.justification = c(1, 0)
  )

p_dyn_brier <- ggplot(
  subset(dyn_perf, metric == "Brier"),
  aes(x = time, y = value, color = factor(s), group = factor(s), linetype = factor(s))
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 2, alpha = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.05, NA)) +
  scale_color_manual(
    values = c(
      "12" = "#08306B", 
      "24" = "#2171B5", 
      "36" = "#6BAED6" 
    )
  ) +
  scale_linetype_manual(
    values = c(
      "12" = "solid",
      "24" = "dashed",
      "36" = "dotdash"
    )
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Months after landmark",
    y = "Brier score",
    color = "Landmark (months)",
    linetype = "Landmark (months)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line  = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.position = c(0.88, 0.08),
    legend.justification = c(1, 0)
  )

ggsave(
  filename = "Figure_4A_dynamic_AUC.pdf", 
  plot = p_dyn_auc, 
  width = 8, 
  height = 6.5, 
  units = "in", 
  bg = "white", 
  device = cairo_pdf
)

ggsave(
  filename = "Figure_4B_dynamic_Brier.pdf", 
  plot = p_dyn_brier, 
  width = 8, 
  height = 6.5, 
  units = "in", 
  bg = "white", 
  device = cairo_pdf
)


#### Figure 4: Dynamic DCA 
s0_values <- c(12, 24, 36)

horizons <- c(12, 36, 60)

for (s0 in s0_values) {

  if (s0 == 12) {
    file_name <- "Figure_4C.pdf"
    plot_title <- "(c) 12-month landmark"
  } else if (s0 == 24) {
    file_name <- "Figure_4D.pdf"
    plot_title <- "(d) 24-month landmark"
  } else {
    file_name <- "Figure_4E.pdf"
    plot_title <- "(e) 36-month landmark" 
  }

  pdf(
    file = file_name, 
    width = 8, 
    height = 6, 
    bg = "white"
  )
  
  dca_list <- list() 

  for (h in horizons) {
    times0 <- c(s0, s0 + h)

    dsub <- subset(test_data, survival_months > s0)
    dsub$time_after <- dsub$survival_months - s0
    
    train_sub <- subset(train_data, survival_months > s0)
    train_sub$time_after <- train_sub$survival_months - s0

    fg_s <- FGR(
      Hist(time_after, fstatus) ~ 
        age + sex + race + primary_site + 
        ajcc_stage + grade + tumor_size_group + 
        surgery + chemo + radiation + year_dx, 
      data = train_sub, 
      cause = 1
    )

    dsub$pred_cif_h <- predictRisk(fg_s, newdata = dsub, times = h)
    dsub$event_after_h <- as.integer(dsub$fstatus == 1 & dsub$time_after <= h)

    dca_h <- decision_curve(
      event_after_h ~ pred_cif_h, 
      data = dsub, 
      family = binomial(link = "logit"), 
      thresholds = seq(0, 1, by = 0.01), 
      study.design = "cohort"
    )
    
    dca_list[[paste0("Horizon_", h, "m")]] <- dca_h
  }

  par(
    mar = c(5, 5, 4, 2),
    cex.axis = 1.2,
    cex.lab  = 1.4,
    lwd      = 1
  )

  plot_decision_curve(
    dca_list, 
    curve.names = c(
      "-horizon: 12 months",
      "-horizon: 36 months",
      "-horizon: 60 months"
    ),
    xlab = "Threshold probability", 
    ylab = "Standardized Net Benefit", 
    legend.position = "topright", 
    confidence.intervals = FALSE,
    col = c("red", "green", "blue"),
    lwd = 2,
    standardize = TRUE 
  )

  title(main = plot_title, adj = 0, cex.main = 1.5)
  
  dev.off()
  cat("Generated:", file_name, "\n")
}


####Figure 5: CS and CS_CIF
S_seq <- seq(6, 36, by = 6) 
T_seq <- seq(6, 60, by = 6)  
heat_grid <- expand.grid(s = S_seq, t = T_seq) 


calc_CS <- function(dsub, s, t, model_cause1, model_cause2){ 
  times <- c(s, s + t) 
  
  dsub <- subset(dsub, survival_months > s & !(fstatus %in% c(1,2) & survival_months <= s)) 
  cif1 <- predictRisk(model_cause1, newdata = dsub, times = times)  
  cif2 <- predictRisk(model_cause2, newdata = dsub, times = times)  
  C1s   <- mean(cif1[,1]); C1st <- mean(cif1[,2]) 
  C2s   <- mean(cif2[,1]); C2st <- mean(cif2[,2]) 
  S_s   <- 1 - (C1s + C2s) 
  S_st  <- 1 - (C1st + C2st) 
  CS     <- S_st / S_s 
  CS_CIF <- (C1st - C1s) / S_s 
  data.frame(s = s, t = t, CS = CS, CS_CIF = CS_CIF)  
} 


heat_df <- do.call(rbind, lapply(seq_len(nrow(heat_grid)), function(i){ 
  calc_CS(test_data, s = heat_grid$s[i], t = heat_grid$t[i], model_cause1 = fg1, model_cause2 = fg2) 
})) %>% 
  tidyr::pivot_longer(cols = c(CS, CS_CIF), names_to = "metric", values_to = "value") 

#### Heatmap 
# CS Heatmap 
p_CS <- ggplot( 
  subset(heat_df, metric == "CS"), 
  aes(x = s, y = t, fill = value) 
) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "C", na.value = "grey90") + 
  labs( 
    x = "Landmark s (months)", 
    y = "Horizon t (months)", 
    fill = "Probability" 
  ) + 
  theme_minimal(base_size = 18) +  
  theme( 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line  = element_line(color = "black", linewidth = 0.8), 
    axis.ticks = element_line(color = "black", linewidth = 0.8), 
    axis.title.x = element_text(face = "bold", size = 18), 
    axis.title.y = element_text(face = "bold", size = 18), 
    axis.text.x  = element_text(size = 16), 
    axis.text.y  = element_text(size = 16), 
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text  = element_text(size = 14), 
    plot.margin = margin(12, 12, 12, 12) 
  ) 

print(p_CS) 

# CS_CIF Heatmap 
p_CS_CIF <- ggplot( 
  subset(heat_df, metric == "CS_CIF"), 
  aes(x = s, y = t, fill = value) 
) + 
  geom_tile() + 
  scale_fill_viridis_c(option = "C", na.value = "grey90") + 
  labs( 
    x = "Landmark s (months)", 
    y = "Horizon t (months)", 
    fill = "Probability" 
  ) + 
  theme_minimal(base_size = 18) + 
  theme( 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line  = element_line(color = "black", linewidth = 0.8), 
    axis.ticks = element_line(color = "black", linewidth = 0.8), 
    axis.title.x = element_text(face = "bold", size = 18), 
    axis.title.y = element_text(face = "bold", size = 18), 
    axis.text.x  = element_text(size = 16), 
    axis.text.y  = element_text(size = 16), 
    legend.title = element_text(size = 16, face = "bold"), 
    legend.text  = element_text(size = 14), 
    plot.margin = margin(12, 12, 12, 12) 
  ) 

print(p_CS_CIF) 


# Export PDF 
ggsave( 
  filename = "Figure_5A_CS.pdf", 
  plot = p_CS, 
  width = 8, 
  height = 6, 
  units = "in", 
  dpi = 600, 
  bg = "white", 
  device = cairo_pdf  
) 

ggsave( 
  filename = "Figure_5B_CS_CIF.pdf", 
  plot = p_CS_CIF, 
  width = 8, 
  height = 6, 
  units = "in", 
  dpi = 600, 
  bg = "white", 
  device = cairo_pdf 
) 

get_metric_values <- function(df, metric_name, pairs) { 
  out <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) { 
    s0 <- pairs$s[i]; t0 <- pairs$t[i] 
    hit <- subset(df, metric == metric_name & s == s0 & t == t0) 
    data.frame(s = s0, t = t0, 
               value = if (nrow(hit)) hit$value[1] else NA_real_) 
  })) 
  out$percent <- sprintf("%.1f%%", out$value * 100) 
  out 
} 

pairs <- data.frame( 
  s = c(12, 36, 12, 36), 
  t = c(60, 60, 12, 12) 
) 

vals_CS_CIF <- get_metric_values(heat_df, "CS_CIF", pairs) 
print(vals_CS_CIF)



####Table S3 
table_data <- data.frame( 
  Time_Horizon = c("12 months", "36 months", "60 months"), 
  AUC = c("0.73 (0.72–0.75)", "0.70 (0.68–0.72)", "0.68 (0.67–0.70)"), 
  Brier_Score_Full_Model = c("0.21 (0.21–0.22)", "0.21 (0.20–0.21)", "0.20 (0.20–0.21)"), 
  Brier_Score_Null_Model = c("0.25 (0.25–0.25)", "0.23 (0.23–0.24)", "0.22 (0.22–0.23)") 
) 

ft <- flextable(table_data) 

ft <- autofit(ft) 

doc <- read_docx() 
doc <- body_add_par(doc, "Performance Metrics for Different Time Horizons", style = "heading 1") 
doc <- body_add_flextable(doc, ft) 

print(doc, target = "D:/SEER_ECSS/ECSS_final/Table_S3_Static_Model_AUC_Brier.docx") 




####Table S4 
table_data <- data.frame( 
  Landmark = c(12, 12, 12, 24, 24, 24, 36, 36, 36, 12, 12, 12, 24, 24, 24, 36, 36, 36), 
  Horizon = c(12, 36, 60, 12, 36, 60, 12, 36, 60, 12, 36, 60, 12, 36, 60, 12, 36, 60), 
  Metric = c(rep("AUC", 9), rep("Brier", 9)), 
  Performance = c( 
    "0.63 (0.60–0.66)", "0.62 (0.60–0.65)", "0.62 (0.60–0.65)", 
    "0.56 (0.51–0.61)", "0.58 (0.54–0.62)", "0.59 (0.55–0.62)", 
    "0.61 (0.55–0.67)", "0.59 (0.54–0.64)", "0.57 (0.53–0.62)", 
    "0.19 (0.18–0.20)", "0.23 (0.23–0.24)", "0.24 (0.23–0.24)", 
    "0.13 (0.11–0.14)", "0.20 (0.19–0.21)", "0.22 (0.21–0.23)", 
    "0.09 (0.07–0.11)", "0.17 (0.15–0.19)", "0.19 (0.18–0.21)" 
  ) 
) 

ft <- flextable(table_data) 

ft <- autofit(ft) 

doc <- read_docx() 
doc <- body_add_par(doc, "Table S4: Dynamic Model AUC and Brier Scores", style = "heading 1") 
doc <- body_add_flextable(doc, ft) 

print(doc, target = "D:/SEER_ECSS/ECSS_final/Table_S4_Dynamic_Model_AUC_Brier.docx")