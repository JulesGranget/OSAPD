# Required libraries install.packages("lme4")

library(readxl)
library(Matrix)
library(lme4)
library(lmerTest)  # For p-values
library(emmeans)
library(ggplot2)
library(e1071)  # For skewness and kurtosis
library(pbkrtest)
library(sjPlot)
library(broom.mixed)
library(writexl)
library(dplyr)

library(patchwork)





################################
######## Pxx COND VCCP ########
################################



fit_model_with_fallback <- function(df, tol_sing = 1e-4, maxfun = 2e5) {
  
  forms <- list(
    lmer_sujet_chan = Pxx ~ A * R * C + (1 | chan),                
    lmer_sujet  = Pxx ~ A * R * C + (1 | chan),                     
    lm_fixed    = Pxx ~ A * R * C                                    
  )
  
  ctrl <- lmerControl(
    optimizer   = "bobyqa",
    optCtrl     = list(maxfun = maxfun),
    calc.derivs = FALSE
  )
  
  bad_fit <- function(mod) {
    # singular OR lme4 convergence messages (often in optinfo)
    isSingular(mod, tol = tol_sing) ||
      !is.null(mod@optinfo$conv$lme4$messages) ||
      !is.null(mod@optinfo$conv$messages)
  }
  
  get_conv_msgs <- function(mod) {
    msgs <- c()
    if (!is.null(mod@optinfo$conv$lme4$messages)) msgs <- c(msgs, mod@optinfo$conv$lme4$messages)
    if (!is.null(mod@optinfo$conv$messages))      msgs <- c(msgs, mod@optinfo$conv$messages)
    if (length(msgs) == 0) return(NA_character_)
    paste(unique(msgs), collapse = " | ")
  }
  
  for (nm in names(forms)) {
    f <- forms[[nm]]
    
    if (nm == "lm_fixed") {
      mod <- lm(f, data = df)
      attr(mod, "fit_stage")   <- nm
      attr(mod, "conv_msgs")   <- NA_character_
      attr(mod, "is_singular") <- NA
      return(mod)
    }
    
    mod <- suppressWarnings(
      suppressMessages(
        try(lmer(f, data = df, control = ctrl), silent = TRUE)
      )
    )
    
    if (inherits(mod, "try-error")) next
    
    if (!bad_fit(mod)) {
      attr(mod, "fit_stage")   <- nm
      attr(mod, "conv_msgs")   <- get_conv_msgs(mod)
      attr(mod, "is_singular") <- isSingular(mod, tol = tol_sing)
      return(mod)
    }
  }
  
  # absolute last-resort safeguard
  mod <- lm(forms$lm_fixed, data = df)
  attr(mod, "fit_stage")   <- "lm_fixed"
  attr(mod, "conv_msgs")   <- NA_character_
  attr(mod, "is_singular") <- NA
  mod
}


root = "/home/jules/Documents/VCCP_JULES/Analyses/precompute/EXPORT_DF/patient_wise/Pxx"
outputdir_fig = "/home/jules/Documents/VCCP_JULES/Analyses/results/Pxx/LMM/diagnosis"
outputdir_df_lmm = "/home/jules/Documents/VCCP_JULES/Analyses/precompute/EXPORT_DF/patient_wise/Pxx/LMM_res/COND"

sujet = "LH018"
protocol_type = "VCCP"

band_list <- c("theta", "alpha", "beta", "gamma")
phase_list <- c("inspi", "expi")

phase_cycle_sel = phase_list[1]
band = band_list[1]

for (band in band_list) {
  
  filename <- paste0("/", sujet, "_df_Pxx_", band, "_VCCP_R.xlsx")
  df_raw <- read_excel(paste(root, filename, sep = "/"))
  
  df_raw$chan <- paste0(df_raw$sujet, "_", df_raw$chan)
  
  ROI_list <- unique(df_raw$ROI)
  
  for (phase_cycle_sel in phase_list) {
    
    df_phase <- subset(df_raw, phase_cycle == phase_cycle_sel)
    
    ROI_sel = ROI_list[17]
    
    for (ROI_sel in ROI_list) {
      
      tryCatch({
        
        print(ROI_sel)
        
        filename_export_diagnostic <- paste(protocol_type, band, phase_cycle_sel, ROI_sel, sep = "_")
        
        df_oneROI <- subset(df_phase, ROI == ROI_sel)
        
        # Optional: counts per subject (quick QA)
        df_count <- df_oneROI %>%
          group_by(sujet) %>%
          summarise(
            n_chan  = n_distinct(chan),
            n_cycle = n_distinct(cycle),
            .groups = "drop"
          )
        print(df_count)
        
        # Keep only required columns
        df <- df_oneROI[, c("sujet", "chan", "cycle", "Pxx", "A", "R", "C")]
        
        # Factors / reference levels
        df$sujet <- factor(df$sujet)
        df$chan  <- factor(df$chan)
        df$A  <- factor(df$A)
        df$R <- factor(df$R)
        df$C <- factor(df$C)
        
        df$A  <- relevel(df$A,  ref = "o")
        df$R  <- relevel(df$R,  ref = "o")
        df$C  <- relevel(df$C,  ref = "o")
        
        # ---- MODEL
        model <- fit_model_with_fallback(df)
        
        fit_stage <- attr(model, "fit_stage")
        conv_msgs <- attr(model, "conv_msgs")
        is_sing   <- attr(model, "is_singular")
        
        message("✅ Model fit stage: ", fit_stage,
                if (!is.na(is_sing)) paste0(" | singular=", is_sing) else "",
                if (!is.na(conv_msgs)) paste0(" | conv=", conv_msgs) else "")
        
        # ---- stats for histogram subtitle
        skew_chan <- round(skewness(df$Pxx), 2)
        kurt_chan <- round(kurtosis(df$Pxx), 2)
        
        # ------------------------------------------------------------
        # PLOT 1: subject-wise boxplot (ggplot)
        # ------------------------------------------------------------
        p_box <- ggplot(df, aes(x = sujet, y = Pxx, color = sujet, fill = sujet)) +
          geom_boxplot(width = .2, alpha = .5, outlier.alpha = 0,
                       position = position_dodge(.9)) +
          stat_summary(fun = median, geom = "point", size = 2,
                       position = position_dodge(.9), color = "white") +
          labs(
            title = paste(filename_export_diagnostic, "Pxx — Subject-wise", sep = " | ")
          ) +
          theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none"
          )
        
        # ------------------------------------------------------------
        # PLOT 2: histogram (ggplot)
        # ------------------------------------------------------------
        p_hist <- ggplot(df, aes(x = Pxx)) +
          geom_histogram(bins = 30, fill = "lightblue", color = "white") +
          labs(
            title    = "Histogram",
            subtitle = paste("kurtosis:", kurt_chan, "| skewness:", skew_chan),
            x        = "Pxx values",
            y        = "Count"
          ) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
        
        # ------------------------------------------------------------
        # PLOT 3: QQ plot of residuals (ggplot)
        # ------------------------------------------------------------
        res_df <- data.frame(res = resid(model))
        
        p_qq <- ggplot(res_df, aes(sample = res)) +
          stat_qq() +
          stat_qq_line() +
          labs(
            title = "QQ plot (residuals)",
            subtitle = paste0("fit_stage=", fit_stage,
                              if (!is.na(is_sing)) paste0(" | singular=", is_sing) else "")
          ) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
        
        # ------------------------------------------------------------
        # PLOT 4: Residuals vs Fitted (homoskedasticity check)
        # Equivalent to: plot(fitted(m7), resid(m7))
        # ------------------------------------------------------------
        
        diag_df <- data.frame(
          fitted = as.numeric(fitted(model)),
          resid  = as.numeric(resid(model))
        )
        
        p_homo <- ggplot(diag_df, aes(x = fitted, y = resid)) +
          geom_point(alpha = 0.35, size = 1) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          geom_smooth(method = "loess", se = FALSE, color = "black") +
          labs(
            title = "Residuals vs Fitted",
            subtitle = "Homoskedasticity check (look for constant spread, no funnel)",
            x = "Fitted values",
            y = "Residuals"
          ) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
        
        # ------------------------------------------------------------
        # COMBINE into ONE figure (patchwork)
        # Layout: boxplot on top, hist + qq below
        # ------------------------------------------------------------
        p_all <- (p_box + p_hist) / (p_qq + p_homo) +
          plot_annotation(
            title = paste("Diagnostics:", filename_export_diagnostic),
            theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
          )
        
        print(p_all)
        # ------------------------------------------------------------
        # SAVE one single PNG
        # ------------------------------------------------------------
        file_diag <- paste("DIAGNOSTIC", filename_export_diagnostic, "Pxx.png", sep = "_")
        
        ggsave(
          filename = paste(outputdir_fig, file_diag, sep = "/"),
          plot     = p_all,
          width    = 12,
          height   = 8,
          dpi      = 150
        )
        
        
        # ------------------------------------------------------------
        # EXPORT MODEL RESULTS
        # - tab_model for HTML-like report
        # - broom.mixed fixed effects + add fit diagnostics
        # ------------------------------------------------------------
        
        tab_model(model, show.re.var = TRUE, show.icc = TRUE, show.r2 = TRUE, show.se = TRUE)
        
        model_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
          mutate(
            fit_stage   = fit_stage,
            singular    = if (!is.na(is_sing)) is_sing else NA,
            conv_msgs   = if (!is.na(conv_msgs)) conv_msgs else NA,
            band        = band,
            phase_cycle = phase_cycle_sel,
            ROI         = ROI_sel
          )
        
        filesxlsx_ROI <- paste("RES", filename_export_diagnostic, "Pxx.xlsx", sep = "_")
        writexl::write_xlsx(model_df, paste(outputdir_df_lmm, filesxlsx_ROI, sep = "/"))
        
      }, error = function(e) {
        message("❌ Error for ", band, " / ", phase_cycle_sel, " / ", ROI_sel, " : ", conditionMessage(e))
      })
      
    }
  }
}










########################################
######## Pxx COND attention ########
########################################



fit_model_with_fallback <- function(df, tol_sing = 1e-4, maxfun = 2e5) {
  
  forms <- list(
    lmer_sujet_chan = Pxx ~ cond + (1 | chan),                
    lmer_sujet  = Pxx ~ cond + (1 | chan),                     
    lm_fixed    = Pxx ~ cond                                    
  )
  
  ctrl <- lmerControl(
    optimizer   = "bobyqa",
    optCtrl     = list(maxfun = maxfun),
    calc.derivs = FALSE
  )
  
  bad_fit <- function(mod) {
    # singular OR lme4 convergence messages (often in optinfo)
    isSingular(mod, tol = tol_sing) ||
      !is.null(mod@optinfo$conv$lme4$messages) ||
      !is.null(mod@optinfo$conv$messages)
  }
  
  get_conv_msgs <- function(mod) {
    msgs <- c()
    if (!is.null(mod@optinfo$conv$lme4$messages)) msgs <- c(msgs, mod@optinfo$conv$lme4$messages)
    if (!is.null(mod@optinfo$conv$messages))      msgs <- c(msgs, mod@optinfo$conv$messages)
    if (length(msgs) == 0) return(NA_character_)
    paste(unique(msgs), collapse = " | ")
  }
  
  for (nm in names(forms)) {
    f <- forms[[nm]]
    
    if (nm == "lm_fixed") {
      mod <- lm(f, data = df)
      attr(mod, "fit_stage")   <- nm
      attr(mod, "conv_msgs")   <- NA_character_
      attr(mod, "is_singular") <- NA
      return(mod)
    }
    
    mod <- suppressWarnings(
      suppressMessages(
        try(lmer(f, data = df, control = ctrl), silent = TRUE)
      )
    )
    
    if (inherits(mod, "try-error")) next
    
    if (!bad_fit(mod)) {
      attr(mod, "fit_stage")   <- nm
      attr(mod, "conv_msgs")   <- get_conv_msgs(mod)
      attr(mod, "is_singular") <- isSingular(mod, tol = tol_sing)
      return(mod)
    }
  }
  
  # absolute last-resort safeguard
  mod <- lm(forms$lm_fixed, data = df)
  attr(mod, "fit_stage")   <- "lm_fixed"
  attr(mod, "conv_msgs")   <- NA_character_
  attr(mod, "is_singular") <- NA
  mod
}


root = "/home/jules/Documents/VCCP_JULES/Analyses/precompute/EXPORT_DF/patient_wise/Pxx"
outputdir_fig = "/home/jules/Documents/VCCP_JULES/Analyses/results/Pxx/LMM/diagnosis"
outputdir_df_lmm = "/home/jules/Documents/VCCP_JULES/Analyses/precompute/EXPORT_DF/patient_wise/Pxx/LMM_res/COND"

sujet = "LH018"
protocol_type = "attention"

band_list <- c("theta", "alpha", "beta", "gamma")
phase_list <- c("inspi", "expi")

phase_cycle_sel = phase_list[1]
band = band_list[1]

for (band in band_list) {
  
  filename <- paste0("/", sujet, "_df_Pxx_", band, "_attention_R.xlsx")
  df_raw <- read_excel(paste(root, filename, sep = "/"))
  
  df_raw$chan <- paste0(df_raw$sujet, "_", df_raw$chan)
  
  ROI_list <- unique(df_raw$ROI)
  
  for (phase_cycle_sel in phase_list) {
    
    df_phase <- subset(df_raw, phase_cycle == phase_cycle_sel)
    
    ROI_sel = ROI_list[1]
    
    for (ROI_sel in ROI_list) {
      
      tryCatch({
        
        print(ROI_sel)
        
        filename_export_diagnostic <- paste(protocol_type, band, phase_cycle_sel, ROI_sel, sep = "_")
        
        df_oneROI <- subset(df_phase, ROI == ROI_sel)
        
        # Optional: counts per subject (quick QA)
        df_count <- df_oneROI %>%
          group_by(sujet) %>%
          summarise(
            n_chan  = n_distinct(chan),
            n_cycle = n_distinct(cycle),
            .groups = "drop"
          )
        print(df_count)
        
        # Keep only required columns
        df <- df_oneROI[, c("sujet", "chan", "cycle", "Pxx", "cond")]
        
        # Factors / reference levels
        df$sujet <- factor(df$sujet)
        df$chan  <- factor(df$chan)
        df$cond  <- factor(df$cond)

        df$cond  <- relevel(df$cond,  ref = "RB")

        # ---- MODEL
        model <- fit_model_with_fallback(df)
        
        fit_stage <- attr(model, "fit_stage")
        conv_msgs <- attr(model, "conv_msgs")
        is_sing   <- attr(model, "is_singular")
        
        message("✅ Model fit stage: ", fit_stage,
                if (!is.na(is_sing)) paste0(" | singular=", is_sing) else "",
                if (!is.na(conv_msgs)) paste0(" | conv=", conv_msgs) else "")
        
        # ---- stats for histogram subtitle
        skew_chan <- round(skewness(df$Pxx), 2)
        kurt_chan <- round(kurtosis(df$Pxx), 2)
        
        # ------------------------------------------------------------
        # PLOT 1: subject-wise boxplot (ggplot)
        # ------------------------------------------------------------
        p_box <- ggplot(df, aes(x = sujet, y = Pxx, color = sujet, fill = sujet)) +
          geom_boxplot(width = .2, alpha = .5, outlier.alpha = 0,
                       position = position_dodge(.9)) +
          stat_summary(fun = median, geom = "point", size = 2,
                       position = position_dodge(.9), color = "white") +
          labs(
            title = paste(filename_export_diagnostic, "Pxx — Subject-wise", sep = " | ")
          ) +
          theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none"
          )
        
        # ------------------------------------------------------------
        # PLOT 2: histogram (ggplot)
        # ------------------------------------------------------------
        p_hist <- ggplot(df, aes(x = Pxx)) +
          geom_histogram(bins = 30, fill = "lightblue", color = "white") +
          labs(
            title    = "Histogram",
            subtitle = paste("kurtosis:", kurt_chan, "| skewness:", skew_chan),
            x        = "Pxx values",
            y        = "Count"
          ) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
        
        # ------------------------------------------------------------
        # PLOT 3: QQ plot of residuals (ggplot)
        # ------------------------------------------------------------
        res_df <- data.frame(res = resid(model))
        
        p_qq <- ggplot(res_df, aes(sample = res)) +
          stat_qq() +
          stat_qq_line() +
          labs(
            title = "QQ plot (residuals)",
            subtitle = paste0("fit_stage=", fit_stage,
                              if (!is.na(is_sing)) paste0(" | singular=", is_sing) else "")
          ) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
        
        # ------------------------------------------------------------
        # PLOT 4: Residuals vs Fitted (homoskedasticity check)
        # Equivalent to: plot(fitted(m7), resid(m7))
        # ------------------------------------------------------------
        
        diag_df <- data.frame(
          fitted = as.numeric(fitted(model)),
          resid  = as.numeric(resid(model))
        )
        
        p_homo <- ggplot(diag_df, aes(x = fitted, y = resid)) +
          geom_point(alpha = 0.35, size = 1) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          geom_smooth(method = "loess", se = FALSE, color = "black") +
          labs(
            title = "Residuals vs Fitted",
            subtitle = "Homoskedasticity check (look for constant spread, no funnel)",
            x = "Fitted values",
            y = "Residuals"
          ) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
        
        # ------------------------------------------------------------
        # COMBINE into ONE figure (patchwork)
        # Layout: boxplot on top, hist + qq below
        # ------------------------------------------------------------
        p_all <- (p_box + p_hist) / (p_qq + p_homo) +
          plot_annotation(
            title = paste("Diagnostics:", filename_export_diagnostic),
            theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
          )
        
        print(p_all)
        # ------------------------------------------------------------
        # SAVE one single PNG
        # ------------------------------------------------------------
        file_diag <- paste("DIAGNOSTIC", filename_export_diagnostic, "Pxx.png", sep = "_")
        
        ggsave(
          filename = paste(outputdir_fig, file_diag, sep = "/"),
          plot     = p_all,
          width    = 12,
          height   = 8,
          dpi      = 150
        )
        
        
        # ------------------------------------------------------------
        # EXPORT MODEL RESULTS
        # - tab_model for HTML-like report
        # - broom.mixed fixed effects + add fit diagnostics
        # ------------------------------------------------------------
        
        tab_model(model, show.re.var = TRUE, show.icc = TRUE, show.r2 = TRUE, show.se = TRUE)
        
        model_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
          mutate(
            fit_stage   = fit_stage,
            singular    = if (!is.na(is_sing)) is_sing else NA,
            conv_msgs   = if (!is.na(conv_msgs)) conv_msgs else NA,
            band        = band,
            phase_cycle = phase_cycle_sel,
            ROI         = ROI_sel
          )
        
        filesxlsx_ROI <- paste("RES", filename_export_diagnostic, "Pxx.xlsx", sep = "_")
        writexl::write_xlsx(model_df, paste(outputdir_df_lmm, filesxlsx_ROI, sep = "/"))
        
      }, error = function(e) {
        message("❌ Error for ", band, " / ", phase_cycle_sel, " / ", ROI_sel, " : ", conditionMessage(e))
      })
      
    }
  }
}






################################
######## Pxx COND HV ########
################################



fit_model_with_fallback <- function(df, tol_sing = 1e-4, maxfun = 2e5) {
  
  forms <- list(
    lmer_sujet_chan = Pxx ~ cond + (1 | chan),                
    lmer_sujet  = Pxx ~ cond + (1 | chan),                     
    lm_fixed    = Pxx ~ cond                                    
  )
  
  ctrl <- lmerControl(
    optimizer   = "bobyqa",
    optCtrl     = list(maxfun = maxfun),
    calc.derivs = FALSE
  )
  
  bad_fit <- function(mod) {
    # singular OR lme4 convergence messages (often in optinfo)
    isSingular(mod, tol = tol_sing) ||
      !is.null(mod@optinfo$conv$lme4$messages) ||
      !is.null(mod@optinfo$conv$messages)
  }
  
  get_conv_msgs <- function(mod) {
    msgs <- c()
    if (!is.null(mod@optinfo$conv$lme4$messages)) msgs <- c(msgs, mod@optinfo$conv$lme4$messages)
    if (!is.null(mod@optinfo$conv$messages))      msgs <- c(msgs, mod@optinfo$conv$messages)
    if (length(msgs) == 0) return(NA_character_)
    paste(unique(msgs), collapse = " | ")
  }
  
  for (nm in names(forms)) {
    f <- forms[[nm]]
    
    if (nm == "lm_fixed") {
      mod <- lm(f, data = df)
      attr(mod, "fit_stage")   <- nm
      attr(mod, "conv_msgs")   <- NA_character_
      attr(mod, "is_singular") <- NA
      return(mod)
    }
    
    mod <- suppressWarnings(
      suppressMessages(
        try(lmer(f, data = df, control = ctrl), silent = TRUE)
      )
    )
    
    if (inherits(mod, "try-error")) next
    
    if (!bad_fit(mod)) {
      attr(mod, "fit_stage")   <- nm
      attr(mod, "conv_msgs")   <- get_conv_msgs(mod)
      attr(mod, "is_singular") <- isSingular(mod, tol = tol_sing)
      return(mod)
    }
  }
  
  # absolute last-resort safeguard
  mod <- lm(forms$lm_fixed, data = df)
  attr(mod, "fit_stage")   <- "lm_fixed"
  attr(mod, "conv_msgs")   <- NA_character_
  attr(mod, "is_singular") <- NA
  mod
}


root = "/home/jules/Documents/VCCP_JULES/Analyses/precompute/EXPORT_DF/patient_wise/Pxx"
outputdir_fig = "/home/jules/Documents/VCCP_JULES/Analyses/results/Pxx/LMM/diagnosis"
outputdir_df_lmm = "/home/jules/Documents/VCCP_JULES/Analyses/precompute/EXPORT_DF/patient_wise/Pxx/LMM_res/COND"

sujet = "LH018"
protocol_type = "HV"

band_list <- c("theta", "alpha", "beta", "gamma")
phase_list <- c("inspi", "expi")

phase_cycle_sel = phase_list[1]
band = band_list[1]

for (band in band_list) {
  
  filename <- paste0("/", sujet, "_df_Pxx_", band, "_HV_R.xlsx")
  df_raw <- read_excel(paste(root, filename, sep = "/"))
  
  df_raw$chan <- paste0(df_raw$sujet, "_", df_raw$chan)
  
  ROI_list <- unique(df_raw$ROI)
  
  for (phase_cycle_sel in phase_list) {
    
    df_phase <- subset(df_raw, phase_cycle == phase_cycle_sel)
    
    ROI_sel = ROI_list[1]
    
    for (ROI_sel in ROI_list) {
      
      tryCatch({
        
        print(ROI_sel)
        
        filename_export_diagnostic <- paste(protocol_type, band, phase_cycle_sel, ROI_sel, sep = "_")
        
        df_oneROI <- subset(df_phase, ROI == ROI_sel)
        
        # Optional: counts per subject (quick QA)
        df_count <- df_oneROI %>%
          group_by(sujet) %>%
          summarise(
            n_chan  = n_distinct(chan),
            n_cycle = n_distinct(cycle),
            .groups = "drop"
          )
        print(df_count)
        
        # Keep only required columns
        df <- df_oneROI[, c("sujet", "chan", "cycle", "Pxx", "cond")]
        
        # Factors / reference levels
        df$sujet <- factor(df$sujet)
        df$chan  <- factor(df$chan)
        df$cond  <- factor(df$cond)
        
        df$cond  <- relevel(df$cond,  ref = "RB")
        
        # ---- MODEL
        model <- fit_model_with_fallback(df)
        
        fit_stage <- attr(model, "fit_stage")
        conv_msgs <- attr(model, "conv_msgs")
        is_sing   <- attr(model, "is_singular")
        
        message("✅ Model fit stage: ", fit_stage,
                if (!is.na(is_sing)) paste0(" | singular=", is_sing) else "",
                if (!is.na(conv_msgs)) paste0(" | conv=", conv_msgs) else "")
        
        # ---- stats for histogram subtitle
        skew_chan <- round(skewness(df$Pxx), 2)
        kurt_chan <- round(kurtosis(df$Pxx), 2)
        
        # ------------------------------------------------------------
        # PLOT 1: subject-wise boxplot (ggplot)
        # ------------------------------------------------------------
        p_box <- ggplot(df, aes(x = sujet, y = Pxx, color = sujet, fill = sujet)) +
          geom_boxplot(width = .2, alpha = .5, outlier.alpha = 0,
                       position = position_dodge(.9)) +
          stat_summary(fun = median, geom = "point", size = 2,
                       position = position_dodge(.9), color = "white") +
          labs(
            title = paste(filename_export_diagnostic, "Pxx — Subject-wise", sep = " | ")
          ) +
          theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "none"
          )
        
        # ------------------------------------------------------------
        # PLOT 2: histogram (ggplot)
        # ------------------------------------------------------------
        p_hist <- ggplot(df, aes(x = Pxx)) +
          geom_histogram(bins = 30, fill = "lightblue", color = "white") +
          labs(
            title    = "Histogram",
            subtitle = paste("kurtosis:", kurt_chan, "| skewness:", skew_chan),
            x        = "Pxx values",
            y        = "Count"
          ) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
        
        # ------------------------------------------------------------
        # PLOT 3: QQ plot of residuals (ggplot)
        # ------------------------------------------------------------
        res_df <- data.frame(res = resid(model))
        
        p_qq <- ggplot(res_df, aes(sample = res)) +
          stat_qq() +
          stat_qq_line() +
          labs(
            title = "QQ plot (residuals)",
            subtitle = paste0("fit_stage=", fit_stage,
                              if (!is.na(is_sing)) paste0(" | singular=", is_sing) else "")
          ) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
        
        # ------------------------------------------------------------
        # PLOT 4: Residuals vs Fitted (homoskedasticity check)
        # Equivalent to: plot(fitted(m7), resid(m7))
        # ------------------------------------------------------------
        
        diag_df <- data.frame(
          fitted = as.numeric(fitted(model)),
          resid  = as.numeric(resid(model))
        )
        
        p_homo <- ggplot(diag_df, aes(x = fitted, y = resid)) +
          geom_point(alpha = 0.35, size = 1) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
          geom_smooth(method = "loess", se = FALSE, color = "black") +
          labs(
            title = "Residuals vs Fitted",
            subtitle = "Homoskedasticity check (look for constant spread, no funnel)",
            x = "Fitted values",
            y = "Residuals"
          ) +
          theme(
            plot.title = element_text(hjust = 0.5)
          )
        
        # ------------------------------------------------------------
        # COMBINE into ONE figure (patchwork)
        # Layout: boxplot on top, hist + qq below
        # ------------------------------------------------------------
        p_all <- (p_box + p_hist) / (p_qq + p_homo) +
          plot_annotation(
            title = paste("Diagnostics:", filename_export_diagnostic),
            theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
          )
        
        print(p_all)
        # ------------------------------------------------------------
        # SAVE one single PNG
        # ------------------------------------------------------------
        file_diag <- paste("DIAGNOSTIC", filename_export_diagnostic, "Pxx.png", sep = "_")
        
        ggsave(
          filename = paste(outputdir_fig, file_diag, sep = "/"),
          plot     = p_all,
          width    = 12,
          height   = 8,
          dpi      = 150
        )
        
        
        # ------------------------------------------------------------
        # EXPORT MODEL RESULTS
        # - tab_model for HTML-like report
        # - broom.mixed fixed effects + add fit diagnostics
        # ------------------------------------------------------------
        
        tab_model(model, show.re.var = TRUE, show.icc = TRUE, show.r2 = TRUE, show.se = TRUE)
        
        model_df <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
          mutate(
            fit_stage   = fit_stage,
            singular    = if (!is.na(is_sing)) is_sing else NA,
            conv_msgs   = if (!is.na(conv_msgs)) conv_msgs else NA,
            band        = band,
            phase_cycle = phase_cycle_sel,
            ROI         = ROI_sel
          )
        
        filesxlsx_ROI <- paste("RES", filename_export_diagnostic, "Pxx.xlsx", sep = "_")
        writexl::write_xlsx(model_df, paste(outputdir_df_lmm, filesxlsx_ROI, sep = "/"))
        
      }, error = function(e) {
        message("❌ Error for ", band, " / ", phase_cycle_sel, " / ", ROI_sel, " : ", conditionMessage(e))
      })
      
    }
  }
}










