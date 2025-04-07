# GMM on all features
gmm_all_ <- function(x, y, ...) {
  x_nms <- names(x)
  df_temp <- data.frame() # record error rates to choose the best feature.
  fit_rec <- vector("list", length(x_nms)) # record classification results of all features.
  names(fit_rec) <- x_nms
  model_rec <- vector("list", length(x_nms)) # record model info of all features.
  names(model_rec) <- x_nms
  for (j in 1:length(x_nms)) {
    nm <- x_nms[j]
    print(sprintf("GMM on feature %s ... ", nm))
    tbl_y <- table(y)
    if (min(tbl_y)/max(tbl_y) <= 0.1 & min(tbl_y) == 1) {# give warning for unbalanced data.
      print(sprintf("Warning: unbalanced data: min(tbl_y)/max(tbl_y) = %s & min(tbl_y) = %s",
                    min(tbl_y)/max(tbl_y), min(tbl_y)))
      print(tbl_y)
    }
    if (length(unique(y)) == 2) {# binary
      mod <- mclust::MclustDA(x[[nm]], y, verbose = FALSE)
    } else {# multi-class
      mod <- mclust::MclustDA(x[[nm]], y, modelType = "EDDA", verbose = FALSE)
    }
    # result
    pred_temp <- predict(mod, x[[nm]])
    error_temp <- classError(pred_temp$classification, y)$errorRate
    model_rec[[nm]] <- mod
    fit_rec[[nm]] <- pred_temp
    temp <-  data.frame(name = nm, error = error_temp)
    df_temp <- rbind(df_temp, temp)
  }
  return (list(df_temp = df_temp, fit_rec = fit_rec, model_rec = model_rec))
}

# get significant features from logistic regression model
get_sig_feas_ <- function(model_lr, test_type = "Wald", sig_level = 0.05) {
  if (test_type == "Wald") { # Wald test for binary classification
    z <- summary(model_lr)$coefficients/summary(model_lr)$standard.errors
    p <- (1 - pnorm(abs(z), 0, 1)) * 2
    temp <- names(p)[p < sig_level]
    if ("(Intercept)" %in% temp) temp <- temp[temp != "(Intercept)"]
  }
  if (test_type == "chisq") { # Likelihood ratio test for multi-class classification
    p <- MASS::dropterm(model_lr, trace=FALSE, test="Chisq")["Pr(Chi)"]
    temp <- rownames(p)[p[["Pr(Chi)"]] < sig_level]
    temp <- temp[!is.na(temp)]
  }
  return (temp)
}


# metrics for evaluation
eval_metric_ <- function(row_class, type = "binary", draw_confusion = FALSE, method = "macro", ...) {
  reference <- as.factor(row_class$y_obs)
  #prediction <- as.factor(row_class$y_pred)
  prediction <- factor(row_class$y_pred, levels = levels(reference))
  tbl <- caret::confusionMatrix(table(prediction, reference), reference = reference)
  accuracy <- tbl$overall[["Accuracy"]]
  sen_na_class <- "" # class with NA when calculating sensitivity
  spe_na_class <- "" # class with NA when calculating specificity
  if (type == "binary") {
    df_metric <- data.frame(t(as.numeric(tbl$byClass)))
    names(df_metric) <- names(tbl$byClass)
    auc <- pROC::roc(reference, row_class$y_prob)$auc
  } else {
    if (method == "macro") {
      temp <- apply(tbl$byClass, 2, function(x) mean(x, na.rm = TRUE))
      if (sum(is.na(tbl$byClass[, "Sensitivity"])) > 0) {
        sen_na_class <- paste(gsub("Class: ", "",
                                   names(tbl$byClass[, "Sensitivity"])[is.na(tbl$byClass[, "Sensitivity"])]),
                              collapse = ",")
      }
      if (sum(is.na(tbl$byClass[, "Specificity"])) > 0) {
        spe_na_class <- paste(gsub("Class: ", "",
                                   names(tbl$byClass[, "Specificity"])[is.na(tbl$byClass[, "Specificity"])]),
                              collapse = ",")
      }
      df_metric <- data.frame(t(as.numeric(temp)))
      names(df_metric) <- names(temp)
      auc <- pROC::multiclass.roc(reference, row_class$y_prob)$auc
    }
  }
  df_metric <- data.frame(Accuracy = accuracy, df_metric, AUC = auc,
                          SEN_NA_Class = sen_na_class, SPE_NA_Class = spe_na_class)
  if (draw_confusion) {
    dat <- data.frame(reference, prediction)
    draw_confusion_marix(dat, title = sprintf("Confusion Matrix for %s classification", type))
  }
  return (list(confusion_matrix = tbl, df_metric = df_metric))
}


# draw confusion matrix
draw_confusion_matrix <- function(dat, title, ...) {
  cf_plot <- ggplot(dat, aes(reference, prediction, fill= Freq)) +
    geom_tile() +
    ggtitle(title) +
    geom_text(aes(label=Freq), size = 20) +
    scale_fill_gradient(low="white", high="#009194") +
    theme(plot.title = element_text(lineheight=.8, face="bold", size = 25),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          axis.title=element_text(size=20)) +
    labs(x = "Reference",y = "Prediction")
  print(cf_plot)
}


# Drop columns or rows by name from given data frame
df_drop_ <- function(df, drops, by_col = TRUE) {
  if (by_col)
    df_temp <- df[ , !(colnames(df) %in% drops)]
  else
    df_temp <- df[!(rownames(df) %in% drops), ]
  return (df_temp)
}


#' To split data into train and test sets.
#' @param dat Data in data frame.
#' @param test_size Proportion of test set, default 0.3.
#' @param shuffle Whether or not to shuffle data before splitting, default TRUE.
#' @param stratify If data is split in a stratified fashion, default TRUE.
#' @param target_var Target variable.
#' @param split_var Stratified variable, default target_var.
#' @param seed Seed for shuffling.
#' @return x_train, y_train, x_test, y_test, dat_train, dat_test, seeds
#' @export train_test_split
#' @keywords internal
train_test_split <- function(dat, test_size = 0.2, shuffle = TRUE, stratify = TRUE,
                             target_var = "target", split_var = NULL, seed = NULL) {
  n <- nrow(dat)
  if (length(unique(rownames(dat))) != n) {
    warnings("row names of data must be unique and therefore transformed to 1-n.")
    rownames(dat) <- 1:n
  }
  row_names <- rownames(dat)
  train_size = 1 - test_size
  train_ind <- character(0)
  test_ind <- character(0)
  seeds <- character(0)
  if (stratify) {
    if (is.null(split_var)) split_var <- target_var
    dat_sp <- split(dat, dat[, split_var])
    row_ind <- split(row_names, dat[, split_var])
    for (i in 1:length(dat_sp)) {
      temp_ind <- row_ind[[i]]
      if (is.null(seed)) seed <- round(runif(1, 1, 1000))
      set.seed(seed)
      train_i <- sample(temp_ind, round(length(temp_ind)*train_size))
      test_i <- setdiff(temp_ind, train_i)
      train_ind <- c(train_ind, train_i)
      test_ind <- c(test_ind, test_i)
      seeds <- paste(seeds, seed, sep = "_")
    }
  } else {
    if (is.null(seed)) seed <- round(runif(1, 1, 1000))
    set.seed(seed)
    train_ind <- sample(row_names, round(train_size*n))
    test_ind <- setdiff(row_names, train_ind)
    seeds <- paste(seeds, seed, sep = "_")
  }
  if (shuffle) {
    set.seed(seed)
    train_ind <- sample(train_ind, length(train_ind))
    set.seed(seed)
    test_ind <- sample(test_ind, length(test_ind))
  }
  # results
  dat_train <- dat[train_ind, ]
  x_train <- df_drop_(dat_train, target_var)
  y_train <- dat_train[, target_var]
  dat_test <- dat[test_ind, ]
  x_test <- df_drop_(dat_test, target_var)
  y_test <- dat_test[, target_var]

  return (list(dat_train = dat_train, x_train = x_train, y_train = y_train,
               dat_test = dat_test, x_test = x_test, y_test = y_test,
               seeds = seeds))
}
