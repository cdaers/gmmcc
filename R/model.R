#' Classifier using GMM chains
#'
#' @param x Predictors in data frame.
#' @param y Target variable.
#' @param n size of data, length(y).
#' @param thres Accuracy threshold.
#' @param proc_name In training or testing phase.
#' @param graph_dir Directory for storing graphics.
#' @param model_gmm Trained GMM classifier.
#' @param gmm_args A list of parameters for GMM chains.
#' \itemize{
#'  \item num_chains: Number of GMM chains.
#'  \item chain_length: Maximum iterations.
#'  \item draw_de: If to draw GMM based density estimator of the classifier in GMM chain.
#' }
#' @param ... dots argument.
#'
#' @return A list with elements.
#' \itemize{
#'  \item row_class_gmm: classification info by GMM Classifier of all rows in data frame.
#'  \item model_gmm: GMM Classifier trained.
#' }
#'
#' @import mclust
#' @export gmm_classifier
#' @keywords internal
gmm_classifier <- function(x, y, n, thres,
                           proc_name = "train",
                           graph_dir = sprintf("%s/graphics", getwd()),
                           model_gmm = NULL,
                           gmm_args = list(num_chains = 3L, chain_length = 3L, draw_de = TRUE),
                           ...) {
  if ("num_chains" %in% names(gmm_args)) {
    num_chains <- gmm_args$num_chains
  } else {
    print("Number of GMM chains is not given and set to 3 by default!")
    num_chains = 3L
  }
  if ("chain_length" %in% names(gmm_args)) {
    chain_length <- gmm_args$chain_length
  } else {
    print("Length of GMM chain is not given and set to 3 by default!")
    chain_length = 3L
  }
  if ("draw_de" %in% names(gmm_args)) {
    draw_de <- gmm_args$draw_de
  } else {
    print("Argument draw_de is not given and set to TRUE by default!")
    draw_de <- TRUE
  }

  if (num_chains > dim(x)[2]) {
    print(sprintf("Classifier can not have more chains than the number of features. num_chains set to %s.",
                  dim(x)[2]))
    num_chains <- dim(x)[2]
  }

  if (proc_name == "test") {
    if (is.null(model_gmm)) {
      stop("Trained model must be given in testing phase!")
    }
    num_chains <- length(model_gmm)
  }


  if (proc_name == "train") {
    model_gmm <- vector("list", num_chains)
    res0 <- gmm_all_(x = x, y = y)
  }

  if (draw_de) {
    graph_dir_gmm <- sprintf("%s/gmm", graph_dir)
    if (!dir.exists(graph_dir_gmm)) {
      dir.create(graph_dir_gmm, showWarnings = TRUE, recursive = TRUE)
    }
  }

  row_class_gmm <- data.frame() # record classification info of all rows.

  x_copy <- x
  y_copy <- y

  for (chain_num in 1L:num_chains) {
    print(sprintf("GMM chain %s in %s phase", chain_num, proc_name))
    if (proc_name == "test") {
      chain_length <- length(model_gmm[[chain_num]]) # Maximal iterations for this chain in training phase.
    }
    for (iter in 1:chain_length) {
      iter_next <- if (iter == chain_length) FALSE else TRUE # If next iteration is needed.
      if (iter == 1) {
        x <- x_copy
        y <- y_copy
      }

      print(sprintf("--------- %s-iter ---------", iter))
      if (proc_name == "train") {
        if (iter == 1) {
          res <- res0
          x_nm <- res$df_temp[which(res$df_temp$error == sort(res$df_temp$error)[chain_num])[1], "name"]
        } else {
          res <- gmm_all_(x = x, y = y)
          x_nm <- res$df_temp[which(res$df_temp$error == min(res$df_temp$error))[1], "name"]
        }
        y_pred <- res$fit_rec[[x_nm]]$classification
        y_prob <- apply(res$fit_rec[[x_nm]]$z, 1, function(x) max(x))
        model_gmm[[chain_num]][[iter]] <- res$model_rec[[x_nm]]
        names(model_gmm[[chain_num]])[iter]<- x_nm
      }
      if (proc_name == "test") {
        x_nm <- names(model_gmm[[chain_num]])[iter]
        if (!is.null(model_gmm[[chain_num]][[iter]])) {
          pred_temp <- mclust::predict.MclustDA(model_gmm[[chain_num]][[iter]], x[[x_nm]])
          y_pred <- pred_temp$classification
          y_prob <- apply(pred_temp$z, 1, function(x) max(x))
        }
      }
      gmm_error <- classError(y_pred, y)$errorRate
      y_misclass <- y[y != y_pred] # misclassified points
      print(sprintf("GMM on feature %s with accuracy rate %s at iteration %s.", x_nm,
                    1 - length(y_misclass)/n, iter))
      if (draw_de) {
        grDevices::png(sprintf("%s/%s_%s_%s_%s.png", graph_dir_gmm, x_nm, proc_name, iter, chain_num))
        plot(mclust::densityMclust(x[[x_nm]]), main = sprintf("MBDE on %s in %s set", x_nm, proc_name),
             ylab = "MBDE", xlab = x_nm, what = "density")
        graphics::points(x[[x_nm]], rep(0,length(x[[x_nm]])), pch = 19, col = y + 2)
        grDevices::dev.off()
      }
      tbl_ymisclass <- table(y_misclass)
      if (length(y_misclass)/n <= 1 - thres) {
        print(sprintf("GMM chain suspended: accuracy threshold %s already reached!", thres))
        iter_next <- FALSE
      }
      if (length(tbl_ymisclass) == 1) {# Only points from one single class left!
        print(sprintf("GMM chain suspended: only (%s) points from class %s left.",
                      length(y_misclass), unique(y_misclass)))
        iter_next <- FALSE
      }
      if (min(tbl_ymisclass) < 2) {# empirical
        print(sprintf("Warning: %s point from class %s left.", min(tbl_ymisclass),
                      names(tbl_ymisclass)[tbl_ymisclass == min(tbl_ymisclass)][1]))

        print(tbl_ymisclass)
        if (length(unique(y_misclass)) == 2) {
          print(sprintf("GMM chain suspended since there are %s classes to separate.", length(unique(y_misclass))))
          iter_next <- FALSE
        } else if (mean(tbl_ymisclass) <= 1.5) {
          print(sprintf("GMM chain suspended since average number of points in each class is %s.",
                        mean(tbl_ymisclass)))
          iter_next <- FALSE
        }
      }

      x_correct <- if (iter_next) x[y == y_pred, ] else x
      y_correct <- if (iter_next) y[y == y_pred] else y
      pred_correct <- if (iter_next) y_pred[y == y_pred] else y_pred
      prob_correct <- if (iter_next) y_prob[y == y_pred] else y_prob
      if (length(y_correct) > 0) {
        row_class_iter <- data.frame(chain_num = chain_num,
                                    row_name = rownames(x_correct),
                                    y_obs = y_correct,
                                    y_pred = pred_correct,
                                    y_prob = prob_correct,
                                    error_rate = gmm_error,
                                    model = "gmm",
                                    feature = x_nm,
                                    iter_next = iter_next,
                                    iter = iter,
                                    proc_name = proc_name)
        row_class_gmm <- rbind(row_class_gmm, row_class_iter)
      }

      if (!iter_next) break

      # unclassified points for next iteration
      x <- x[which(y != y_pred), ]
      y <- y[which(y != y_pred)]
    }
  }
  return (list(row_class_gmm = row_class_gmm,
               model_gmm = model_gmm))
}


#' Model training
#'
#' @param x Predictors in data frame.
#' @param y Target variable.
#' @param features_ignored_gmm Features ignored by GMM Classifier, default NULL.
#' @param model_names Model names, default linear regression (lr) and GMM classifier (gmm).
#' @param thres Accuracy threshold, default 0.98.
#' @param min_thres accuracy rate needed, default 0.92.
#' @param graph_dir Directory where graphics are stored.
#'
#' @param lr_params A list of parameters for Logistic Regression model.
#' \itemize{
#'  \item penalty: If select features using LASSO, default FALSE.
#'  \item strict_selection: If select features strictly, e.g. features with non-zero
#'  coefficients in LASSO model or features with p-value <= sig_level in logistic regression,
#' default FALSE.
#'  \item test_type: Significant test for selecting features from LR model, default "Wald".
#' When set to "Wald", a normal approximation is used,
#' otherwise a more accurate likelihood ratio test is carried out.
#'  \item sig_level: Significant level, default 0.05.
#'  \item min_features: Minimal number of features to be selected, default 5L.
#' }
#'
#' @param gmm_args A list of parameters for GMM chains.
#' \itemize{
#'  \item num_chains: Number of GMM chains.
#'  \item chain_length: Maximum iterations.
#' }
#' @param ... dots argument.
#' @return A list of elements.
#' \itemize{
#'  \item model_trained: Mixture model trained.
#'  \item row_class: Classification info of all rows in data frame.
#' }
#' @import caret
#' @import glmnet
#' @import nnet
#' @export model_train
#' @keywords internal
model_train <- function(x, y,
                        features_ignored_gmm = NULL,
                        model_names = c("lr", "gmm"),
                        thres = 0.98,
                        min_thres = 0.92,
                        graph_dir = sprintf("%s/graphics", getwd()),
                        lr_params = list(penalty = FALSE, strict_selection = FALSE, # parameters for LR model
                                         test_type = "Wald", sig_level = 0.05, min_features = 5L),
                        gmm_args = list(num_chains = 3L, chain_length = 3L), # parameters for GMM chain
                        ...) {
  if (!"penalty" %in% names(lr_params)) lr_params$penalty = FALSE
  if (!"strict_selection" %in% names(lr_params)) lr_params$strict_selection = FALSE
  if (!"test_type" %in% names(lr_params)) lr_params$test_type = "Wald"
  if (!"sig_level" %in% names(lr_params)) lr_params$sig_level = 0.05
  if (!"min_features" %in% names(lr_params)) lr_params$min_features = 5
  if (!"trace" %in% names(lr_params)) lr_params$trace = FALSE
  if (!"MaxNWts" %in% names(lr_params)) lr_params$MaxNWts = 2000



  if (length(intersect(model_names, c("lr", "gmm"))) == 0) {
    stop("model_names should be one of c('lr', 'gmm') or both.")
  }
  if ("gmm" %in% model_names) {
    gmm_run <- TRUE
  }
  if (lr_params$min_features > ncol(x)) {
    print("Warning: Number of minimal features exceeds ncol(x). min_features set to %s.")
    lr_params$min_features <- dim(x)[2]
  }

  model_trained <- vector("list", length(model_names))
  names(model_trained) <- model_names
  n <- length(y) # number of rows
  m <- length(x) # number of features
  row_class <- data.frame() # record classification info in data frame.

  x_copy <- x
  y_copy <- y

  if ("lr" %in% model_names) {
    gmm_only <- FALSE
    family <- if (length(unique(y)) == 2) "binomial" else "multinomial"

    # feature selection
    if (lr_params$penalty) {
      if (family == "binomial") {
        model_lr0 <- glmnet::glmnet(as.matrix(x), y, family = family,
                                    lambda=glmnet::cv.glmnet(as.matrix(x), y)$lambda.1se)
        coef_temp <- glmnet::coef.glmnet(model_lr0)[, 1]
        x_selected <- names(coef_temp)[coef_temp != 0]
      } else {
        model_lr0 <- glmnet::glmnet(as.matrix(x), y, family = family, #type.multinomial = "grouped",
                                    lambda=glmnet::cv.glmnet(as.matrix(x), y)$lambda.1se)
        coef_temp0 <- glmnet::coef.glmnet(model_lr0)
        coef_temp <- coef_temp0[1][[1]][, 1]
        x_selected <- names(coef_temp)[coef_temp != 0]
        for (k in 2:length(coef_temp0)) {
          coef_temp <- coef_temp0[k][[1]][, 1]
          x_selected <- union(x_selected, names(coef_temp)[coef_temp != 0])
        }
      }
      x_selected <- x_selected[x_selected != ""]
    } else {
      dat_train <- data.frame(target = as.factor(y), x)
      model_lr0 <- nnet::multinom(target ~. , data = dat_train, trace = lr_params$trace,
                                  MaxNWts = lr_params$MaxNWts)
      x_selected <- get_sig_feas_(model_lr0, test_type = lr_params$test_type,
                                  sig_level = lr_params$sig_level) # get significant features from LR model
    }
    if (!lr_params$strict_selection & length(x_selected) > 0) {
      var_imp <- if (!lr_params$penalty) caret::varImp(model_lr0)
                  else caret::varImp(model_lr0, lambda = model_lr0$lambda)
      x_selected <- rownames(var_imp)[match(sort(var_imp[["Overall"]], decreasing = TRUE),
                                            var_imp[["Overall"]])][1L:lr_params$min_features]
    }
    if (length(x_selected) == 0) {
      print("Coefficients of all variables in LR are 0. Trained Model consists of GMM Classifier only.")
      gmm_only <- TRUE
    }

    # logistic regression model based on selected features
    if (length(x_selected) > 0) {
      x_lr <- x[x_selected]
      df_selected <- data.frame(target = as.factor(y), x_lr)
      model_lr <- nnet::multinom(target ~. , data = df_selected, trace = lr_params$trace,
                                 MaxNWts = lr_params$MaxNWts)
      model_trained[["lr"]] <- model_lr
      model_trained[["x_selected"]] <- x_selected
      y_pred <- stats::predict(model_lr, as.matrix(x_lr), type = "class")
      if (family == "binomial") {
        y_prob <- as.numeric(stats::predict(model_lr, as.matrix(x_lr), type = "probs"))
      } else {
        y_prob <- as.numeric(apply(stats::predict(model_lr, as.matrix(x_lr), type = "probs"), 1, max))
      }
      # identify correctly classified points and calculate accuracy rate
      y_correct <- y[y == y_pred]
      y_misclass <- y[y != y_pred]
      lr_accuracy <- length(y_correct)/length(y)
      print(sprintf("Accuracy of Logistic Regression on training data: %s.", lr_accuracy))
      if (lr_accuracy >= thres) {
        print(sprintf("Training phase finished with Logistic Regression: accuracy rate %s already reached!", thres))
        gmm_run <- FALSE
      } else {
        tbl_ymisclass <- table(y_misclass)
        if (length(tbl_ymisclass) == 1) {
          if (lr_accuracy >= min_thres) {
            print("Minimal accuracy threshold reached. Trained Model consists of Logistic Regression only.")
            gmm_run <- FALSE
          } else {
            print("Minimal accuracy threshold not reached. Trained Model consists of GMM Classifier only.")
            gmm_only <- TRUE
          }
        } else if (min(tbl_ymisclass) < 2) {# Only empirical.
          print(tbl_ymisclass)
          if (length(tbl_ymisclass) == 2 | mean(tbl_ymisclass) <= 1.5) {
            print(sprintf("Only %s point of class %s left for GMM in training phase.",
                          min(tbl_ymisclass), names(tbl_ymisclass)[tbl_ymisclass == min(tbl_ymisclass)][1]))
            if (lr_accuracy >= min_thres) {
              print("Minimal accuracy threshold reached. Trained Model consists of Logistic Regression only.")
              gmm_run <- FALSE
            } else {
              print("Minimal accuracy threshold not reached. Trained Model consists of GMM Classifier only.")
              gmm_only <- TRUE
            }
          }
        }
      }
    }

    if (gmm_only) {
      model_trained[["lr"]] <- NULL
      x_selected <- ""
      x <- x_copy
      y <- y_copy
      gmm_args$num_chains <- 5L
      gmm_args$chain_length <- 7L
    } else {
      x_correct <- if (gmm_run) x[y == y_pred, ] else x
      y_correct <- if (gmm_run) y[y == y_pred] else y
      pred_correct <- if (gmm_run) y_pred[y == y_pred] else y_pred
      prob_correct <- if (gmm_run) y_prob[y == y_pred] else y_prob
      row_class_lr <- data.frame(chain_num = 0,
                                 row_name = rownames(x_correct),
                                 y_obs = y_correct,
                                 y_pred = pred_correct,
                                 y_prob = prob_correct,
                                 error_rate = 1 - lr_accuracy,
                                 model = "lr",
                                 feature = paste(x_selected, collapse = ","),
                                 iter_next = TRUE,
                                 iter = 0,
                                 proc_name = "train")
      row_class <- rbind(row_class, row_class_lr)
      # misclassified points
      x <- x[which(y != y_pred), ]
      y <- y[which(y != y_pred)]
    }
  }

  if (!gmm_run) {# if criterion already fulfilled, return trained model and classification info.
    return (list(model_trained = model_trained, row_class = row_class))
  }

  if ("gmm" %in% model_names) {
    if (!is.null(features_ignored_gmm)) {
        x <- df_drop_(x, features_ignored_gmm)
    }
    model_temp <- gmm_classifier(x = x, y = y, n = n, thres = thres, graph_dir = graph_dir, gmm_args = gmm_args)
    model_trained[["gmm"]] <- model_temp$model_gmm
    row_class_gmm <- model_temp$row_class_gmm
    row_class <- rbind(row_class, row_class_gmm)
  }
  return (list(model_trained = model_trained, row_class = row_class))
}


#' Model testing
#'
#' @param x Predictors in data frame.
#' @param y Target variable.
#' @param features_ignored_gmm Features ignored by GMM Classifier, default NULL.
#' @param thres Accuracy threshold, default 0.98.
#' @param graph_dir Directory where graphics are stored.
#' @param model_trained Trained model from training phase.
#' @param ... dots argument.
#' @return row_class Classification info of all rows in data frame
#' @import glmnet
#' @import nnet
#' @export model_test
#' @keywords internal
model_test <- function(x, y,
                       features_ignored_gmm = NULL,
                       thres = 0.98,
                       model_trained = NULL,
                       graph_dir = sprintf("%s/graphics", getwd()),
                       ...) {
  if (is.null(model_trained)) {
    stop("Trained model must be given in testing phase!")
  }
  model_gmm <- model_trained[["gmm"]]
  gmm_run <- if (!is.null(model_gmm)) TRUE else FALSE

  n <- length(y) # number of rows
  m <- length(x) # number of features
  row_class <- data.frame() # record classification info of all rows in data frame.

  model_lr <- model_trained[["lr"]]
  family <- if (length(unique(y)) == 2) "binomial" else "multinomial"
  if (!is.null(model_lr)) {
    x_lr <- x[, model_trained[["x_selected"]]]
    y_pred <- stats::predict(model_lr, as.matrix(x_lr), type = "class")
    if (family == "binomial") {
      y_prob <- as.numeric(stats::predict(model_lr, as.matrix(x_lr), type = "probs"))
    } else {
      y_prob <- as.numeric(apply(stats::predict(model_lr, as.matrix(x_lr), type = "probs"), 1, max))
    }
    # identify correctly classified points and calculate accuracy rate.
    y_correct <- y[y == y_pred]
    y_misclass <- y[y != y_pred]
    lr_accuracy <- length(y_correct)/length(y)
    print(sprintf("Accuracy of Logistic Regression on test dataset: %s.", lr_accuracy))
    if (lr_accuracy >= thres) {
      print(sprintf("Testing phase finished with Logistic Regression: accuracy rate %s already reached!", thres))
      gmm_run <- FALSE
    }
    tbl_ymisclass <- table(y_misclass)
    if (length(tbl_ymisclass) == 1) {# only points from one single class left, no need to run gmm_classifier!
      print(sprintf("Testing phase finished with Logistic Regression since only (%s) points from class %s left.",
                    length(y_misclass), unique(y_misclass)))
      gmm_run <- FALSE
    } else {
      if (min(tbl_ymisclass) < 2) {
        print(sprintf("Warning: %s point of class %s left for gmm_classifier in testing phase.",
                      min(tbl_ymisclass), names(tbl_ymisclass)[tbl_ymisclass == min(tbl_ymisclass)][1]))
        print(tbl_ymisclass)
        if (length(unique(y_misclass)) == 2) gmm_run <- FALSE
      }
    }
    x_correct <- if (gmm_run) x[y == y_pred, ] else x
    y_correct <- if (gmm_run) y[y == y_pred] else y
    pred_correct <- if (gmm_run) y_pred[y == y_pred] else y_pred
    prob_correct <- if (gmm_run) y_prob[y == y_pred] else y_prob
    row_class_lr <- data.frame(chain_num = 0,
                               row_name = rownames(x_correct),
                               y_obs = y_correct,
                               y_pred = pred_correct,
                               y_prob = prob_correct,
                               error_rate = 1 - lr_accuracy,
                               model = "lr",
                               feature = paste(model_trained[["x_selected"]], collapse = ","),
                               iter_next = TRUE,
                               iter = 0,
                               proc_name = "test")
    row_class <- rbind(row_class, row_class_lr)
    # misclassified points
    x <- x[which(y != y_pred), ]
    y <- y[which(y != y_pred)]
  }

  if (!gmm_run) {# return classification info if criterion fulfilled, .
    return (list(row_class = row_class))
  }

  if (!is.null(features_ignored_gmm)) {
    x <- df_drop_(x, features_ignored_gmm)
  }

  model_temp <- gmm_classifier(x = x, y = y, n = n, thres = thres, model_gmm = model_gmm,
                               proc_name = "test", graph_dir = graph_dir)
  row_class_gmm <- model_temp$row_class_gmm
  row_class <- rbind(row_class, row_class_gmm)

  return (list(row_class = row_class))
}


#' Model run
#'
#' Run LR-GMM model length(rep_labels), default 10 times on given data with random train/test split.
#'
#' @param dat Data to be classified in data frame.
#' @param dat_nm Name of data, default "data".
#'
#' @param target_var Name of target variable, default "target".
#' @param features_ignored_gmm Features ignored by GMM Classifier, default NULL.
#'
#' @param test_size Proportion of test set, default 0.3.
#' @param stratify If data is split in a stratified fashion, default TRUE.
#'
#' @param model_names Model names, default linear regression (lr) and GMM classifier (gmm).
#' @param thres Accuracy threshold, default 0.98.
#' @param min_thres accuracy rate needed, default 0.92.
#'
#' @param rep_labels Labels for the repeated experiments.
#' @param seed Seed for generating random number.
#' @param model_saved If trained model to be saved for stability test, default FALSE.
#' @param row_class_saved If Classification info of all rows saved in data frame, default FALSE.
#'
#' @param result_dir Directory where results are stored.
#' @param graph_dir Directory where graphics are stored.
#'
#' @param opt_metric Metric based on which the optimal chain is selected, default "Accuracy".
#' @param opt_phase Which phase is used for selecting optimal chain, default "test".
#'
#' @param lr_params A list of parameters for Logistic Regression.
#' \itemize{
#'  \item penalty: If select features using LASSO, default FALSE.
#'  \item strict_selection: If select features strictly, e.g. features with non-zero
#' coefficients in LASSO model or features with p-value <= sig_level in logistic regression,
#' default FALSE.
#'  \item test_type: Significant test for selecting features from LR model, default "Wald".
#' When set to "Wald", a normal approximation is used,
#' otherwise a more accurate likelihood ratio test is carried out.
#'  \item sig_level: Significant level, default 0.05.
#'  \item min_features: Minimal number of features to be selected, default 5.
#' }
#'
#' @param gmm_args A list of parameters for GMM Chains.
#' \itemize{
#'  \item num_chains: Number of GMM chains.
#'  \item chain_length: Maximum iterations.
#' }

#' @return A list of elements
#' \itemize{
#'  \item metric_rec: Model metrics like accuracy, precision, recall, f1, roc in data frame.
#'  \item template_rec: Feature templates models returned.
#'  \item result_dir: Directory where results are stored.
#' }
#' @param ... dots argument.
#' @examples
#'
#' ################### demo 1: Diagnosis of breast cancer
#' ## Binary classification
#' library(gmmcc)
#' data(brcancer_diag)
#' # save the current working directory
#' wd <- getwd()
#' # changes it to its parent directory
#' setwd(dirname(getwd()))
#' model_result <- model_run(brcancer_diag, dat_nm = "brcancer_diag",
#'                          lr_params = list(penalty = TRUE, min_features = 5),
#'                          model_saved = TRUE, opt_metric = "AUC")
#' # accuracy rates from training phase vs. testing phase
#' accuracy <- model_result$metric_rec$Accuracy
#' phase <- model_result$metric_rec$phase
#' boxplot(accuracy ~ phase)
#'
#' ################### demo 2: Sales Forecasting
#' ## Multiclass classification
#' data(sale)
#' model_result <- model_run(sale, dat_nm = "sale",
#'                           gmm_args = list(num_chains = 3),
#'                           model_saved = TRUE, opt_metric = "AUC")
#' # accuracy rates from training phase vs. testing phase
#' accuracy <- model_result$metric_rec$Accuracy
#' phase <- model_result$metric_rec$phase
#' boxplot(accuracy ~ phase)
#' # resets the working directory
#' setwd(wd)
#'
#' @export model_run
model_run <- function(dat, dat_nm = "data", # data
                      target_var = "target", features_ignored_gmm = NULL, # features
                      test_size = 0.3, stratify = TRUE, # parameters for train_test_split
                      model_names = c("lr", "gmm"), thres = 0.98, min_thres = 0.92, # model
                      rep_labels = 1L:10L, seed = 2025, model_saved = FALSE, row_class_saved = FALSE,
                      result_dir = sprintf("%s/result/%s", getwd(), dat_nm),
                      graph_dir = sprintf("%s/graphics/%s", getwd(), dat_nm),
                      opt_metric = "Accuracy", opt_phase = "test", # optimization
                      lr_params = list(penalty = FALSE, strict_selection = FALSE, # parameters for LR model
                                       test_type = "Wald", sig_level = 0.05, min_features = 5L),
                      gmm_args = list(num_chains = 3L, chain_length = 3L), # parameters for GMM chains
                      ...) {

  if (length(list(...)) > 0) {
    extra_args <- list(...)
    #graph_dir <- if ("graph_dir" %in% names(extra_args)) extra_args$graph_dir else NULL
  }

  if (!dir.exists(result_dir)) {
    dir.create(result_dir, showWarnings = TRUE, recursive = TRUE)
  }
  metric_rec <- data.frame() # record accuracy, precision, recall, f1, roc.
  template_rec <- data.frame() # record templates returned.
  if (gmm_args$num_chains > (dim(dat)[2] - 1)) {
    print(sprintf("Classifier can not have more chains than the number of features. num_chains reset to %s.",
                  dim(dat)[2] - 1))
    gmm_args$num_chains <- dim(dat)[2] - 1
  }
  type <- if (length(unique(dat[[target_var]])) > 2) "multiple" else "binary"
  # run rep times
  for (rep_label in rep_labels) {
    print(sprintf("--------- run [%s] ---------", rep_label))
    ptm <- proc.time()
    train_test <- train_test_split(dat = dat, test_size = test_size, stratify = stratify, seed = rep_label)
    print(sprintf("---------  training [%s] started ---------", rep_label))
    train_result <- model_train(x = train_test$x_train, y = train_test$y_train,
                                features_ignored_gmm = features_ignored_gmm,
                                lr_params = lr_params, gmm_args = gmm_args,
                                graph_dir = sprintf("%s/rep_%s", graph_dir, rep_label))

    if (model_saved) saveRDS(train_result$model_trained,
                             sprintf("%s/%s_model_%s.rds", result_dir, dat_nm, rep_label))
    if (row_class_saved) utils::write.csv(train_result$row_class,
                             sprintf("%s/%s_model_train_%s.csv", result_dir, dat_nm, rep_label),
                             row.names = F)
    metric_train <- data.frame()
    template_train <- data.frame()
    for (num in 1L:gmm_args$num_chains) {
      row_class_temp <- train_result$row_class[train_result$row_class$chain_num %in% c(0, num), ]
      metric_temp_train <- data.frame(eval_metric_(row_class_temp, type = type)$df_metric,
                                      phase = "train", rep = rep_label, chain_num = num)
      metric_train <- rbind(metric_train, metric_temp_train)
      template_temp_train <- data.frame(feature = paste(unique(row_class_temp$feature), collapse = " ;  "),
                                        phase = "train", rep = rep_label, chain_num = num)
      template_train <- rbind(template_train, template_temp_train)
    }
    if (opt_phase == "train") {
      opt_chain_ind <- which(metric_train[[opt_metric]] == max(metric_train[[opt_metric]]))[1]
    }
    print(sprintf("---------  training [%s] finished ---------", rep_label))
    print(sprintf("---------  testing [%s] started ---------", rep_label))
    test_result <- model_test(x = train_test$x_test, y = train_test$y_test,
                              features_ignored_gmm = features_ignored_gmm,
                              model_trained = train_result$model_trained,
                              graph_dir = sprintf("%s/rep_%s", graph_dir, rep_label))
    if (row_class_saved) utils::write.csv(test_result$row_class,
                                   sprintf("%s/%s_model_test_%s.csv", result_dir, dat_nm, rep_label),
                                   row.names = F)
    metric_test <- data.frame()
    template_test <- data.frame()
    for (num in 1L:gmm_args$num_chains) {
      row_class_temp <- test_result$row_class[test_result$row_class$chain_num %in% c(0, num), ]
      metric_temp_test <- data.frame(eval_metric_(row_class_temp, type = type)$df_metric,
                                     phase = "test", rep = rep_label, chain_num = num)
      metric_test <- rbind(metric_test, metric_temp_test)
      template_temp_test <- data.frame(feature = paste(unique(row_class_temp$feature), collapse = " ;  "),
                                       phase = "test", rep = rep_label, chain_num = num)
      template_test <- rbind(template_test, template_temp_test)
    }
    if (opt_phase == "test") {
      opt_chain_ind <- which(metric_test[[opt_metric]] == max(metric_test[[opt_metric]]))[1]
    }
    print(sprintf("run [%s] costs %s seconds.", rep_label, proc.time()[3] - ptm[3]))
    print(sprintf("---------  test [%s] finished ---------", rep_label))
    metric_rec <- rbind(metric_rec, metric_train[opt_chain_ind, ], metric_test[opt_chain_ind, ])
    template_rec <- rbind(template_rec, template_train[opt_chain_ind, ], template_test[opt_chain_ind, ])
  }
  print("# ! --------- metrics returned --------- ! #")
  print(metric_rec)
  print("# ! --------- templates returned --------- ! #")
  print(template_rec)
  utils::write.csv(metric_rec, sprintf("%s/%s_metric_rec.csv", result_dir, dat_nm), row.names= F)
  utils::write.csv(template_rec, sprintf("%s/%s_template_rec.csv", result_dir, dat_nm), row.names= F)
  return (list(metric_rec = metric_rec, template_rec = template_rec, result_dir = result_dir))
}


#' Model stability test on different train/test sets
#'
#' @param dat Data to be classified in data frame.
#' @param dat_nm Name of data, default "data".
#' @param target_var Name of target variable, default "target".
#' @param rep Number of repeat experiments, default 10.
#' @param seed Seed for generating random number.
#' @param test_sizes Different test sizes, default c(0.2, 0.3, 0.4).
#' @param stratify If data is split in a stratified fashion, default TRUE.
#' @param features_ignored_gmm Features ignored by GMM Classifier.
#' @param result_dir Directory where results are stored.
#' @param model_chosen Model to be tested.
#' @param no_chain Which chain to be tested.
#' @param ... dots argument.
#' @return stab_test_result Result of stability test.
#' @export model_stab_test
#' @keywords internal
model_stab_test <- function(dat, dat_nm = "data", target_var = "target",
                            rep = 10, seed = 2025,
                            test_sizes = c(0.2, 0.3, 0.4), stratify = TRUE,
                            features_ignored_gmm = NULL, result_dir = paste0(getwd(), "/result"),
                            model_chosen, no_chain,
                            ...) {
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, showWarnings = TRUE, recursive = TRUE)
  }
  if (!is.null(model_chosen[["gmm"]])) {
    num_chains = length(model_chosen[["gmm"]])
    if (sum(1L:num_chains == no_chain) == 0) {
      print(sprintf("%s chain does not exist in model given.", no_chain))
    } else {
      model_chosen[["gmm"]] <- model_chosen[["gmm"]][no_chain]
    }
  }
  stab_test_result <- data.frame() # record accuracy, precision, recall, f1, roc.
  type <- if (length(unique(dat[[target_var]])) > 2) "multiple" else "binary"
  # run rep times
  for (i in 1:rep) {
    print(sprintf("--------- Stability test: run [%s] ---------", i))
    ptm <- proc.time()
    for (test_size in test_sizes) {
      print(sprintf("--------- test size [%s] ---------", test_size))
      train_test <- train_test_split(dat = dat, test_size = test_size, stratify = stratify, seed = i)
      test_result <- model_test(x = train_test$x_test, y = train_test$y_test,
                                features_ignored_gmm = features_ignored_gmm, model_trained = model_chosen)
      row_class_temp <- test_result$row_class
      metric_temp <- data.frame(eval_metric_(row_class_temp, type = type)$df_metric,
                                     test_size = test_size, rep = i, chain_num = no_chain)
      stab_test_result <- rbind(stab_test_result, metric_temp)
    }
    print(sprintf("Stability test: run [%s] costs %s seconds.", i, proc.time()[3] - ptm[3]))
    print(sprintf("---------  Stability test: run [%s] finished ---------", i))
  }
  utils::write.csv(stab_test_result, sprintf("%s/%s_model_stab_test.csv", result_dir, dat_nm), row.names= F)
  return (stab_test_result)
}


#' Run model stability test
#'
#' Run stability test on the best model from the result of model_run.
#'
#' @param dat Data to be classified in data frame.
#' @param dat_nm Name of data, default "data".
#' @param model_result Result from model_run.
#' @param ref_phase Reference phase, only "test" supported.
#' @param eval_inds Evaluation indicators used to choose reference model.
#' @param model_ref Reference model, default NULL.
#' @param ... dots argument.
#' @return A list of elements
#' \itemize{
#'  \item result: Result of stability tests.
#'  \item template_chosen: Feature template of the chosen model.
#' }
#'
#' @examples
#'
#' ################### demo 1: Diagnosis of breast cancer
#' ## Binary classification
#' library(gmmcc)
#' data(brcancer_diag)
#' # save the current working directory
#' wd <- getwd()
#' # changes it to its parent directory
#' setwd(dirname(getwd()))
#' model_result <- model_run(brcancer_diag, dat_nm = "brcancer_diag",
#'                          lr_params = list(penalty = TRUE, min_features = 5),
#'                          model_saved = TRUE, opt_metric = "AUC")
#' # model stability test
#' test_result <- stab_test_run(brcancer_diag, "brcancer_diag", model_result)
#' accuracy <- test_result$result$Accuracy
#' test_size <- test_result$result$test_size
#' boxplot(accuracy ~ test_size)
#'
#' ################### demo 2: Sales Forecasting
#' ## Multiclass classification
#' data(sale)
#' model_result <- model_run(sale, dat_nm = "sale",
#'                           gmm_args = list(num_chains = 3),
#'                           model_saved = TRUE, opt_metric = "AUC")
#' # model stability test
#' test_result <- stab_test_run(sale, "sale", model_result)
#' accuracy <- test_result$result$Accuracy
#' test_size <- test_result$result$test_size
#' boxplot(accuracy ~ test_size)
#' # resets the working directory
#' setwd(wd)
#'
#' @export stab_test_run
stab_test_run <- function(dat, dat_nm,
                          model_result,
                          ref_phase = "test",
                          eval_inds = list(eval1 = "Accuracy", eval2 = "AUC"),
                          model_ref = NULL,
                          ...) {
  if (ref_phase != "test") {
    print("only test phase supported! Reference phase is reset to test phase.")
    ref_phase <- "test"
  }
  template_rec <- model_result$template_rec
  result_dir <- model_result$result_dir
  if (is.null(model_ref)) {
    metric_rec <- model_result$metric_rec
    model_ref <- metric_rec[which(metric_rec$phase == ref_phase), ]
    for (eval_ind in eval_inds) {
      model_ref <- model_ref[model_ref[[eval_ind]] == max(model_ref[[eval_ind]]), ]
    }
  }
  if (nrow(model_ref) > 1) model_ref <- model_ref[sample(nrow(model_ref), 1), ]
  model_chosen <- readRDS(sprintf("%s/%s_model_%s.rds", model_result$result_dir, dat_nm, model_ref$rep))
  template_chosen <- template_rec$feature[template_rec$phase == ref_phase &
                                            template_rec$rep == model_ref$rep &
                                            template_rec$chain_num == model_ref$chain_num]
  result <- model_stab_test(dat, dat_nm, model_chosen = model_chosen, no_chain = model_ref$chain_num)
  return (list(result = result, template_chosen = template_chosen))
}
