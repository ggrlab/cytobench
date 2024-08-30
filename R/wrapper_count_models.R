#' @title Make hyperparameter-optimized models for all outcomes and feature sets
#' @description
#' Use the list of dataframes in df_list to create hyperparameter-optimized models for all potential dependent variables
#' in dvs_potential. Only one of the dataframes is used in the final model based on the validation performance.
#' The idea is to have a list of dataframes with different independent variables (e.g. "cluster" and "metaCluster")
#' and to select the best set of features based on the validation performance.
#'
#' The independent variables are given either as character vector "ivs" or as a regular expression
#' "ivs_regex" to select the columns from df_list[[1]].
#'
#' The dvs are separated into binary, multi-class and continuous outcomes. Multi-class variables
#' have to be defined in dvs_multiclass.
#'
#' Based on the column "tvt_col" in each dataframe, the data is split into training, validation and test set.
#'
#' The training and validation set are used for hyperparameter optimization. "hparam_n_evaluations" random
#' hyperparameters are tested for each learner. The best model is selected according to the
#' default performance measures for each task type.
#'
#' IMPORTANT: Only ONE of the df_list dataframes is used in the final model. Which model it
#' is is also selected based on the performance on the validation set.
#'
#' After HP optimization, a final model is
#' trained on the training and validation set together.
#'
#' learners_classification holds all learners for classification tasks, learners_regression for regression tasks.
#'
#' The final model is then applied to all the samples together.
#'
#' @param df_list
#' A list of dataframes with the same columns. Each dataframe is used to create a task.
#' Only one of the dataframes is used in the final model, based on the validation performance.
#' @param tvt_col
#' The column in each dataframe that defines the training, validation and test set.
#' Training is used to train models on each task, validation to select the best model and df of df_list.
#' The final model is applied to all samples. Default is "tvt".
#' @param outdir
#' The directory where the results are saved. Default is "res_count_models".
#' If NA, no results are saved.
#' @param dvs_potential
#' The potential dependent variables in the dataframes. For each of them, one final
#' model is trained. Default is the first column name of the first dataframe.
#' @param ivs
#' The independent variables in the dataframes. Default is all columns except the first one.
#' If ivs_regex is given, this is ignored.
#' @param ivs_regex
#' A regular expression to select the independent variables WITHIN EACH dataframe. This is in contrast to ivs, which must be the same for all dataframes.
#' If given, ivs is ignored. Default is "[cC]luster".
#' @param dvs_multiclass
#' The dependent variables that are multi-class. Default is c("ABO", "A_AB", "B_AB", "abo_removedAB").
#' @param seed
#' The seed for the random hyperparameter optimization. Default is 2372876.
#' @param hparam_n_evaluations
#' The number of random hyperparameters to test for each learner. Default is 1000.
#' @param learners_classification
#' A list of classification learners. Default is a ranger learner with 2-20 depth and 500, 1000, 1500, 2000 trees.
#' @param learners_regression
#' A list of regression learners. Default is a ranger learner with 2-20 depth and 500, 1000, 1500, 2000 trees.
#' @export
wrapper_count_models <- function(df_list,
                                 tvt_col = "tvt",
                                 outdir = "res_count_models",
                                 dvs_potential = colnames(df_list[[1]])[1],
                                 ivs = colnames(df_list[[1]])[-1],
                                 ivs_regex = "[cC]luster",
                                 dvs_multiclass = c("ABO", "A_AB", "B_AB", "abo_removedAB"),
                                 seed = 2372876,
                                 hparam_n_evaluations = 1000,
                                 learners_classification = list(
                                     mlr3::lrn(
                                         "classif.ranger",
                                         predict_type = "prob", predict_sets = c("train", "test"),
                                         max.depth = paradox::to_tune(2, 20), # minimum and maximum depth
                                         num.trees = paradox::to_tune(c(500, 1000, 1500, 2000))
                                     )
                                     # Glmnet fails with too little samples in each group. I will skip for this reason
                                     # ,mlr3::lrn(
                                     #     "classif.cv_glmnet",
                                     #     predict_type = "prob", predict_sets = c("train", "test"),
                                     #     alpha = paradox::to_tune(0, 1), s = "lambda.1se"
                                     # )
                                 ),
                                 learners_regression = list(
                                     mlr3::lrn(
                                         "regr.ranger",
                                         predict_sets = c("train", "test"),
                                         max.depth = paradox::to_tune(2, 20), # minimum and maximum depth
                                         num.trees = paradox::to_tune(c(500, 1000, 1500, 2000))
                                     )
                                     # ,mlr3::lrn(
                                     #     "regr.cv_glmnet",
                                     #     predict_sets = c("train", "test"),
                                     #     alpha = paradox::to_tune(0, 1), s = "lambda.1se"
                                     # )
                                 ),
                                 dv_class_positive = NULL) {
    pacman::p_load("mlr3learners")
    pacman::p_load("glmnet")
    pacman::p_load("ranger")
    if (all(is.null(names(learners_classification)))) {
        names(learners_classification) <- sapply(learners_classification, function(x) x$id)
    }
    if (all(is.null(names(learners_regression)))) {
        names(learners_regression) <- sapply(learners_regression, function(x) x$id)
    }

    dir.create(outdir, recursive = TRUE, showWarnings = TRUE)
    if (any(sapply(df_list, function(x) !all(dvs_potential %in% colnames(x))))) {
        stop("Not all dataframes have the potential dependent variable columns")
    }
    outcome_types <- c()
    for (outcome_x in dvs_potential) {
        # Check if binary, multiclass or continuous
        unique_vals <- unique(df_list[[1]][[outcome_x]])
        if (length(unique_vals) == 2) {
            cat("Binary outcome: ", outcome_x, "\n")
            outcome_type <- "binary"
        } else if (length(unique_vals) == 1) {
            cat("Single value outcome: ", outcome_x, "    ------ skipping\n")
            next
        } else if (outcome_x %in% dvs_multiclass) {
            cat("Multiclass outcome: ", outcome_x, "\n")
            outcome_type <- "multiclass"
        } else {
            cat("Continuous outcome: ", outcome_x, "\n")
            outcome_type <- "continuous"
        }
        if (outcome_type != "continuous") {
            if (any(table(df_list[[1]][, c(tvt_col, outcome_x)]) < 2)) {
                warning("Outcome  ", outcome_x, "  does not have at least 2 samples in each tvt group, skipping\n")
                next
            }
        }
        outcome_types <- c(outcome_types, outcome_type)
        names(outcome_types)[length(outcome_types)] <- outcome_x
    }
    # # The following warnings are expected:
    # 1: In eval(ei, envir) :
    #   Outcome  ABO  does not have at least 2 samples in each tvt group, skipping
    # 2: In eval(ei, envir) :
    #   Outcome  A_AB  does not have at least 2 samples in each tvt group, skipping
    # 3: In eval(ei, envir) :
    #   Outcome  B_AB  does not have at least 2 samples in each tvt group, skipping
    # 4: In eval(ei, envir) :
    #   Outcome  APH_HB  does not have at least 2 samples in each tvt group, skipping
    # 5: In eval(ei, envir) :
    #   Outcome  APH_MCHC  does not have at least 2 samples in each tvt group, skipping
    # 6: In eval(ei, envir) :
    #   Outcome  APH_RETI  does not have at least 2 samples in each tvt group, skipping
    # 7: In eval(ei, envir) :
    #   Outcome  abo_ab  does not have at least 2 samples in each tvt group, skipping
    # For each outcome and every feature set (cluster and metaCluster), create a task

    tasklist <- sapply(
        names(outcome_types), function(x) {
            sapply(
                df_list, function(counts_x) {
                    if (!is.na(ivs_regex)) {
                        ivs <- grep(ivs_regex, colnames(counts_x), value = TRUE)
                    }
                    counts_x <- counts_x[, c(tvt_col, x, ivs)]
                    counts_x <- na.omit(counts_x)
                    if (outcome_types[[x]] != "continuous") {
                        positive_label <- NULL
                        if (!all(is.null(dv_class_positive))) {
                            # If the positive label is defined, use it, otherwise infer from the data
                            positive_label <- dv_class_positive[[x]]
                        }
                        new_task <- mlr3::as_task_classif(
                            counts_x[, c(x, ivs)],
                            target = x,
                            positive = positive_label
                        )
                    } else {
                        new_task <- mlr3::as_task_regr(
                            counts_x[, c(x, ivs)],
                            target = x
                        )
                    }
                    new_task$id <- x
                    return(list("task" = new_task, "tvt" = counts_x[[tvt_col]]))
                },
                simplify = FALSE
            )
        },
        simplify = FALSE
    )


    learners_tuned <- sapply(names(tasklist), function(task_XX) {
        task_x <- tasklist[[task_XX]]
        sapply(names(task_x), function(cluster_metacluster) {
            cat("\n\n\n\n\nStarting task  ", task_XX, " on ", cluster_metacluster, "\n")
            task_x <- task_x[[cluster_metacluster]]
            tvt_x <- task_x$tvt
            task_x <- task_x$task
            # Here you can define measures. We go directly with the defaults
            # https://mlr3tuning.mlr-org.com/reference/ti.html#default-measures
            # https://mlr3tuning.mlr-org.com/reference/TuningInstanceBatchMultiCrit.html
            measures <- NULL
            if (task_x$task_type == "classif") {
                current_learners <- learners_classification
                # if (task_x$properties != "multiclass") {
                #     measures <- c(mlr3::msr("classif.auc"))
                # } else {
                #     measures <- c(mlr3::msr("classif.mauc_au1p"))
                # }
            } else {
                current_learners <- learners_regression
                # measures <- c(mlr3::msr("regr.rmse"), mlr3::msr("regr.mse"), mlr3::msr("regr.rsq"))
            }
            resampling_train_validation <- mlr3::rsmp("custom")
            resampling_train_validation$instantiate(
                task_x,
                train = list(which(tvt_x == "train")),
                test = list(which(tvt_x == "validation"))
            )
            # if (
            #     length(resampling_train_validation$test_set(1)) > 50 ||
            #         length(resampling_train_validation$train_set(1)) > 50) {
            #     stop("I expected a maximum of 50 samples in the training and validation set")
            # }
            lapply(current_learners, function(lrn_x) {
                tuning_instance <- mlr3tuning::ti(
                    task = task_x,
                    learner = lrn_x,
                    resampling = resampling_train_validation,
                    terminator = mlr3tuning::trm("evals", n_evals = hparam_n_evaluations),
                    measures = measures
                )
                tuner <- mlr3tuning::tnr("random_search", batch_size = 10)
                cat("Tuning '", format(task_x$formula()), "'  with  ", lrn_x$id, "\n")
                set.seed(seed)
                best_tuned <- tuner$optimize(tuning_instance)
                optimal_learner <- lrn_x$clone()
                optimal_learner$param_set$values <- best_tuned$learner_param_vals[[1]]
                return(list(
                    optimal_learner,
                    learner_performance = best_tuned
                ))
            })
        }, simplify = FALSE)
    }, simplify = FALSE)
    if (!is.null(outdir)) {
        qs::qsave(learners_tuned, file = file.path(outdir, "learners_tuned.qs"))
    }
    ### Select "cluster" or "metaCluster" based on the performance

    ### Train the optimal learner on training + validation data together
    learners_tuned_whichCluster <- lapply(learners_tuned, function(outcome_x) {
        performances <- sapply(outcome_x, function(x) {
            lapply(x, function(y) {
                as.matrix(y$learner_performance[, ncol(y$learner_performance), with = FALSE][[1]])
            })
        }, simplify = FALSE)
        bound <- data.table::rbindlist(performances, idcol = "cluster")
        bound_mat <- as.matrix(bound[, 2:ncol(bound)])
        rownames(bound_mat) <- bound$cluster
        bound_mat <- t(bound_mat)
        perf_max <- apply(bound_mat, 1, function(x) colnames(bound_mat)[which.min(x)])
        selected_clustering <- lapply(names(perf_max), function(model_x) {
            cluster_x <- perf_max[model_x]
            outcome_x[[cluster_x]][[model_x]][[1]]
        })
        names(selected_clustering) <- unname(perf_max)
        selected_clustering
    })
    # Note that you might have multiple list elements called "cluster" or "metaCluster" now.


    ### Now train everything using this hyperparameter-optimized learner
    final_models <- sapply(names(learners_tuned_whichCluster), function(task_XX) {
        task_x <- learners_tuned_whichCluster[[task_XX]]
        models <- list()
        models[[task_XX]] <- list()
        for (i in seq_along(task_x)) {
            clustering_x <- names(task_x)[[i]]
            lrn_x_hparam_optimized <- task_x[[i]]$clone()
            task_data <- tasklist[[task_XX]][[clustering_x]]
            cat("\n\nStarting task  ", task_XX, " on ", clustering_x, " with learner ", lrn_x_hparam_optimized$id, "\n")
            lrn_x_hparam_optimized$train(
                task_data$task,
                row_ids = which(task_data[["tvt"]] %in% c("train", "validation"))
            )
            # n_glmnet <- lrn_x_hparam_optimized$model$glmnet.fit$nobs
            # n_ranger <- lrn_x_hparam_optimized$model$num.samples
            # n_samples <- max(n_glmnet, n_ranger)
            # if (n_samples > 100 || n_samples < 75) {
            #     stop("I expected to use most of the samples in the training and validation set TOGETHER. This must be less or equal to 100, potentially less, but also probably more than 75.")
            # }
            models[[task_XX]][[lrn_x_hparam_optimized$id]] <- list(lrn_x_hparam_optimized$clone())
            names(models[[task_XX]][[lrn_x_hparam_optimized$id]]) <- clustering_x
        }
        return(models[[1]])
    }, simplify = FALSE)
    if (!is.null(outdir)) {
        qs::qsave(final_models, file = file.path(outdir, "final_models.qs"))
    }

    ### Finally, apply the models to all samples - including the test set.
    predictions <- sapply(names(final_models), function(outcome_xx) {
        sapply(names(final_models[[outcome_xx]]), function(model_xx) {
            sapply(names(final_models[[outcome_xx]][[model_xx]]), function(cluster_xx) {
                current_model <- final_models[[outcome_xx]][[model_xx]][[cluster_xx]]
                current_task <- tasklist[[outcome_xx]][[cluster_xx]][["task"]]
                dt_pred <- current_model$predict(current_task) |>
                    data.table::as.data.table()
                orig_cn <- colnames(dt_pred)
                which_probs <- grepl("prob", orig_cn)
                prob_names <- paste0(sub("prob.", "", orig_cn[which_probs]), collapse = "__,__")

                new_cn <- sub("prob.*", "prob", colnames(dt_pred))
                new_cn[grepl("prob", new_cn)] <- paste0("prob.level_", 1:sum(grepl("prob", new_cn)))
                colnames(dt_pred) <- new_cn
                dt_pred[["tvt"]] <- tasklist[[outcome_xx]][[cluster_xx]][["tvt"]]
                dt_pred[["task_type"]] <- current_model$task_type
                dt_pred[["outcome_levels"]] <- prob_names
                return(dt_pred)
            }, simplify = FALSE) |> data.table::rbindlist(idcol = "cluster")
        }, simplify = FALSE) |> data.table::rbindlist(idcol = "model")
    }, simplify = FALSE) |> data.table::rbindlist(idcol = "outcome", fill = TRUE)
    if (!is.null(outdir)) {
        data.table::fwrite(predictions, file = file.path(outdir, "predictions.csv"))
    }
    return(
        list(
            final_models = final_models,
            predictions = predictions,
            tvt_table = table(df_list[[tvt_col]])
        )
    )
}
