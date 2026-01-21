#' @title Make hyperparameter-optimized models for all outcomes and feature sets
#' @description
#' Use the list of dataframes in df_list to create hyperparameter-optimized models for all potential dependent variables
#' in dvs_potential. Only one of the dataframes is used in the final model based on the validation performance.
#' The idea is to have a list of dataframes with different independent variables (e.g. "cluster" and "metaCluster")
#' and to select the best set of features based on the validation performance.
#'
#' The independent variables are given either as character vector "ivs" or as a regular expression
#' "ivs_regex" to select the columns from `df_list[[1]]`.
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
#' @param df_list A list of dataframes with the same OUTCOME columns but otherwise different columns are allowed. Each dataframe is used to create a task.
#' Only one of the dataframes is used in the final model, based on the validation performance.
#' @param tvt_col
#' The column in each dataframe that defines the training, validation and test set.
#' Training is used to train models on each task, validation to select the best model and df of df_list.
#' The final model is applied to all samples. Default is "tvt".
#' @param outdir
#' The directory where the results are saved. Default is "res_count_models". If NA, no results are saved.
#' @param dvs_potential The potential dependent variables. For each of them, one final model is trained. Default is the first column name of the first dataframe.
#' @param ivs
#' The independent variables in the dataframes. Default is all columns except the first one.
#' If ivs_regex is given, this is ignored.
#' @param ivs_regex
#' A regular expression to select the independent variables WITHIN EACH dataframe.
#' This is in contrast to ivs, which must be the same for all dataframes.
#' If given, ivs is ignored. Default is `"[cC]luster"` (This would fit columns like
#' 'cluster_1', 'cluster_2', 'metaCluster_1', 'metaCluster_2', etc.).
#' @param dvs_multiclass Dependent variables treated as multi-class. Default: c("ABO", "A_AB", "B_AB", "abo_removedAB").
#' @param seed Random seed for reproducibility. Default: 2372876.
#' @param hparam_n_evaluations Number of random hyperparameter configurations. Default: 1000.
#' @param learners_classification
#' A list of classification learners. Default is a ranger learner with 2-20 depth and 500, 1000, 1500, 2000 trees.
#' @param learners_regression
#' A list of regression learners. Default is a ranger learner with 2-20 depth and 500, 1000, 1500, 2000 trees.
#' @param dv_class_positive
#' A vector of positive labels for the classification tasks. If NULL, the positive label is inferred from the data.
#' @param measures
#' The performance measures for the hyperparameter optimization. If NULL, the default measures are used.
#' @param hpoptimized_final_trainsets Sets used for final model training. Default: c("train", "validation").
#' @param verbose Logical. If TRUE, prints progress messages. Default: TRUE.
#' @export
#' @keywords models
#' @keywords cytometry
#' @examples
#' ff_example <- example_processed()
#' fsom <- do_flowsom_TESTING(ff_example)
#' res <- flowSOM_predict(
#'     flowsom_result = fsom,
#'     flowset = flowCore::flowSet(ff_example)
#' )
#' fake_clusterings <- lapply(res[["ncells_per_x"]], function(x) {
#'     tmp <- lapply(1:25, function(y) {
#'         x[1, -1] <- as.list(sample(1000, ncol(x) - 1))
#'         return(x)
#'     }) |> do.call(what = rbind)
#'     set.seed(40)
#'     tmp[["tvt"]] <- sample(
#'         c("train", "validation", "test", "prospective"),
#'         nrow(tmp),
#'         replace = TRUE
#'     )
#'     tmp[["sample"]] <- sample(c("A", "B"), nrow(tmp), replace = TRUE)
#'     return(tmp)
#' })
#' lapply(fake_clusterings, function(x) {
#'     dplyr::count(x, tvt, sample)
#' })
#' allowed_clusterings <- c("cluster", "metaCluster")
#' finalmodels_predictions <- wrapper_count_models(
#'     df_list = fake_clusterings[allowed_clusterings],
#'     tvt_col = "tvt",
#'     # outdir = file.path(outdir_base, basename(metaclustering_dir), training_device),
#'     outdir = local_tempdir_time(),
#'     dvs_potential = c("sample"),
#'     dvs_multiclass = c(),
#'     ivs_regex = "[cC]luster",
#'     hparam_n_evaluations = 3,
#'     seed = 1372873,
#'     learners_classification = list(
#'         mlr3::lrn(
#'             "classif.ranger",
#'             predict_type = "prob", predict_sets = c("train", "test"),
#'             max.depth = paradox::to_tune(2, 20), # minimum and maximum depth
#'             num.trees = paradox::to_tune(c(500, 1000, 1500, 2000)),
#'             importance = "impurity"
#'         )
#'     ),
#'     dv_class_positive = c("sample" = "B"),
#'     measures = mlr3::msr("classif.logloss")
#' )
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
                                         # paradox::to_tune(minimum depth, maximum depth)
                                         max.depth = paradox::to_tune(2, 20),
                                         num.trees = paradox::to_tune(c(500, 1000, 1500, 2000)),
                                         use_weights = "error"
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
                                         max.depth = paradox::to_tune(2, 20),
                                         num.trees = paradox::to_tune(c(500, 1000, 1500, 2000)),
                                         use_weights = "error"
                                     )
                                 ),
                                 dv_class_positive = NULL,
                                 measures = NULL,
                                 hpoptimized_final_trainsets = c("train", "validation"),
                                 verbose = TRUE) {
    row_ids <- NULL # to avoid R CMD check note about undefined global variable
    # Load required learners
    pacman::p_load("mlr3learners", "glmnet", "ranger")

    # Assign names to learners if unnamed
    if (all(is.null(names(learners_classification)))) {
        names(learners_classification) <- sapply(learners_classification, function(x) x$id)
    }
    if (all(is.null(names(learners_regression)))) {
        names(learners_regression) <- sapply(learners_regression, function(x) x$id)
    }

    # Create output directory if required
    dir.create(outdir, recursive = TRUE, showWarnings = TRUE)

    # Check for required columns
    if (any(sapply(df_list, function(x) !all(dvs_potential %in% colnames(x))))) {
        stop("Not all dataframes have the potential dependent variable columns")
    }

    # Convert input data.tables to data.frames
    df_list <- lapply(df_list, function(x) {
        if (data.table::is.data.table(x) || !is.data.frame(x)) {
            return(as.data.frame(x))
        } else {
            return(x)
        }
    })

    outcome_types <- c()
    for (outcome_x in dvs_potential) {
        # Determine outcome type: binary, multiclass, or continuous
        unique_vals <- unique(df_list[[1]][[outcome_x]])
        if (length(unique_vals) == 2) {
            if (verbose) message(paste0("Binary outcome: ", outcome_x, "\n"))
            outcome_type <- "binary"
        } else if (length(unique_vals) == 1) {
            if (verbose) message(paste0("Single value outcome: ", outcome_x, "    ------ skipping\n"))
            next
        } else if (outcome_x %in% dvs_multiclass) {
            if (verbose) message(paste0("Multiclass outcome: ", outcome_x, "\n"))
            outcome_type <- "multiclass"
        } else {
            if (verbose) message(paste0("Continuous outcome: ", outcome_x, "\n"))
            outcome_type <- "continuous"
        }

        # Check for sufficient samples in train/val for classification
        if (outcome_type != "continuous") {
            if (any(table(df_list[[1]][, c(tvt_col, outcome_x)]) < 2)) {
                warning("Outcome ", outcome_x, " does not have at least 2 samples per tvt group, skipping")
                next
            }
        }
        outcome_types <- c(outcome_types, outcome_type)
        names(outcome_types)[length(outcome_types)] <- outcome_x
    }

    if (length(outcome_types) == 0) {
        stop("No valid outcomes found in dvs_potential")
    }

    # Create mlr3 tasks for each outcome and feature set
    tasklist <- sapply(
        names(outcome_types), function(x) {
            sapply(
                df_list, function(counts_x) {
                    # Select independent variables using regex if specified
                    if (!is.na(ivs_regex)) {
                        ivs <- grep(ivs_regex, colnames(counts_x), value = TRUE)
                    }

                    # Keep only relevant columns: tvt_col, outcome, and ivs
                    counts_x <- counts_x[, c(tvt_col, x, ivs)]

                    # Add a row identifier for later merging predictions
                    if (all(is.null(rownames(counts_x)))) {
                        counts_x[["rowname_original"]] <- seq_len(nrow(counts_x))
                    } else {
                        counts_x <- tibble::rownames_to_column(counts_x, var = "rowname_original")
                    }

                    # Remove rows with NA values
                    counts_x_noNA <- stats::na.omit(counts_x)
                    if (nrow(counts_x_noNA) < nrow(counts_x)) {
                        warnings("Removed ", nrow(counts_x) - nrow(counts_x_noNA), " rows with NAs")
                    }

                    # Create classification or regression task
                    if (outcome_types[[x]] != "continuous") {
                        positive_label <- NULL
                        if (!all(is.null(dv_class_positive))) {
                            # If the positive label is defined, use it, otherwise infer from the data
                            positive_label <- dv_class_positive[[x]]
                        }
                        new_task <- mlr3::as_task_classif(
                            counts_x_noNA[, c("rowname_original", x, ivs)],
                            target = x,
                            positive = positive_label
                        )
                    } else {
                        new_task <- mlr3::as_task_regr(
                            counts_x_noNA[, c("rowname_original", x, ivs)],
                            target = x
                        )
                    }

                    new_task$id <- x
                    new_task$set_col_roles("rowname_original", roles = "name")
                    return(list("task" = new_task, "tvt" = counts_x_noNA[[tvt_col]]))
                },
                simplify = FALSE
            )
        },
        simplify = FALSE
    )

    # Perform hyperparameter tuning for each task and feature set
    learners_tuned <- sapply(names(tasklist), function(task_XX) {
        task_x <- tasklist[[task_XX]]
        sapply(names(task_x), function(cluster_metacluster) {
            if (verbose) message(paste0("\n\n\n\n\nStarting task  ", task_XX, " on ", cluster_metacluster, "\n"))
            task_x <- task_x[[cluster_metacluster]]
            tvt_x <- task_x$tvt
            task_x <- task_x$task
            # Here you can define measures. We go directly with the defaults
            # https://mlr3tuning.mlr-org.com/reference/ti.html#default-measures
            # https://mlr3tuning.mlr-org.com/reference/TuningInstanceBatchMultiCrit.html

            if (all(is.null(measures))) {
                measures <- NULL
            }
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
            # Define resampling split for tuning
            resampling_train_validation <- mlr3::rsmp("custom")
            resampling_train_validation$instantiate(
                task_x,
                train = list(which(tvt_x == "train")),
                test = list(which(tvt_x == "validation"))
            )

            # Tune each learner and store best configuration
            lapply(current_learners, function(lrn_x) {
                tuning_instance <- mlr3tuning::ti(
                    task = task_x,
                    learner = lrn_x,
                    resampling = resampling_train_validation,
                    terminator = mlr3tuning::trm("evals", n_evals = hparam_n_evaluations),
                    measures = measures
                )
                tuner <- mlr3tuning::tnr("random_search", batch_size = 10)
                if (verbose) message(paste0("Tuning '", format(task_x$formula()), "'  with  ", lrn_x$id, "\n"))
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

    # Save tuned learners
    if (!is.null(outdir)) {
        qs2::qsave(learners_tuned, file = file.path(outdir, "learners_tuned.qs"))
    }

    # Select best-performing strategy for each outcome
    learners_tuned_whichCluster <- lapply(learners_tuned, function(outcome_x) {
        performances <- sapply(outcome_x, function(x) {
            lapply(x, function(y) {
                as.matrix(y$learner_performance[, ncol(y$learner_performance), with = FALSE][[1]])
            })
        }, simplify = FALSE)

        # Compare models using last performance column (e.g. loss, AUC, etc.)
        bound <- data.table::rbindlist(performances, idcol = "cluster")
        bound_mat <- as.matrix(bound[, 2:ncol(bound)])
        rownames(bound_mat) <- bound$cluster
        bound_mat <- t(bound_mat)
        perf_max <- apply(bound_mat, 1, function(x) colnames(bound_mat)[which.min(x)])

        # Select learner with best performance
        selected_clustering <- lapply(names(perf_max), function(model_x) {
            cluster_x <- perf_max[model_x]
            outcome_x[[cluster_x]][[model_x]][[1]]
        })
        names(selected_clustering) <- unname(perf_max)
        selected_clustering
    })
    # Note that you might have multiple list elements called "cluster" or "metaCluster" now.

    # Establish the optimal cutoffs based on the validation set
    # BEFORE retraining on training and validation set together. I am unsure if this
    # is the perfect choice, but after training on the training and validation set together,
    # the performance is often perfect(ly overfit). Thus any established cutoff on
    # the following performance is likely to be useless.

    ### Now train on hpoptimized_final_trainsets (often "train" or "train" + "validation") using this hyperparameter-optimized learner
    final_models <- sapply(names(learners_tuned_whichCluster), function(task_XX) {
        task_x <- learners_tuned_whichCluster[[task_XX]]
        models <- list()
        models[[task_XX]] <- list()
        for (i in seq_along(task_x)) {
            clustering_x <- names(task_x)[[i]]
            lrn_hpoptim_train <- task_x[[i]]$clone()
            lrn_hpoptim_final <- task_x[[i]]$clone()
            task_data <- tasklist[[task_XX]][[clustering_x]]

            # Train on
            #    1) train alone
            #    2) then on hpoptimized_final_trainsets (often train + validation)
            lrn_hpoptim_train$train(task_data$task, row_ids = which(task_data[["tvt"]] == "train"))
            lrn_hpoptim_final$train(task_data$task, row_ids = which(task_data[["tvt"]] %in% hpoptimized_final_trainsets))

            # Compute ROC threshold on the "validation" data, regardless of where
            # the final model was trained on.
            final_model <- lrn_hpoptim_final$clone()
            data_validation <- data.table::as.data.table(task_data[["task"]])[task_data[["tvt"]] == "validation"]
            pred_fromtrain_toval <- lrn_hpoptim_train$predict_newdata(data_validation) |>
                data.table::as.data.table()
            if ("LearnerClassif" %in% class(lrn_hpoptim_final)) {
                tasklabels <- names(pred_fromtrain_toval)[4:5]
                names(tasklabels) <- c("positive", "negative") # MLR3 convention
                tasklabels <- sub("prob.", "", tasklabels)
                res_roc <- pROC::roc(
                    response = pred_fromtrain_toval[["truth"]],
                    predictor = pred_fromtrain_toval[[4]],
                    levels = tasklabels[c("negative", "positive")], # pROC convention
                    # If I change the direction here, I have to change it in the predictions as well!!
                    direction = "<"
                ) |>
                    pROC::coords("best", ret = "all", best.method = "closest.topleft")
                res_roc_first <- res_roc[1, ]
                # I cannot change "final_model" directly (R6 class)
                # but I can change the underlying model. This is a bit hacky, but it works.
                final_model[["model"]][["threshold_proc_closest.topleft"]] <- res_roc_first[["threshold"]]
            }
            models[[task_XX]][[final_model$id]] <- list(final_model)
            names(models[[task_XX]][[final_model$id]]) <- clustering_x
        }
        return(models[[1]]) # Todo: only in case of only one single model?
    }, simplify = FALSE)

    # Save final models
    if (!is.null(outdir)) {
        qs2::qsave(final_models, file = file.path(outdir, "final_models.qs"))
    }

    ### Finally, apply the models to all samples - including the test set.
    predictions <- sapply(names(final_models), function(outcome_xx) {
        sapply(names(final_models[[outcome_xx]]), function(model_xx) {
            sapply(names(final_models[[outcome_xx]][[model_xx]]), function(cluster_xx) {
                current_model <- final_models[[outcome_xx]][[model_xx]][[cluster_xx]]
                current_task <- tasklist[[outcome_xx]][[cluster_xx]][["task"]]
                dt_pred <- current_model$predict(current_task) |>
                    data.table::as.data.table()

                # Because of the potential na.omit() when creating the task, the
                # row_ids do NOT match the actual samples, shown in rowname_original
                dt_pred[, sample := c(current_task$data(cols = "rowname_original"))]
                dt_pred[, row_ids := NULL]

                # Postprocess classification predictions
                if ("TaskClassif" %in% class(current_task)) {
                    dt_pred[["class_truth"]] <- dt_pred[["truth"]]
                    dt_pred[["class_response"]] <- dt_pred[["response"]]
                    # Set the predicted response according to the threshold identified through
                    # the ROC curve on the validation set.
                    # dt_pred[[4]] is the positive class

                    if (!"multiclass" %in% current_task$properties) {
                        predicted_negative <- dt_pred[[4]] < current_model[["model"]][["threshold_proc_closest.topleft"]]
                        ll_pred <- levels(dt_pred[["response"]])
                        dt_pred[["class_response"]] <- factor(
                            ifelse(predicted_negative, ll_pred[2], ll_pred[1]),
                            levels = ll_pred
                        )
                    }

                    dt_pred[, c("truth", "response") := NULL]
                    # I do the following label substitution to allow for the same column names after data.table::rbindlist
                    orig_cn <- colnames(dt_pred)
                    which_probs <- grepl("prob", orig_cn)
                    prob_names <- paste0(sub("prob.", "", orig_cn[which_probs]), collapse = "__,__")

                    # Normalize column names for binding
                    new_cn <- sub("prob.*", "prob", orig_cn)
                    new_cn[grepl("prob", new_cn)] <- paste0("prob.level_", 1:sum(grepl("prob", new_cn)))
                    colnames(dt_pred) <- new_cn
                    dt_pred[["outcome_levels"]] <- prob_names
                }

                dt_pred[["tvt"]] <- tasklist[[outcome_xx]][[cluster_xx]][["tvt"]]
                dt_pred[["task_type"]] <- current_model$task_type
                return(dt_pred)
            }, simplify = FALSE) |> data.table::rbindlist(idcol = "cluster")
        }, simplify = FALSE) |> data.table::rbindlist(idcol = "model")
    }, simplify = FALSE) |> data.table::rbindlist(idcol = "outcome", fill = TRUE)

    # Save predictions
    if (!is.null(outdir)) {
        data.table::fwrite(predictions, file = file.path(outdir, "predictions.csv"))
    }

    return(list(
        final_models = final_models,
        predictions = predictions,
        tvt_table = table(df_list[[tvt_col]])
    ))
}
