####################################################################################################
#                                    ImmunoResponse Predictor 
#                     
####################################################################################################

library(shiny)
library(glmnet)
library(data.table)
library(org.Hs.eg.db)
library(AnnotationDbi)

#You can adjust the size to deploy the app to ShinyApp.io (better under 1GB memory)
#options(shiny.maxRequestSize = 40 * 1024^2)

loadTable <- function(file, transpose = FALSE, convertToMatrix = TRUE,
                      sep = ",", header = TRUE) {
  data <- read.csv(file, sep = sep, header = header,
                   row.names = 1, check.names = FALSE)
  if (transpose) data <- t(data)
  if (convertToMatrix) data <- as.matrix(data)
  return(data)
}

# Remove NA values and replace with 0
scrub_na <- function(M, fill = 0) {
  M[is.na(M)] <- fill
  M
}

# Make sure Gene IDs are in the correct format
clean_gene_ids <- function(ids) {
  ids2 <- trimws(ids)
  ids2 <- sub("\\.\\d+$", "", ids2)          
  ids2 <- sub("^X([0-9]+)$", "\\1", ids2)   
  toupper(ids2)
}

# Create a matrix with samples as rows and genes as columns
# Return log2(TPM+1) normalized matrix
calculate_log2TPMp1 <- function(counts) {
# for our datasets, the most we encounter is 2000, adjust as needed
  gene_length <- 2000

  cdl <- counts / gene_length
  
  # Calculate TPM and normalize by row sums (per sample)
  sums <- rowSums(cdl, na.rm = TRUE)
  sums[sums == 0] <- 1  
  tpm <- cdl / sums * 1e6
  
  # Apply log2(TPM + 1) on the dataset (row)
  log2TPMp1 <- log1p(tpm) / log(2)
  return(log2TPMp1)
}

###############################################################
# Mapping Gene IDs to Model IDs
###############################################################

# Map everything to match the mUC model IDs
build_mUC_mapping <- function(mUC_entrez, mUC_model_ids) {
  keys <- as.character(mUC_entrez)  # we need to convert Entrez IDs to character
  sym   <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = keys,
    column = "SYMBOL", keytype = "ENTREZID",
    multiVals = "first"
  )
  ensg  <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = keys,
    column = "ENSEMBL", keytype = "ENTREZID",
    multiVals = "first"
  )
  alias_list <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = keys,
    column = "ALIAS", keytype = "ENTREZID",
    multiVals = "list"
  )
  
  syn_list <- vector("list", length(keys))
  names(syn_list) <- mUC_model_ids
  
  for (i in seq_along(keys)) {
    ent <- keys[i]
    mid <- mUC_model_ids[i]   # this is the model ID, "X9700"

    s <- c(
      mid,           
      ent
    )
    if (!is.na(sym[i]))  s <- c(s, sym[i])
    if (!is.na(ensg[i])) s <- c(s, ensg[i])
    
    if (!is.null(alias_list[[i]]) && length(alias_list[[i]]) > 0) {
      s <- c(s, alias_list[[i]])
    }
    
    syn_list[[i]] <- unique(clean_gene_ids(s))
  }
  
  all_syn <- unlist(syn_list, use.names = FALSE)
  all_mid <- rep(names(syn_list), times = lengths(syn_list))
  keep    <- !duplicated(all_syn)
  
  synonym_to_model <- all_mid[keep]
  names(synonym_to_model) <- all_syn[keep]
  
  list(
    synonym_to_model = synonym_to_model,
    model_ids        = mUC_model_ids
  )
}

# Map everything to the mRCC model IDs
build_mRCC_mapping <- function(mRCC_model_ids) {
  base <- sub("\\.\\d+$", "", mRCC_model_ids)  
  
  # Use ENSEMBLTRANS to go transcript → gene
  sym <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = base,
    column = "SYMBOL", keytype = "ENSEMBLTRANS",
    multiVals = "first"
  )
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = base,
    column = "ENTREZID", keytype = "ENSEMBLTRANS",
    multiVals = "first"
  )
  ensg <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = base,
    column = "ENSEMBL", keytype = "ENSEMBLTRANS",
    multiVals = "first"
  )
  
  alias_list <- rep(list(NULL), length(base))
  names(alias_list) <- base
  valid_entrez <- !is.na(entrez)
  if (any(valid_entrez)) {
    alias_tmp <- AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = entrez[valid_entrez],
      column = "ALIAS", keytype = "ENTREZID",
      multiVals = "list"
    )
    alias_list[valid_entrez] <- alias_tmp
  }
  
  syn_list <- vector("list", length(mRCC_model_ids))
  names(syn_list) <- mRCC_model_ids
  
  for (i in seq_along(mRCC_model_ids)) {
    mid <- mRCC_model_ids[i]   
    b   <- base[i]           
    
    s <- c(mid, b)
    
    if (!is.na(sym[i]))   s <- c(s, sym[i])
    if (!is.na(ensg[i]))  s <- c(s, ensg[i])
    if (!is.na(entrez[i])) s <- c(s, as.character(entrez[i]))
    
    if (!is.null(alias_list[[i]]) && length(alias_list[[i]]) > 0) {
      s <- c(s, alias_list[[i]])
    }
    
    syn_list[[i]] <- unique(clean_gene_ids(s))
  }
  
  all_syn <- unlist(syn_list, use.names = FALSE)
  all_mid <- rep(names(syn_list), times = lengths(syn_list))
  keep    <- !duplicated(all_syn)
  
  synonym_to_model <- all_mid[keep]
  names(synonym_to_model) <- all_syn[keep]
  
  list(
    synonym_to_model = synonym_to_model,
    model_ids        = mRCC_model_ids
  )
}

# This is to normalize the column names to the model IDs, does not change the values
normalize_colnames_generic <- function(cols, mapping) {
  cleaned <- clean_gene_ids(cols)
  out     <- cols
  hits    <- cleaned %in% names(mapping$synonym_to_model)
  out[hits] <- mapping$synonym_to_model[ cleaned[hits] ]
  out
}

###############################################################
# GUI Panel Instructions
###############################################################
# You can edit mostly this part to anything you want
ui <- fluidPage(
  titlePanel("ImmunoResponse Predictor"),
  sidebarLayout(
    sidebarPanel(
      br(), br(),
      fileInput("testFile", "Upload Test.csv (For mUC or mRCC)", accept = c(".csv")),
      br(),
      selectInput(
        "modelFile", "Select trained model", 
        choices = c(
          "select model" = "",
          "mUC model"  = "logistic-Model-train-muc-test-muc.rds", 
          "mRCC model" = "logistic-Model-train-rcc-test-rcc.rds"
        )
      ),
      br(),
      actionButton("predictButton", "Make predictions"),
      br(), br(), br(), br(), br(),
      uiOutput("downloadUI")
    ),
    mainPanel(
      h4("Instructions:"),
      p("1. ", strong("Upload Test Data:"), "Upload a CSV file containing gene-expression matrix."),
      p("2. ", strong("File Structure:"), "The rows of the file represent individual samples, while the columns correspond to gene expression data."),
      p("3. ", strong("Gene IDs:"), "User can use one of the following gene identifiers in the columns:"),
      tags$ul(
        tags$li(strong("Gene Symbols"), "(e.g., TP53)"),
        tags$li(strong("Entrez Gene IDs"), "(e.g., 7157)"),
        tags$li(strong("Ensembl Gene IDs"), "(e.g., ENSG00000141510)"),
        tags$li(strong("Ensembl Transcript IDs"), "(e.g., ENST00000269305 or ENST00000269305.3)")
      ),
      p("4. ", strong("Select Model:"), "Choose a pre-trained model (either mUC or mRCC) from the dropdown."),
      p("5. ", strong("Generate Predictions:"), "Click the 'Generate predictions' button to process the data."),
      p("6. ", strong("Applicability Metric:"), "After predictions, the system calculates a score based on how closely your data matches the model’s trained biological structure."),
      p("7. ", strong("LogitDA Score Cutoff:"), 
        "LogitDA scores >0.50 (<0.29) were predicted to R(NR), and intermediate scores (0.29 < score < 0.50) were classified as NA."),
      p("8. ", strong("Download Predictions:"), "After predictions are made, you can download a CSV file containing:"),
      tags$ul(
        tags$li(strong("Sample ID")),
        tags$li(strong("Cosine Distances from Rs and NRs groups")), 
        tags$li(strong("LogitDA_Score")), 
        tags$li(strong("LogitDA_score_label")), 
        tags$li(strong("iCosinDist_label")), 
        tags$li(strong("% of applicability"), " (Number of Matches / Total Number of Samples) × 100")
      ),
      br(),
      
      br(),
      textOutput("status"),
      textOutput("predictionsCount")
    )
  )
)

###############################################################
# Connect to the backend of ShinyApp
###############################################################
server <- function(input, output, session) {
  predictions <- reactiveVal(NULL)
  
  # Read uploaded test file
  testData <- reactive({
    req(input$testFile)
    expr.test0 <- fread(
      input$testFile$datapath,
      header = TRUE,
      sep = ",",
      na.strings = c("", "NA")
    )
    sampleIDs <- expr.test0[[1]]                # First column as sample IDs
    expr.test0 <- expr.test0[, -1, with = FALSE]  # Remove sample ID column
    expr.test0_matrix <- as.matrix(expr.test0)
    rownames(expr.test0_matrix) <- sampleIDs
    expr.test0_matrix[is.na(expr.test0_matrix)] <- 0
    list(data = expr.test0_matrix, sampleIDs = sampleIDs)
  })
  
  # This part: we need it to compute the correct cosine distance for SINGLE SAMPLE
  # This will make sure that the single samples do not get NA values for Cosince Distance2 NRs, Rs, etc...
  compute_average_cosine_distances <- function(y1, y2, x) {
    # It transforms the input to a numeric vector
    x <- as.vector(as.numeric(x))
    # Handle NA values
    if (any(is.na(x))) {
      x[is.na(x)] <- 0
    }
    cosine_distance <- function(a, b) {
      # Ensure both are numeric vectors
      a <- as.vector(as.numeric(a))
      b <- as.vector(as.numeric(b))
      # Handle NA values
      if (any(is.na(a))) a[is.na(a)] <- 0
      if (any(is.na(b))) b[is.na(b)] <- 0
      
      norm_a <- sqrt(sum(a^2, na.rm = TRUE))
      norm_b <- sqrt(sum(b^2, na.rm = TRUE))
      if (is.na(norm_a) || is.na(norm_b) || norm_a == 0 || norm_b == 0) return(NA_real_)
      1 - sum(a * b, na.rm = TRUE) / (norm_a * norm_b)
    }
    
    avg_distance_R  <- if (nrow(y1) > 0 && ncol(y1) > 0) {
      distances_R <- apply(y1, 1, function(sample) {
        cosine_distance(x, as.vector(as.numeric(sample)))
      })
      result <- mean(distances_R, na.rm = TRUE)
      if (is.na(result) || is.nan(result)) NA_real_ else result
    } else {
      NA_real_
    }
    
    avg_distance_NR <- if (nrow(y2) > 0 && ncol(y2) > 0) {
      distances_NR <- apply(y2, 1, function(sample) {
        cosine_distance(x, as.vector(as.numeric(sample)))
      })
      result <- mean(distances_NR, na.rm = TRUE)
      if (is.na(result) || is.nan(result)) NA_real_ else result
    } else {
      NA_real_
    }
    
    c(avg_distance_R, avg_distance_NR)
  }
  
  # Apply ML standardization and compute cosine distances for iCosinDist_label
  # mUC uses IMvigor210, mRCC uses IMmotion150
  compute_ML_standardized_cosine_distances <- function(x_test, model_type, model_gene_ids, sampleIDs) {
    # Load ML-standardized training data based on model type
    if (model_type == "mUC") {
      # Use IMvigor210 for mUC
      ml_train_path <- file.path("standardized_ML_log2TPMp1_train.csv.gz")
      # Try to find IMvigor210-specific file, fallback to general ML file
      if (file.exists("standardized_ML_log2TPMp1_train_IMvigor210.csv.gz")) {
        ml_train_path <- file.path("standardized_ML_log2TPMp1_train_IMvigor210.csv.gz")
      }
    } else if (model_type == "mRCC") {
      # Use IMmotion150 for mRCC
      ml_train_path <- file.path("standardized_ML_log2TPMp1_train.csv.gz")
      # Try to find IMmotion150-specific file, fallback to general ML file
      if (file.exists("standardized_ML_log2TPMp1_train_IMmotion150.csv.gz")) {
        ml_train_path <- file.path("standardized_ML_log2TPMp1_train_IMmotion150.csv.gz")
      }
    } else {
      stop("Unknown model type for ML standardization")
    }
    
    if (!file.exists(ml_train_path)) {
      stop("ML-standardized training data not found: ", ml_train_path)
    }
    
    # Load ML-standardized training data
    ml_train_data <- loadTable(
      file = gzfile(ml_train_path, "rt"),
      transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
    )
    
    # For ML standardization: use the ML-standardized training data statistics
    # to standardize the test data to match the same distribution
    train_means_ml <- colMeans(ml_train_data, na.rm = TRUE)
    train_sds_ml <- apply(ml_train_data, 2, sd, na.rm = TRUE)
    train_sds_ml[train_sds_ml == 0 | is.na(train_sds_ml)] <- 1
    
    # Align columns
    test_cols_clean <- sub("^X", "", colnames(x_test))
    train_cols_clean <- sub("^X", "", colnames(ml_train_data))
    common_cols <- intersect(test_cols_clean, train_cols_clean)
    
    if (length(common_cols) == 0) {
      stop("No common columns between test data and ML training data")
    }
    
    test_common_idx <- match(common_cols, test_cols_clean)
    train_common_idx <- match(common_cols, train_cols_clean)
    
    test_subset <- x_test[, test_common_idx[!is.na(test_common_idx)], drop = FALSE]
    train_means_subset <- train_means_ml[train_common_idx[!is.na(train_common_idx)]]
    train_sds_subset <- train_sds_ml[train_common_idx[!is.na(train_common_idx)]]
    
    # Standardize test data using ML training statistics
    ml_test_standardized <- sweep(test_subset, 2, train_means_subset, "-")
    ml_test_standardized <- sweep(ml_test_standardized, 2, train_sds_subset, "/")
    ml_test_standardized[is.na(ml_test_standardized)] <- 0
    
    # Use ML-standardized training data directly (already standardized)
    ml_train_standardized <- ml_train_data[, train_common_idx[!is.na(train_common_idx)], drop = FALSE]
    ml_train_standardized[is.na(ml_train_standardized)] <- 0
    
    # We need to load the true labels for the training data to get the responder/non-responder groups
    if (model_type == "mUC") {
      annot_path <- file.path("models/mUC_response_train.zip")
    } else if (model_type == "mRCC") {
      annot_path <- file.path("models/RCC_response_train.zip")
    } else {
      stop("Unknown model type for annotation loading")
    }
    
    if (!file.exists(annot_path)) {
      stop("Annotation file not found: ", annot_path)
    }
    
    sampleAnnot.train <- read.csv(unz(annot_path, "response_train.csv"))
    
    # Match training samples
    train_sampleIDs <- rownames(ml_train_standardized)
    if (!"RNASEQ_SAMPLE_ID" %in% colnames(sampleAnnot.train)) {
      if ("SampleID" %in% colnames(sampleAnnot.train)) {
        sampleAnnot.train$RNASEQ_SAMPLE_ID <- sampleAnnot.train$SampleID
      } else {
        stop("No sample ID column in annotations")
      }
    }
    
    common_samples <- intersect(train_sampleIDs, sampleAnnot.train$RNASEQ_SAMPLE_ID)
    if (length(common_samples) == 0) {
      stop("No common samples between ML training data and annotations")
    }
    
    train_idx <- match(common_samples, train_sampleIDs)
    ml_train_standardized <- ml_train_standardized[train_idx, , drop = FALSE]
    sampleAnnot.train <- sampleAnnot.train[match(common_samples, sampleAnnot.train$RNASEQ_SAMPLE_ID), ]
    
    # Get responder and non-responder groups from ML-standardized data
    if (!"RESPONSE" %in% colnames(sampleAnnot.train)) {
      if ("Responder" %in% colnames(sampleAnnot.train)) {
        sampleAnnot.train$RESPONSE <- ifelse(sampleAnnot.train$Responder == "R" | sampleAnnot.train$Responder == 1, 1, 0)
      } else {
        stop("No RESPONSE column in annotations")
      }
    }
    
    y1_ml <- ml_train_standardized[sampleAnnot.train$RESPONSE == 1 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
    y2_ml <- ml_train_standardized[sampleAnnot.train$RESPONSE == 0 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
    
    # Compute ML-standardized cosine distances for test samples
    ml_cosdist_list <- vector("list", nrow(ml_test_standardized))
    for (i in seq_len(nrow(ml_test_standardized))) {
      # Extract as vector and convert to numeric (single samples)
      x_sample_ml <- as.vector(as.numeric(ml_test_standardized[i, , drop = TRUE]))
      # Handle any remaining NA values
      if (any(is.na(x_sample_ml))) {
        x_sample_ml[is.na(x_sample_ml)] <- 0
      }
      avg_distances_ml <- compute_average_cosine_distances(y1_ml, y2_ml, x_sample_ml)
      ml_cosdist_list[[i]] <- as.numeric(avg_distances_ml[1])  # CosDist_2_Rs for ML-standardized
    }
    
    ml_cosdist_2_rs <- unlist(ml_cosdist_list)
    
    list(ml_cosdist_2_rs = ml_cosdist_2_rs, ml_test_standardized = ml_test_standardized)
  }
  
  # ORR prior with ML standardization
  apply_orr_prior_ML <- function(results_df, ml_cosdist_2_rs, orr = 0.2) {
    # Order by ML-standardized cosine distances
    results_df <- results_df[order(ml_cosdist_2_rs), ]
    total_samples <- nrow(results_df)
    x <- round(total_samples * orr)
    
    # Build prior labels based on ML-standardized cosine distances
    results_df$CosineDist_prior <- c(rep("R", x), rep("NR", total_samples - x))
    rs_matching <- sum(results_df$LogitDA_pred_label[1:x] == results_df$CosineDist_prior[1:x])
    nrs_matching <- sum(results_df$LogitDA_pred_label[(x + 1):total_samples] == 
                          results_df$CosineDist_prior[(x + 1):total_samples])
    total_matching <- rs_matching + nrs_matching
    prior_final_percentage <- round((total_matching / total_samples) * 100)
    
    # Rename CosineDist_prior column in the final results to iCosinDist_label
    names(results_df)[names(results_df) == "CosineDist_prior"] <- "iCosinDist_label"
    
    list(results_df = results_df, prior_final_percentage = prior_final_percentage)
  }
  
  apply_orr_prior <- function(results_df, orr = 0.2) {
    results_df <- results_df[order(results_df$CosDist_2_Rs), ]
    total_samples <- nrow(results_df)
    x <- round(total_samples * orr)
    
    # Build the prior labels and compare to model labels
    results_df$CosineDist_prior <- c(rep("R", x), rep("NR", total_samples - x))
    rs_matching <- sum(results_df$LogitDA_pred_label[1:x] == results_df$CosineDist_prior[1:x])
    nrs_matching <- sum(results_df$LogitDA_pred_label[(x + 1):total_samples] == 
                          results_df$CosineDist_prior[(x + 1):total_samples])
    total_matching <- rs_matching + nrs_matching
    prior_final_percentage <- round((total_matching / total_samples) * 100)
    
    # Rename CosineDist_prior column in the final results to iCosinDist_label
    names(results_df)[names(results_df) == "CosineDist_prior"] <- "iCosinDist_label"
    
    list(results_df = results_df, prior_final_percentage = prior_final_percentage)
  }
  
  # Calculate % of applicability using ML standardization (Remember: we use ML results to compute the applicability)
  calculate_ML_standardization_applicability <- function(results_df, test_data_matrix, model_gene_ids, sampleIDs) {
    tryCatch({
      # Load ML-standardized training data
      ml_train_path <- file.path("standardized_ML_log2TPMp1_train.csv.gz")
      if (!file.exists(ml_train_path)) {
        agreement <- sum(results_df$LogitDA_score_label == results_df$iCosinDist_label, na.rm = TRUE)
        total_valid <- sum(!is.na(results_df$LogitDA_score_label) & !is.na(results_df$iCosinDist_label))
        if (total_valid > 0) {
          return(round((agreement / total_valid) * 100))
        } else {
          return(0)
        }
      }
      
      ml_train_data <- loadTable(
        file = gzfile(ml_train_path, "rt"),
        transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
      )
      
      # Load ML-standardized test data if available
      ml_test_path <- file.path("standardized_ML_log2TPMp1_test.csv.gz")
      if (file.exists(ml_test_path)) {
        ml_test_data <- loadTable(
          file = gzfile(ml_test_path, "rt"),
          transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
        )
        
        # Match test samples with uploaded sampleIDs
        if (all(sampleIDs %in% rownames(ml_test_data))) {
          ml_test_standardized <- ml_test_data[sampleIDs, , drop = FALSE]
        } else {
          # If exact match not found, use first n rows
          ml_test_standardized <- ml_test_data[1:min(length(sampleIDs), nrow(ml_test_data)), , drop = FALSE]
        }
      } else {
        # If test file doesn't exist, standardize test data using ML training statistics
        train_means_ml <- colMeans(ml_train_data, na.rm = TRUE)
        train_sds_ml <- apply(ml_train_data, 2, sd, na.rm = TRUE)
        train_sds_ml[train_sds_ml == 0 | is.na(train_sds_ml)] <- 1
        
        # Align test data columns with ML training data
        test_cols_clean <- sub("^X", "", colnames(test_data_matrix))
        train_cols_clean <- sub("^X", "", colnames(ml_train_data))
        common_cols <- intersect(test_cols_clean, train_cols_clean)
        
        if (length(common_cols) == 0) {
          # Fall back to agreement-based calculation
          agreement <- sum(results_df$LogitDA_score_label == results_df$iCosinDist_label, na.rm = TRUE)
          total_valid <- sum(!is.na(results_df$LogitDA_score_label) & !is.na(results_df$iCosinDist_label))
          if (total_valid > 0) {
            return(round((agreement / total_valid) * 100))
          } else {
            return(0)
          }
        }
        
        test_common_idx <- match(common_cols, test_cols_clean)
        train_common_idx <- match(common_cols, train_cols_clean)
        
        test_subset <- test_data_matrix[, test_common_idx[!is.na(test_common_idx)], drop = FALSE]
        train_means_subset <- train_means_ml[train_common_idx[!is.na(train_common_idx)]]
        train_sds_subset <- train_sds_ml[train_common_idx[!is.na(train_common_idx)]]
        
        # Standardize test data using ML training statistics
        ml_test_standardized <- sweep(test_subset, 2, train_means_subset, "-")
        ml_test_standardized <- sweep(ml_test_standardized, 2, train_sds_subset, "/")
        ml_test_standardized[is.na(ml_test_standardized)] <- 0
      }
      
      # Definition of applicability: total matches / total number of samples
      agreement <- sum(results_df$LogitDA_score_label == results_df$iCosinDist_label, na.rm = TRUE)
      total_valid <- sum(!is.na(results_df$LogitDA_score_label) & !is.na(results_df$iCosinDist_label))
      
      if (total_valid > 0) {
        ml_applicability <- round((agreement / total_valid) * 100)
      } else {
        ml_applicability <- 0
      }
      
      ml_applicability
    }, error = function(e) {
      agreement <- sum(results_df$LogitDA_score_label == results_df$iCosinDist_label, na.rm = TRUE)
      total_valid <- sum(!is.na(results_df$LogitDA_score_label) & !is.na(results_df$iCosinDist_label))
      if (total_valid > 0) {
        return(round((agreement / total_valid) * 100))
      } else {
        return(0)
      }
    })
  }
  
  save_and_report_results <- function(results_df, prior_final_percentage) {
    results_df$`% of applicability` <- ""
    results_df$`% of applicability`[1] <- prior_final_percentage
    results_df
  }
  
  observeEvent(input$predictButton, {
    req(input$testFile, input$modelFile)
    tryCatch({
      withProgress(message = "Generating predictions", value = 0.1, {
        # ---------------- Load model ----------------
        model_path <- file.path("models", input$modelFile)
        if (!file.exists(model_path)) stop("Model file not found: ", model_path)
        
        bestModel <- readRDS(model_path)
        model_gene_ids <- bestModel$beta@Dimnames[[1]]
        
        # Build mapping depending on model
        if (input$modelFile == "logistic-Model-train-muc-test-muc.rds") {
          mUC_entrez <- as.numeric(sub("^X", "", model_gene_ids))
          mapping <- build_mUC_mapping(mUC_entrez, model_gene_ids)
        } else if (input$modelFile == "logistic-Model-train-rcc-test-rcc.rds") {
          mapping <- build_mRCC_mapping(model_gene_ids)
        } else {
          stop("Please select a valid model.")
        }
        
        # ---------------- Test data ----------------
        test_data <- testData()
        expr.test0_matrix <- test_data$data
        sampleIDs <- test_data$sampleIDs
        
        # Normalize uploaded gene IDs to model IDs
        colnames(expr.test0_matrix) <- normalize_colnames_generic(colnames(expr.test0_matrix), mapping)
        
        gene_ids_clean       <- sub("^X", "", model_gene_ids)
        test_colnames_clean  <- sub("^X", "", colnames(expr.test0_matrix))
        match_idx            <- match(gene_ids_clean, test_colnames_clean)
        missing_mask         <- is.na(match_idx)
        
        # ------------- Branch by model: load TRAIN + align ----------------
        if (input$modelFile == "logistic-Model-train-muc-test-muc.rds") {
          # ----- mUC block (original ZIP-based approach) -----
          train_zip_path <- file.path("models/mUC_log2TPMp1_train.zip")
          if (!file.exists(train_zip_path)) stop("Training zip file not found: ", train_zip_path)
          train_data_raw_mUC <- loadTable(
            file = unz(train_zip_path, "log2TPMp1_train.csv"),
            transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
          )
          if (!is.matrix(train_data_raw_mUC)) stop("mUC training data is not a matrix")
          sampleIDs_train <- rownames(train_data_raw_mUC)
          
          # Handle missing genes in the test data
          if (any(missing_mask)) {
            add_mat <- matrix(0,
                              nrow = nrow(expr.test0_matrix),
                              ncol = sum(missing_mask))
            colnames(add_mat) <- gene_ids_clean[missing_mask]
            expr.test0_matrix <- cbind(expr.test0_matrix, add_mat)
            test_colnames_clean <- c(test_colnames_clean, gene_ids_clean[missing_mask])
            match_idx <- match(gene_ids_clean, test_colnames_clean)
          }
          
          x.test <- expr.test0_matrix[, match_idx, drop = FALSE]
          
          train_cols_order <- match(
            sub("^X", "", colnames(x.test)),
            sub("^X", "", colnames(train_data_raw_mUC))
          )
          if (any(is.na(train_cols_order))) {
            missing_in_train <- colnames(x.test)[is.na(train_cols_order)]
            stop("Model genes missing in mUC training data: ", paste(missing_in_train, collapse = ", "))
          }
          trainM <- train_data_raw_mUC[, train_cols_order, drop = FALSE]
          testM  <- x.test
          
          train_data <- scale(trainM)
          if (nrow(testM) == 1) {
            train_means <- attr(train_data, "scaled:center")
            train_sds <- attr(train_data, "scaled:scale")
            train_sds[train_sds == 0 | is.na(train_sds)] <- 1  
            x.test <- sweep(testM, 2, train_means, "-")
            x.test <- sweep(x.test, 2, train_sds, "/")
            x.test <- as.matrix(x.test)
            rownames(x.test) <- rownames(testM)
            colnames(x.test) <- colnames(testM)
          } else {
            x.test <- scale(testM)
          }
          train_data[is.na(train_data)] <- 0
          x.test[is.na(x.test) | is.nan(x.test)] <- 0
          
          annot_zip_path <- file.path("models/mUC_response_train.zip")
          if (!file.exists(annot_zip_path)) stop("Annotation zip file not found: ", annot_zip_path)
          sampleAnnot.train <- read.csv(unz(annot_zip_path, "response_train.csv"))
          
          if (!"RNASEQ_SAMPLE_ID" %in% colnames(sampleAnnot.train)) {
            stop("No RNASEQ_SAMPLE_ID column in mUC annotations")
          }
          if (!any(sampleIDs_train %in% sampleAnnot.train$RNASEQ_SAMPLE_ID)) {
            stop("No common sample IDs between mUC training data and annotations")
          }
          common_samples <- intersect(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID)
          sample_idx <- which(sampleIDs_train %in% common_samples)
          sampleIDs_train <- sampleIDs_train[sample_idx]
          train_data <- train_data[sampleIDs_train, , drop = FALSE]
          sampleAnnot.train <- sampleAnnot.train[match(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID), ]
          
          train_data_ordered <- train_data[
            match(sampleAnnot.train$RNASEQ_SAMPLE_ID, rownames(train_data)),
            , drop = FALSE
          ]
          y1 <- train_data_ordered[sampleAnnot.train$RESPONSE == 1 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
          y2 <- train_data_ordered[sampleAnnot.train$RESPONSE == 0 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
          
        } else if (input$modelFile == "logistic-Model-train-rcc-test-rcc.rds") {
          train_zip_path <- file.path("models/standardized_QN_TPM_train.csv.gz")
          if (!file.exists(train_zip_path)) stop("Training csv.gz file not found: ", train_zip_path)
          train_data_raw_mRCC <- loadTable(
            file = gzfile(train_zip_path, "rt"),
            transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
          )
          if (!is.matrix(train_data_raw_mRCC)) stop("RCC training data is not a matrix")
          sampleIDs_train <- rownames(train_data_raw_mRCC)
          
          train_means_vec <- colMeans(
            train_data_raw_mRCC[, intersect(colnames(train_data_raw_mRCC), gene_ids_clean), drop = FALSE],
            na.rm = TRUE
          )
          train_means_full <- setNames(rep(NA_real_, length(gene_ids_clean)), gene_ids_clean)
          common_genes <- intersect(names(train_means_vec), gene_ids_clean)
          train_means_full[common_genes] <- train_means_vec[common_genes]
          
          if (any(missing_mask)) {
            add_vals <- train_means_full[missing_mask]
            add_vals[is.na(add_vals)] <- 0
            add_mat <- matrix(
              rep(add_vals, each = nrow(expr.test0_matrix)),
              nrow = nrow(expr.test0_matrix),
              byrow = FALSE
            )
            colnames(add_mat) <- gene_ids_clean[missing_mask]
            expr.test0_matrix <- cbind(expr.test0_matrix, add_mat)
            test_colnames_clean <- c(test_colnames_clean, gene_ids_clean[missing_mask])
            match_idx <- match(gene_ids_clean, test_colnames_clean)
          }
          
          x.test <- expr.test0_matrix[, match_idx, drop = FALSE]
          train_cols_order <- match(
            sub("^X", "", colnames(x.test)),
            sub("^X", "", colnames(train_data_raw_mRCC))
          )
          if (any(is.na(train_cols_order))) {
            missing_in_train <- colnames(x.test)[is.na(train_cols_order)]
            stop("Model genes missing in RCC training data: ", paste(missing_in_train, collapse = ", "))
          }
          train_data <- train_data_raw_mRCC[, train_cols_order, drop = FALSE]
          
          annot_zip_path <- file.path("models/RCC_response_train.zip")
          if (!file.exists(annot_zip_path)) stop("Annotation zip file not found: ", annot_zip_path)
          sampleAnnot.train <- read.csv(unz(annot_zip_path, "response_train.csv"))
          
          if (!"RNASEQ_SAMPLE_ID" %in% colnames(sampleAnnot.train)) {
            stop("No RNASEQ_SAMPLE_ID column in RCC annotations")
          }
          if (!any(sampleIDs_train %in% sampleAnnot.train$RNASEQ_SAMPLE_ID)) {
            stop("No common sample IDs between RCC training data and annotations")
          }
          common_samples <- intersect(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID)
          sample_idx <- which(sampleIDs_train %in% common_samples)
          sampleIDs_train <- sampleIDs_train[sample_idx]
          train_data <- train_data[sampleIDs_train, , drop = FALSE]
          sampleAnnot.train <- sampleAnnot.train[match(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID), ]
          
          train_data_ordered <- train_data[
            match(sampleAnnot.train$RNASEQ_SAMPLE_ID, rownames(train_data)),
            , drop = FALSE
          ]
          y1 <- train_data_ordered[sampleAnnot.train$RESPONSE == 1 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
          y2 <- train_data_ordered[sampleAnnot.train$RESPONSE == 0 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
        } else {
          stop("Please select a valid model.")
        }
        
        if (anyNA(x.test) || any(is.nan(x.test))) {
          x.test[is.na(x.test) | is.nan(x.test)] <- 0
          output$status <- renderText("Some NA/NaN values detected after preprocessing; replaced with 0 to proceed.")
        }
        
        x.test <- as.matrix(x.test)
        if (any(is.na(x.test)) || any(is.nan(x.test))) stop("NA/NaN values detected in x.test before prediction")
        beta <- as.matrix(bestModel$beta)
        a0   <- bestModel$a0
        if (ncol(x.test) != nrow(beta)) {
          stop("Dimension mismatch: x.test has ", ncol(x.test), " columns, beta has ", nrow(beta), " rows")
        }
        pred_prob  <- 1 / (1 + exp(-(x.test %*% beta + a0)))
        pred_prob  <- as.vector(pred_prob)
        pred_class <- ifelse(pred_prob > 0.5, 1, 0)
        pred_labels <- ifelse(pred_class == 1, "R", "NR")
        
        # ---------------- Cosine distance block ----------------
        results_list <- vector("list", nrow(x.test))
        for (i in seq_len(nrow(x.test))) {
          # Extract as vector and convert to numeric (critical for single samples)
          x_sample <- as.vector(as.numeric(x.test[i, , drop = TRUE]))
          # Handle any remaining NA values
          if (any(is.na(x_sample))) {
            x_sample[is.na(x_sample)] <- 0
          }
          avg_distances <- compute_average_cosine_distances(y1, y2, x_sample)
          # Make sure all dataframes columns are numeric
          results_list[[i]] <- data.frame(
            CosDist_2_Rs          = as.numeric(avg_distances[1]),
            CosDist_2_NRs         = as.numeric(avg_distances[2]),
            LogitDA_Score         = as.numeric(pred_prob[i]),
            LogitDA_pred_label    = as.character(pred_labels[i]),
            stringsAsFactors = FALSE
          )
        }
        results_df <- do.call(rbind, results_list)
        results_df$sampleID <- sampleIDs
        
        # Nathan's thershold code
        results_df$LogitDA_score_label <- sapply(results_df$LogitDA_Score, function(s) {
          if (s >= 0.5002) "R" else if (s < 0.2895) "NR" else NA_character_
        })
        
        results_df <- results_df[match(sampleIDs, results_df$sampleID), ]
        results_df <- results_df[, c("sampleID", setdiff(names(results_df), "sampleID"))]
        stopifnot(all(results_df$sampleID == sampleIDs))
        
        # Revision Made by Shen-Han Chiu, we don't call Kim anymore, it's UNC-108 now. 
        # Same goes for Moreno. 
        results_df$sampleID <- gsub("Moreno", "The53", results_df$sampleID)
        results_df$sampleID <- gsub("Kim", "UNC-108", results_df$sampleID)
        
        # Determine model type
        model_type <- if (input$modelFile == "logistic-Model-train-muc-test-muc.rds") {
          "mUC"
        } else if (input$modelFile == "logistic-Model-train-rcc-test-rcc.rds") {
          "mRCC"
        } else {
          stop("Unknown model type")
        }
        
        
        if (model_type == "mUC") {
          # Apply ML standardization and compute cosine distances for iCosinDist_label
          # mUC uses IMvigor210
          ml_standardization_result <- compute_ML_standardized_cosine_distances(
            x_test = x.test,
            model_type = model_type,
            model_gene_ids = model_gene_ids,
            sampleIDs = sampleIDs
          )
          
          orr_results <- apply_orr_prior_ML(
            results_df = results_df,
            ml_cosdist_2_rs = ml_standardization_result$ml_cosdist_2_rs,
            orr = 0.2
          )
        } else {
          # For mRCC: Use original ST results (original cosine distances)
          orr_results <- apply_orr_prior(results_df, orr = 0.2)
        }
        results_df <- orr_results$results_df
        
        # mUC and mRCC still use ST appraoch to compute cosine distance, but ML standardization is for applicability results
        if (model_type == "mUC") {
          ml_applicability <- calculate_ML_standardization_applicability(
            results_df = results_df,
            test_data_matrix = expr.test0_matrix,
            model_gene_ids = model_gene_ids,
            sampleIDs = sampleIDs
          )
          applicability_percentage <- ml_applicability
        } else {
          applicability_percentage <- orr_results$prior_final_percentage
        }
        
        results_df  <- save_and_report_results(
          results_df = results_df,
          prior_final_percentage = applicability_percentage
        )
        
        if ("LogitDA_pred_label" %in% names(results_df)) {
          results_df$LogitDA_pred_label <- NULL
        }
        
        col_order <- c("sampleID", "CosDist_2_Rs", "CosDist_2_NRs", "LogitDA_Score", 
                       "LogitDA_score_label", "iCosinDist_label", "% of applicability")
        remaining_cols <- setdiff(names(results_df), col_order)
        col_order <- c(col_order, remaining_cols)
        col_order <- col_order[col_order %in% names(results_df)]
        results_df <- results_df[, col_order, drop = FALSE]
        predictions(results_df)
        
        # ---------------- Generate Results ----------------
        msg_bits <- c("Predictions generated successfully!")
        if (any(missing_mask)) {
          listed <- paste(head(gene_ids_clean[missing_mask], 8), collapse = ", ")
          suffix <- if (sum(missing_mask) > 8) " ..." else ""
          msg_bits <- c(msg_bits)
        } else {
          msg_bits <- c(msg_bits, "No missing genes; exact match to model panel.")
        }
        msg_bits <- c(msg_bits, sprintf("%% of applicability: %.2f%%",
                                        applicability_percentage))
        output$status <- renderText(paste(msg_bits, collapse = " | "))
        
        output$downloadUI <- renderUI({ downloadButton("downloadPredictions", "Download predictions") })
      })
    }, error = function(e) {
      output$status <- renderText(paste("Error generating predictions:", e$message))
    })
  })
  
  output$predictionsCount <- renderText({
    preds <- predictions()
    if (is.null(preds)) return("Rows in predictions: 0")
    paste("Rows in predictions:", nrow(preds))
  })
  
  output$downloadPredictions <- downloadHandler(
    filename = function() paste("predictions-", Sys.Date(), ".csv", sep = ""),
    content = function(file) {
      preds <- predictions()
      if (is.null(preds)) stop("No predictions available to download. Please generate predictions first.")
      # If only a single sample was uploaded, export the core prediction columns
      if (nrow(preds) == 1) {
        keep_cols <- c("sampleID", "CosDist_2_Rs", "CosDist_2_NRs", "LogitDA_Score", "LogitDA_score_label")
        missing_cols <- setdiff(keep_cols, colnames(preds))
        if (length(missing_cols) > 0) {
          stop("Missing expected columns in predictions: ", paste(missing_cols, collapse = ", "))
        }
        fwrite(preds[, keep_cols, drop = FALSE], file, row.names = FALSE, na = "NA")
      } else {
        fwrite(preds, file, row.names = FALSE, na = "NA")
      }
    },
    # output CSV file
    contentType = "text/csv"
  )
  
  output$status <- renderText("")
}

shinyApp(ui = ui, server = server)
