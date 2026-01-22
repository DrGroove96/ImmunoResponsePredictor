####################################################################################################
#                                    ImmunoResponse Predictor 
#                     
####################################################################################################

library(shiny)
library(glmnet)
library(data.table)

## --- Bioconductor annotation libs for ID mapping ---
library(org.Hs.eg.db)
library(AnnotationDbi)

options(shiny.maxRequestSize = 40 * 1024^2)

###############################################################
# Helpers
###############################################################

loadTable <- function(file, transpose = FALSE, convertToMatrix = TRUE,
                      sep = ",", header = TRUE) {
  data <- read.csv(file, sep = sep, header = header,
                   row.names = 1, check.names = FALSE)
  if (transpose) data <- t(data)
  if (convertToMatrix) data <- as.matrix(data)
  return(data)
}

# last-resort NA scrubber
scrub_na <- function(M, fill = 0) {
  M[is.na(M)] <- fill
  M
}

# Generic cleaning of gene IDs:
# - trim, strip transcript version, X1234 -> 1234, uppercase
clean_gene_ids <- function(ids) {
  ids2 <- trimws(ids)
  ids2 <- sub("\\.\\d+$", "", ids2)          # strip .version at end
  ids2 <- sub("^X([0-9]+)$", "\\1", ids2)    # X1234 -> 1234
  toupper(ids2)
}

# Calculate log2(TPM+1) from counts
# counts: matrix with samples as rows, genes as columns
# Returns: log2(TPM+1) normalized matrix
calculate_log2TPMp1 <- function(counts) {
  # Use constant gene length (2000 bp) for all genes
  gene_length <- 2000
  
  # Calculate counts per length (cdl)
  cdl <- counts / gene_length
  
  # Calculate TPM: normalize by row sums (per sample)
  sums <- rowSums(cdl, na.rm = TRUE)
  sums[sums == 0] <- 1  # Avoid division by zero
  tpm <- cdl / sums * 1e6
  
  # Apply log2(TPM + 1)
  log2TPMp1 <- log1p(tpm) / log(2)
  
  return(log2TPMp1)
}

###############################################################
# Build synonym-to-model mapping tables (for renaming cols)
###############################################################

# Map everything → mUC Entrez panel (X<entrez>)
build_mUC_mapping <- function(mUC_entrez, mUC_model_ids) {
  keys <- as.character(mUC_entrez)  # Entrez IDs as character
  
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
    mid <- mUC_model_ids[i]   # e.g. "X9700"
    
    s <- c(
      mid,           # exact model ID (X9700)
      ent            # plain Entrez ("9700")
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
    synonym_to_model = synonym_to_model,   # cleaned synonym -> "X9700"
    model_ids        = mUC_model_ids
  )
}

# Map everything → mRCC ENST panel (ENSTxxx.y) using org.Hs.eg.db
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
  
  # From ENTREZ, pull ALIAS symbols
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

# Generic colname normalizer that only RENAMEs columns, does not touch values
normalize_colnames_generic <- function(cols, mapping) {
  cleaned <- clean_gene_ids(cols)
  out     <- cols
  hits    <- cleaned %in% names(mapping$synonym_to_model)
  out[hits] <- mapping$synonym_to_model[ cleaned[hits] ]
  out
}

###############################################################
# ===================== UI =========================
###############################################################
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
        tags$li(strong("% of applicability (Number of Matches / Total Number of Samples) × 100"))
      ),
      br(),
      textOutput("status"),
      textOutput("predictionsCount")
    )
  )
)

###############################################################
# ===================== SERVER =========================
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
  
  # Cosine-distance helpers
  compute_average_cosine_distances <- function(y1, y2, x) {
    # Convert x to vector if it's a matrix
    if (is.matrix(x)) {
      x <- as.vector(x)
    }
    
    cosine_distance <- function(a, b) {
      # Ensure both are vectors
      a <- as.vector(a)
      b <- as.vector(b)
      
      # Check if vectors have same length
      if (length(a) != length(b)) {
        warning("Vector length mismatch in cosine_distance: a=", length(a), ", b=", length(b))
        return(NA_real_)
      }
      
      norm_a <- sqrt(sum(a^2, na.rm = TRUE))
      norm_b <- sqrt(sum(b^2, na.rm = TRUE))
      
      # Handle zero or NA norms - return NA for invalid cases
      if (is.na(norm_a) || is.na(norm_b) || norm_a == 0 || norm_b == 0) {
        return(NA_real_)
      }
      
      # Calculate cosine similarity, then convert to distance
      cos_sim <- sum(a * b, na.rm = TRUE) / (norm_a * norm_b)
      # Clamp cosine similarity to [-1, 1] to avoid numerical issues
      cos_sim <- max(-1, min(1, cos_sim))
      1 - cos_sim
    }
    
    # Check if groups are empty
    if (nrow(y1) == 0) {
      warning("y1 (responder group) is empty")
      avg_distance_R <- NA_real_
    } else {
      # Ensure x and y1 have matching dimensions
      if (length(x) != ncol(y1)) {
        warning("Dimension mismatch: x length=", length(x), ", y1 ncol=", ncol(y1))
        avg_distance_R <- NA_real_
      } else {
        distances_R <- apply(y1, 1, function(sample) cosine_distance(x, sample))
        avg_distance_R <- mean(distances_R, na.rm = TRUE)
        # If all distances are NA, return NA
        if (is.nan(avg_distance_R) || all(is.na(distances_R))) avg_distance_R <- NA_real_
      }
    }
    
    if (nrow(y2) == 0) {
      warning("y2 (non-responder group) is empty")
      avg_distance_NR <- NA_real_
    } else {
      # Ensure x and y2 have matching dimensions
      if (length(x) != ncol(y2)) {
        warning("Dimension mismatch: x length=", length(x), ", y2 ncol=", ncol(y2))
        avg_distance_NR <- NA_real_
      } else {
        distances_NR <- apply(y2, 1, function(sample) cosine_distance(x, sample))
        avg_distance_NR <- mean(distances_NR, na.rm = TRUE)
        # If all distances are NA, return NA
        if (is.nan(avg_distance_NR) || all(is.na(distances_NR))) avg_distance_NR <- NA_real_
      }
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
    
    # Load annotations to get responder/non-responder groups
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
      x_sample_ml <- ml_test_standardized[i, , drop = FALSE]
      avg_distances_ml <- compute_average_cosine_distances(y1_ml, y2_ml, x_sample_ml)
      ml_cosdist_list[[i]] <- avg_distances_ml[1]  # CosDist_2_Rs for ML-standardized
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
  
  # ORR prior (original - kept for backward compatibility if needed)
  apply_orr_prior <- function(results_df, orr = 0.2) {
    results_df <- results_df[order(results_df$CosDist_2_Rs), ]
    total_samples <- nrow(results_df)
    x <- round(total_samples * orr)
    
    # Same logic as original code: build prior labels and compare to model labels
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
  
  # Calculate % of applicability using ML standardization
  calculate_ML_standardization_applicability <- function(results_df, test_data_matrix, model_gene_ids, sampleIDs) {
    tryCatch({
      # Load ML-standardized training data
      ml_train_path <- file.path("standardized_ML_log2TPMp1_train.csv.gz")
      if (!file.exists(ml_train_path)) {
        # If ML training file doesn't exist, fall back to agreement-based calculation
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
      
      # Calculate ML applicability: percentage of agreement between LogitDA_score_label and iCosinDist_label
      # This represents how well the ML-standardized predictions align with cosine distance labels
      agreement <- sum(results_df$LogitDA_score_label == results_df$iCosinDist_label, na.rm = TRUE)
      total_valid <- sum(!is.na(results_df$LogitDA_score_label) & !is.na(results_df$iCosinDist_label))
      
      if (total_valid > 0) {
        ml_applicability <- round((agreement / total_valid) * 100)
      } else {
        ml_applicability <- 0
      }
      
      ml_applicability
    }, error = function(e) {
      # On error, fall back to agreement-based calculation
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
        
        # Model gene panel (strip any X prefixes)
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
          
          # If some genes missing in test → add zero columns
          if (any(missing_mask)) {
            add_mat <- matrix(0,
                              nrow = nrow(expr.test0_matrix),
                              ncol = sum(missing_mask))
            colnames(add_mat) <- gene_ids_clean[missing_mask]
            expr.test0_matrix <- cbind(expr.test0_matrix, add_mat)
            test_colnames_clean <- c(test_colnames_clean, gene_ids_clean[missing_mask])
            match_idx <- match(gene_ids_clean, test_colnames_clean)
          }
          
          # Build test matrix in model-gene order
          x.test <- expr.test0_matrix[, match_idx, drop = FALSE]
          
          # Align train columns with test order and scale
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
          
          # Scale training data
          train_data <- scale(trainM)
          train_data[is.na(train_data)] <- 0
          
          # Scale test data using training statistics (not independent scaling)
          train_means <- attr(train_data, "scaled:center")
          train_sds <- attr(train_data, "scaled:scale")
          train_sds[train_sds == 0 | is.na(train_sds)] <- 1  # Avoid division by zero
          
          # Apply training statistics to test data
          x.test <- sweep(testM, 2, train_means, "-")
          x.test <- sweep(x.test, 2, train_sds, "/")
          x.test[is.na(x.test)] <- 0
          
          # TRAIN annotations
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
          # ----- mRCC block -----
          train_zip_path <- file.path("models/standardized_QN_TPM_train.csv.gz")
          if (!file.exists(train_zip_path)) stop("Training csv.gz file not found: ", train_zip_path)
          train_data_raw_mRCC <- loadTable(
            file = gzfile(train_zip_path, "rt"),
            transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
          )
          if (!is.matrix(train_data_raw_mRCC)) stop("RCC training data is not a matrix")
          sampleIDs_train <- rownames(train_data_raw_mRCC)
          
          # Per-gene training means (for imputation)
          train_means_vec <- colMeans(
            train_data_raw_mRCC[, intersect(colnames(train_data_raw_mRCC), gene_ids_clean), drop = FALSE],
            na.rm = TRUE
          )
          train_means_full <- setNames(rep(NA_real_, length(gene_ids_clean)), gene_ids_clean)
          common_genes <- intersect(names(train_means_vec), gene_ids_clean)
          train_means_full[common_genes] <- train_means_vec[common_genes]
          
          # Impute missing genes in TEST using training means (fallback 0)
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
          
          # Build matrices in model-gene order (no extra scaling: QN+standardized already)
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
          
          # TRAIN annotations
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
        
        # Validate y1 and y2 are not empty
        if (nrow(y1) == 0) {
          stop("No responder samples found in training data. Cannot compute cosine distances to responder group.")
        }
        if (nrow(y2) == 0) {
          stop("No non-responder samples found in training data. Cannot compute cosine distances to non-responder group.")
        }
        
        # Final guard: ensure no NA in x.test before prediction
        if (anyNA(x.test)) {
          x.test <- scrub_na(x.test, fill = 0)
          output$status <- renderText("Some NA values detected after preprocessing; replaced with 0 to proceed.")
        }
        
        # Check for zero-norm test samples (all zeros)
        test_norms <- sqrt(rowSums(x.test^2, na.rm = TRUE))
        if (any(test_norms == 0, na.rm = TRUE)) {
          warning("Some test samples have zero norm (all zeros). Cosine distances may be NA for these samples.")
        }
        
        # ---------------- Manual logistic prediction ----------------
        x.test <- as.matrix(x.test)
        if (any(is.na(x.test))) stop("NA values detected in x.test before prediction")
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
        # Diagnostic: Check dimensions and sample norms
        cat("y1 dimensions:", nrow(y1), "x", ncol(y1), "\n", file = stderr())
        cat("y2 dimensions:", nrow(y2), "x", ncol(y2), "\n", file = stderr())
        cat("x.test dimensions:", nrow(x.test), "x", ncol(x.test), "\n", file = stderr())
        
        # Check if y1 and y2 have any non-zero values
        y1_norms <- sqrt(rowSums(y1^2, na.rm = TRUE))
        y2_norms <- sqrt(rowSums(y2^2, na.rm = TRUE))
        cat("y1 samples with zero norm:", sum(y1_norms == 0, na.rm = TRUE), "out of", length(y1_norms), "\n", file = stderr())
        cat("y2 samples with zero norm:", sum(y2_norms == 0, na.rm = TRUE), "out of", length(y2_norms), "\n", file = stderr())
        
        results_list <- vector("list", nrow(x.test))
        for (i in seq_len(nrow(x.test))) {
          x_sample <- x.test[i, , drop = FALSE]
          # Convert to vector for cosine distance calculation
          x_sample_vec <- as.vector(x_sample)
          
          # Check if test sample has zero norm
          x_norm <- sqrt(sum(x_sample_vec^2, na.rm = TRUE))
          cat("Sample", i, "- x_norm:", x_norm, "\n", file = stderr())
          if (x_norm == 0) {
            cat("Warning: Test sample", i, "has zero norm (all zeros)\n", file = stderr())
          }
          
          # Verify dimensions match
          if (length(x_sample_vec) != ncol(y1) || length(x_sample_vec) != ncol(y2)) {
            cat("Error: Dimension mismatch - x_sample length:", length(x_sample_vec), 
                ", y1 ncol:", ncol(y1), ", y2 ncol:", ncol(y2), "\n", file = stderr())
            avg_distances <- c(NA_real_, NA_real_)
          } else {
            avg_distances <- compute_average_cosine_distances(y1, y2, x_sample_vec)
            cat("Sample", i, "- CosDist_2_Rs:", avg_distances[1], ", CosDist_2_NRs:", avg_distances[2], "\n", file = stderr())
          }
          # Create data.frame with all columns
          results_list[[i]] <- data.frame(
            CosDist_2_Rs          = avg_distances[1],
            CosDist_2_NRs         = avg_distances[2],
            LogitDA_Score         = pred_prob[i],
            LogitDA_pred_label    = pred_labels[i],
            stringsAsFactors = FALSE
          )
        }
        results_df <- do.call(rbind, results_list)
        results_df$sampleID <- sampleIDs
        
        # Add LogitDA_score_label based on LogitDA_Score thresholds
        results_df$LogitDA_score_label <- sapply(results_df$LogitDA_Score, function(s) {
          if (s >= 0.5002) "R" else if (s < 0.2895) "NR" else NA_character_
        })
        
        results_df <- results_df[match(sampleIDs, results_df$sampleID), ]
        results_df <- results_df[, c("sampleID", setdiff(names(results_df), "sampleID"))]
        stopifnot(all(results_df$sampleID == sampleIDs))
        
        # Replace output names: Moreno -> The53, Kim -> UNC-108
        results_df$sampleID <- gsub("Moreno", "The53", results_df$sampleID)
        results_df$sampleID <- gsub("Kim", "UNC-108", results_df$sampleID)
        
        # ---------------- Compute iCosinDist_label ----------------
        # Determine model type
        model_type <- if (input$modelFile == "logistic-Model-train-muc-test-muc.rds") {
          "mUC"
        } else if (input$modelFile == "logistic-Model-train-rcc-test-rcc.rds") {
          "mRCC"
        } else {
          stop("Unknown model type")
        }
        
        # For mUC: Use ML standardization (IMvigor210)
        # For mRCC: Use original ST results (original cosine distances)
        if (model_type == "mUC") {
          # Apply ML standardization and compute cosine distances for iCosinDist_label
          # mUC uses IMvigor210
          ml_standardization_result <- compute_ML_standardized_cosine_distances(
            x_test = x.test,
            model_type = model_type,
            model_gene_ids = model_gene_ids,
            sampleIDs = sampleIDs
          )
          
          # ORR prior with ML standardization for iCosinDist_label
          orr_results <- apply_orr_prior_ML(
            results_df = results_df,
            ml_cosdist_2_rs = ml_standardization_result$ml_cosdist_2_rs,
            orr = 0.2
          )
        } else {
          # For mRCC: Use original ST results (original cosine distances)
          # ORR prior with original cosine distances
          orr_results <- apply_orr_prior(results_df, orr = 0.2)
        }
        results_df <- orr_results$results_df
        
        # Calculate % of applicability
        # For mUC: Use ML standardization
        # For mRCC: Use original ORR prior percentage (ST approach)
        if (model_type == "mUC") {
          # Calculate % of applicability using ML standardization for mUC
          ml_applicability <- calculate_ML_standardization_applicability(
            results_df = results_df,
            test_data_matrix = expr.test0_matrix,
            model_gene_ids = model_gene_ids,
            sampleIDs = sampleIDs
          )
          applicability_percentage <- ml_applicability
        } else {
          # For mRCC: Use original ORR prior percentage (ST approach)
          applicability_percentage <- orr_results$prior_final_percentage
        }
        
        # Save results with appropriate applicability percentage
        results_df  <- save_and_report_results(
          results_df = results_df,
          prior_final_percentage = applicability_percentage
        )
        
        # Remove internal-only prediction label column from final outputs
        if ("LogitDA_pred_label" %in% names(results_df)) {
          results_df$LogitDA_pred_label <- NULL
        }
        
        # Reorder columns: iCosinDist_label before % of applicability
        col_order <- c("sampleID", "CosDist_2_Rs", "CosDist_2_NRs", "LogitDA_Score", 
                       "LogitDA_score_label", "iCosinDist_label", "% of applicability")
        # Add any remaining columns that might not be in the standard order
        remaining_cols <- setdiff(names(results_df), col_order)
        col_order <- c(col_order, remaining_cols)
        # Keep only columns that exist
        col_order <- col_order[col_order %in% names(results_df)]
        results_df <- results_df[, col_order, drop = FALSE]
        
        predictions(results_df)
        
        # ---------------- UI feedback ----------------
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
      # If only a single sample was uploaded, export a minimal set of columns.
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
    contentType = "text/csv"
  )
  
  output$status <- renderText("")
}

shinyApp(ui = ui, server = server)
