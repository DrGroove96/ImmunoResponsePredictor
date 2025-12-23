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

# Map everything → Melanoma Entrez panel (X<entrez>) using org.Hs.eg.db
build_Melanoma_mapping <- function(melanoma_model_ids) {
  # Remove X prefix to get Entrez IDs
  entrez_ids <- sub("^X", "", melanoma_model_ids)
  keys <- as.character(entrez_ids)
  
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
  names(syn_list) <- melanoma_model_ids
  
  for (i in seq_along(keys)) {
    ent <- keys[i]
    mid <- melanoma_model_ids[i]   # e.g. "X6352"
    
    s <- c(
      mid,           # exact model ID (X6352)
      ent            # plain Entrez ("6352")
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
    synonym_to_model = synonym_to_model,   # cleaned synonym -> "X6352"
    model_ids        = melanoma_model_ids
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
          "mRCC model" = "logistic-Model-train-rcc-test-rcc.rds",
          "Melanoma model" = "logistic-Model-train-Melanoma-test-rcc.rds"
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
      p("7. ", strong("Download Predictions:"), "After predictions are made, you can download a CSV file containing:"),
      tags$ul(
        tags$li(strong("Sample ID")),
        tags$li(strong("Cosine Distances from Rs and NRs groups")),
        tags$li(strong("Predicted Response"), "(R or NR)"),
        tags$li(strong("Applicability Score"), "(Cosine Distance)")
      ),
      p("8. ", strong("Decision Rule for LogitDA_Score:"), "Samples with ", code("LogitDA_Score > 0.5"), " are classified as responders (Rs), and all others as non-responders (NRs)."),
      p("9. ", strong("Model Note (mUC):"), "We predict those with Pcd4989g(mUC)."),
      
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
    cosine_distance <- function(a, b) {
      norm_a <- sqrt(sum(a^2, na.rm = TRUE))
      norm_b <- sqrt(sum(b^2, na.rm = TRUE))
      if (is.na(norm_a) || is.na(norm_b) || norm_a == 0 || norm_b == 0) return(NA_real_)
      1 - sum(a * b, na.rm = TRUE) / (norm_a * norm_b)
    }
    avg_distance_R  <- if (nrow(y1) > 0) mean(apply(y1, 1, function(sample) cosine_distance(x, sample)), na.rm = TRUE) else NA
    avg_distance_NR <- if (nrow(y2) > 0) mean(apply(y2, 1, function(sample) cosine_distance(x, sample)), na.rm = TRUE) else NA
    c(avg_distance_R, avg_distance_NR)
  }
  
  # ORR prior
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
        } else if (input$modelFile == "logistic-Model-train-Melanoma-test-rcc.rds") {
          mapping <- build_Melanoma_mapping(model_gene_ids)
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

          train_data <- scale(trainM)
          x.test     <- scale(testM)
          train_data[is.na(train_data)] <- 0
          x.test[is.na(x.test)]         <- 0
          
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
          # ----- mRCC block (using log2TPMp1 normalization) -----
          # Try multiple possible file paths
          train_zip_path <- NULL
          is_counts <- FALSE
          
          # Check for log2TPMp1 file first
          if (file.exists(file.path("models/log2TPMp1_train.csv.gz"))) {
            train_zip_path <- file.path("models/log2TPMp1_train.csv.gz")
            is_counts <- FALSE
          } else if (file.exists(file.path("models/counts_train.csv.gz"))) {
            train_zip_path <- file.path("models/counts_train.csv.gz")
            is_counts <- TRUE
          } else if (file.exists(file.path("models/standardized_QN_TPM_train.csv.gz"))) {
            # Auto-create log2TPMp1 file from existing standardized file
            cat("Creating log2TPMp1_train.csv.gz from standardized_QN_TPM_train.csv.gz...\n")
            withProgress(message = "Creating log2TPMp1 file", value = 0.3, {
              # Load existing file
              temp_data <- fread(
                file.path("models/standardized_QN_TPM_train.csv.gz"),
                header = TRUE,
                sep = ",",
                na.strings = c("", "NA")
              )
              
              # Get sample IDs
              if ("V1" %in% colnames(temp_data)) {
                sample_ids_temp <- temp_data$V1
                temp_data <- temp_data[, -1, with = FALSE]
              } else {
                sample_ids_temp <- rownames(temp_data)
              }
              
              temp_matrix <- as.matrix(temp_data)
              rownames(temp_matrix) <- sample_ids_temp
              temp_matrix[is.na(temp_matrix)] <- 0
              
              # Check if needs conversion (if median > 100, likely counts; if < 5, already log2TPMp1)
              median_val <- median(temp_matrix, na.rm = TRUE)
              if (is.na(median_val)) {
                # If all values are NA, assume already log2TPMp1
                median_val <- 0
              }
              if (!is.na(median_val) && median_val > 100) {
                temp_matrix <- calculate_log2TPMp1(temp_matrix)
              } else if (!is.na(median_val) && median_val >= 5 && median_val <= 100) {
                # Likely TPM, convert to log2TPMp1
                temp_matrix <- log1p(temp_matrix) / log(2)
              }
              # If median < 5, assume already log2TPMp1
              
              temp_matrix[is.na(temp_matrix)] <- 0
              
              # Save as log2TPMp1 file
              temp_dt <- as.data.table(temp_matrix)
              temp_dt <- cbind(sampleID = rownames(temp_matrix), temp_dt)
              fwrite(temp_dt, file.path("models/log2TPMp1_train.csv.gz"), row.names = FALSE, na = "NA")
              cat("Successfully created models/log2TPMp1_train.csv.gz\n")
            })
            
            train_zip_path <- file.path("models/log2TPMp1_train.csv.gz")
            is_counts <- FALSE
          } else {
            stop("Training csv.gz file not found. Expected one of: models/log2TPMp1_train.csv.gz, models/counts_train.csv.gz, or models/standardized_QN_TPM_train.csv.gz")
          }
          
          train_data_raw_mRCC <- loadTable(
            file = gzfile(train_zip_path, "rt"),
            transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
          )
          if (!is.matrix(train_data_raw_mRCC)) stop("RCC training data is not a matrix")
          sampleIDs_train <- rownames(train_data_raw_mRCC)
          
          # If training data is counts, convert to log2TPMp1
          if (is_counts) {
            cat("Converting training counts to log2TPMp1...\n")
            train_data_raw_mRCC <- calculate_log2TPMp1(train_data_raw_mRCC)
          }
          
          # Align test data columns with model genes
          # Normalize uploaded gene IDs to model IDs (already done above, but ensure consistency)
          test_colnames_clean <- sub("^X", "", colnames(expr.test0_matrix))
          match_idx <- match(gene_ids_clean, test_colnames_clean)
          missing_mask <- is.na(match_idx)
          if (length(missing_mask) == 0) missing_mask <- logical(0)
          
          # If test data appears to be counts (large values), convert to log2TPMp1
          # Otherwise assume it's already log2TPMp1 or similar expression values
          if (length(expr.test0_matrix) > 0 && !all(is.na(expr.test0_matrix))) {
            test_max_val <- max(expr.test0_matrix, na.rm = TRUE)
            if (!is.na(test_max_val) && is.finite(test_max_val) && test_max_val > 100) {
              cat("Converting test counts to log2TPMp1...\n")
              expr.test0_matrix <- calculate_log2TPMp1(expr.test0_matrix)
            }
          }
          
          # Per-gene training means (for imputation of missing genes)
          train_means_vec <- colMeans(
            train_data_raw_mRCC[, intersect(colnames(train_data_raw_mRCC), gene_ids_clean), drop = FALSE],
            na.rm = TRUE
          )
          train_means_full <- setNames(rep(NA_real_, length(gene_ids_clean)), gene_ids_clean)
          common_genes <- intersect(names(train_means_vec), gene_ids_clean)
          train_means_full[common_genes] <- train_means_vec[common_genes]
          
          # Impute missing genes in TEST using training means (fallback 0)
          if (length(missing_mask) > 0 && any(missing_mask, na.rm = TRUE)) {
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
            missing_mask <- is.na(match_idx)
          }
          
          # Build matrices in model-gene order
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
          
          # Standardize both training and test data using training statistics
          cat("Standardizing log2TPMp1 data using training statistics...\n")
          train_means <- colMeans(train_data, na.rm = TRUE)
          train_sds <- apply(train_data, 2, sd, na.rm = TRUE)
          train_sds[train_sds == 0 | is.na(train_sds)] <- 1
          
          train_data <- sweep(train_data, 2, train_means, "-")
          train_data <- sweep(train_data, 2, train_sds, "/")
          
          x.test <- sweep(x.test, 2, train_means, "-")
          x.test <- sweep(x.test, 2, train_sds, "/")
          
          train_data[is.na(train_data)] <- 0
          x.test[is.na(x.test)] <- 0
          
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
        } else if (input$modelFile == "logistic-Model-train-Melanoma-test-rcc.rds") {
          # ----- Melanoma block (using Melanoma-specific ST standardized data) -----
          # Try ST file first, fallback to QN if ST doesn't exist
          train_zip_path <- file.path("models/standardized_ST_TPM_train_Melanoma.csv.gz")
          if (!file.exists(train_zip_path)) {
            train_zip_path <- file.path("models/standardized_QN_TPM_train_Melanoma.csv.gz")
            if (!file.exists(train_zip_path)) {
              stop("Training csv.gz file not found. Expected: models/standardized_ST_TPM_train_Melanoma.csv.gz or models/standardized_QN_TPM_train_Melanoma.csv.gz")
            }
          }
          train_data_raw_Melanoma <- loadTable(
            file = gzfile(train_zip_path, "rt"),
            transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
          )
          if (!is.matrix(train_data_raw_Melanoma)) stop("Melanoma training data is not a matrix")
          sampleIDs_train <- rownames(train_data_raw_Melanoma)
          
          # Per-gene training means (for imputation)
          train_means_vec <- colMeans(
            train_data_raw_Melanoma[, intersect(colnames(train_data_raw_Melanoma), gene_ids_clean), drop = FALSE],
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
          
          # Build matrices in model-gene order
          x.test <- expr.test0_matrix[, match_idx, drop = FALSE]
          train_cols_order <- match(
            sub("^X", "", colnames(x.test)),
            sub("^X", "", colnames(train_data_raw_Melanoma))
          )
          if (any(is.na(train_cols_order))) {
            missing_in_train <- colnames(x.test)[is.na(train_cols_order)]
            stop("Model genes missing in Melanoma training data: ", paste(missing_in_train, collapse = ", "))
          }
          
          # Transform test data using ComBat + QN normalization (following pseudo code)
          # Training data is QN normalized and standardized TPM
          # Test data is log2(TPM+1), so we need to:
          # 1. Back-transform from log2(TPM+1) to TPM: 2^x - 1
          # 2. Apply ComBat for batch correction (train vs test)
          # 3. Apply QN normalization
          # 4. Standardize using QN-normalized training statistics
          
          cat("Transforming test data: back-transforming from log2(TPM+1) to TPM...\n")
          x.test_tpm <- 2^x.test - 1
          x.test_tpm[x.test_tpm < 0] <- 0  # Handle numerical precision issues
          
          # Check if required packages are available
          if (!requireNamespace("preprocessCore", quietly = TRUE)) {
            stop("preprocessCore package is required. Please install it: install.packages('BiocManager'); BiocManager::install('preprocessCore')")
          }
          if (!requireNamespace("sva", quietly = TRUE)) {
            stop("sva package is required. Please install it: install.packages('BiocManager'); BiocManager::install('sva')")
          }
          
          # Load original training data
          cat("Loading original training data for ComBat + QN normalization...\n")
          original_train_path <- file.path("standardized_ST_log2TPMp1_train.csv.gz")
          if (!file.exists(original_train_path)) {
            stop("Original training data file not found: ", original_train_path)
          }
          
          original_train_raw <- loadTable(
            file = gzfile(original_train_path, "rt"),
            transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
          )
          # Back-transform original training data from log2(TPM+1) to TPM
          original_train_tpm <- 2^original_train_raw - 1
          original_train_tpm[original_train_tpm < 0] <- 0
          
          # Get gene columns matching model genes (in same order as train_cols_order)
          model_gene_ids <- colnames(train_data_raw_Melanoma)[train_cols_order]
          original_gene_match <- match(
            sub("^X", "", colnames(original_train_tpm)),
            sub("^X", "", model_gene_ids)
          )
          original_train_subset <- original_train_tpm[, !is.na(original_gene_match), drop = FALSE]
          
          # Reorder to match model gene order
          original_col_order <- match(
            sub("^X", "", model_gene_ids),
            sub("^X", "", colnames(original_train_subset))
          )
          original_train_subset <- original_train_subset[, original_col_order[!is.na(original_col_order)], drop = FALSE]
          
          # Ensure test data has same genes in same order
          test_cols_ordered <- match(
            sub("^X", "", colnames(original_train_subset)),
            sub("^X", "", colnames(x.test_tpm))
          )
          x.test_aligned <- x.test_tpm[, test_cols_ordered[!is.na(test_cols_ordered)], drop = FALSE]
          
          # Store original row names for splitting later
          train_rownames <- rownames(original_train_subset)
          test_rownames <- rownames(x.test_aligned)
          
          # Following pseudo code: transpose to probe (gene) x chip (sample) matrix
          cat("Transposing data for ComBat + QN normalization...\n")
          expr.train0.t <- t(original_train_subset)
          expr.test0.t <- t(x.test_aligned)
          
          # Get sample (batch) size
          train_len <- ncol(expr.train0.t)
          test_len <- ncol(expr.test0.t)
          bat <- c(rep("train", train_len), rep("test", test_len))
          
          # Combine training and test data
          data_all <- cbind(expr.train0.t, expr.test0.t)
          
          # ComBat - for batch correction
          cat("Running ComBat for batch correction...\n")
          combat_data <- tryCatch({
            sva::ComBat(
              dat = data_all,
              batch = bat,
              mod = NULL,
              par.prior = TRUE,
              prior.plots = FALSE
            )
          }, error = function(e) {
            warning("ComBat failed: ", e$message, ". Proceeding without ComBat.")
            data_all
          })
          
          # Quantile Normalization - QN
          cat("Performing Quantile Normalization (QN)...\n")
          expr.all <- preprocessCore::normalize.quantiles(combat_data, copy = TRUE)
          
          # Transpose back to samples x genes
          expr.all <- t(expr.all)
          colnames(expr.all) <- rownames(expr.train0.t)  # Gene names
          rownames(expr.all) <- c(colnames(expr.train0.t), colnames(expr.test0.t))  # Sample names
          
          # Split back to training and test sets
          expr.train <- expr.all[rownames(expr.all) %in% train_rownames, , drop = FALSE]
          expr.test <- expr.all[rownames(expr.all) %in% test_rownames, , drop = FALSE]
          
          # Ensure correct order
          expr.train <- expr.train[match(train_rownames, rownames(expr.train)), , drop = FALSE]
          expr.test <- expr.test[match(test_rownames, rownames(expr.test)), , drop = FALSE]
          
          # Calculate QN-normalized training statistics for standardization
          train_means_qn <- colMeans(expr.train, na.rm = TRUE)
          train_sds_qn <- apply(expr.train, 2, sd, na.rm = TRUE)
          train_sds_qn[train_sds_qn == 0 | is.na(train_sds_qn)] <- 1
          
          # Standardize test data using QN-normalized training statistics
          cat("Standardizing QN-normalized test data using training statistics...\n")
          x.test <- sweep(expr.test, 2, train_means_qn, "-")
          x.test <- sweep(x.test, 2, train_sds_qn, "/")
          
          # Reorder back to match original test gene order
          x.test <- x.test[, match(
            sub("^X", "", colnames(x.test_tpm)),
            sub("^X", "", colnames(x.test))
          ), drop = FALSE]
          
          x.test[is.na(x.test)] <- 0
          
          # TRAIN annotations
          annot_zip_path <- file.path("models/Melanoma_response_train.zip")
          annot_csv_name <- "response_train.csv"
          if (file.exists(annot_zip_path)) {
            # Check what files are in the zip
            zip_files <- unzip(annot_zip_path, list = TRUE)$Name
            if ("Melanoma_response_train.csv" %in% zip_files) {
              annot_csv_name <- "Melanoma_response_train.csv"
            }
          }
          if (!file.exists(annot_zip_path)) stop("Annotation zip file not found: ", annot_zip_path)
          sampleAnnot.train <- read.csv(unz(annot_zip_path, annot_csv_name))
          
          # Map Melanoma annotation column names to match mUC/RCC format
          if ("SampleID" %in% colnames(sampleAnnot.train) && !"RNASEQ_SAMPLE_ID" %in% colnames(sampleAnnot.train)) {
            sampleAnnot.train$RNASEQ_SAMPLE_ID <- sampleAnnot.train$SampleID
          }
          if ("Responder" %in% colnames(sampleAnnot.train) && !"RESPONSE" %in% colnames(sampleAnnot.train)) {
            sampleAnnot.train$RESPONSE <- ifelse(sampleAnnot.train$Responder == "R" | sampleAnnot.train$Responder == 1, 1, 0)
          }
          
          if (!"RNASEQ_SAMPLE_ID" %in% colnames(sampleAnnot.train)) {
            stop("No RNASEQ_SAMPLE_ID column in Melanoma annotations")
          }
          if (!any(sampleIDs_train %in% sampleAnnot.train$RNASEQ_SAMPLE_ID)) {
            stop("No common sample IDs between Melanoma training data and annotations")
          }
          common_samples <- intersect(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID)
          sample_idx <- which(sampleIDs_train %in% common_samples)
          sampleIDs_train <- sampleIDs_train[sample_idx]
          train_data <- train_data_raw_Melanoma[sampleIDs_train, train_cols_order, drop = FALSE]
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
        
        # Final guard: ensure no NA in x.test before prediction
        if (anyNA(x.test)) {
          x.test <- scrub_na(x.test, fill = 0)
          output$status <- renderText("Some NA values detected after preprocessing; replaced with 0 to proceed.")
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
        results_list <- vector("list", nrow(x.test))
        for (i in seq_len(nrow(x.test))) {
          x_sample <- x.test[i, , drop = FALSE]
          avg_distances <- compute_average_cosine_distances(y1, y2, x_sample)
          # Create data.frame with all columns, using list to preserve column name with special character
          results_list[[i]] <- data.frame(
            CosDist_2_Rs          = avg_distances[1],
            CosDist_2_NRs         = avg_distances[2],
            LogitDA_Score         = pred_prob[i],
            LogitDA_pred_label    = pred_labels[i],
            stringsAsFactors = FALSE
          )
          # Add column with proper name (preserving > character)
          results_list[[i]][["LogitDA_score > 0.5"]] <- pred_labels[i]
        }
        results_df <- do.call(rbind, results_list)
        # Ensure column name with special character is preserved
        if ("LogitDA_score > 0.5" %in% names(results_df)) {
          # Column name is already correct
        } else {
          # If somehow lost, restore it
          col_idx <- which(names(results_df) == "LogitDA_pred_label")
          if (length(col_idx) > 0 && col_idx < ncol(results_df)) {
            names(results_df)[col_idx + 1] <- "LogitDA_score > 0.5"
          }
        }
        results_df$sampleID <- sampleIDs
        
        # Add LogitDA_score_label based on LogitDA_Score thresholds
        results_df$LogitDA_score_label <- sapply(results_df$LogitDA_Score, function(s) {
          if (s >= 0.5002) "R" else if (s < 0.2895) "NR" else NA_character_
        })
        
        results_df <- results_df[match(sampleIDs, results_df$sampleID), ]
        results_df <- results_df[, c("sampleID", setdiff(names(results_df), "sampleID"))]
        stopifnot(all(results_df$sampleID == sampleIDs))
        
        # ---------------- ORR prior + final table ----------------
        orr_results <- apply_orr_prior(results_df, orr = 0.2)
        results_df  <- save_and_report_results(
          results_df = orr_results$results_df,
          prior_final_percentage = orr_results$prior_final_percentage
        )
        
        # Remove internal-only prediction label column from final outputs
        if ("LogitDA_pred_label" %in% names(results_df)) {
          results_df$LogitDA_pred_label <- NULL
        }
        
        # Reorder columns: LogitDA_score_label before % of applicability
        col_order <- c("sampleID", "CosDist_2_Rs", "CosDist_2_NRs", "LogitDA_Score", 
                       "LogitDA_score > 0.5", "LogitDA_score_label", "% of applicability", 
                       "iCosinDist_label")
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
                                        orr_results$prior_final_percentage))
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
        keep_cols <- c("sampleID", "LogitDA_Score", "LogitDA_score_label")
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
