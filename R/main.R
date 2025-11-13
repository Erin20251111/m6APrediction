
##' Encode DNA 5-mer strings into factor columns per nucleotide position
##'
##' Converts a vector of equal-length DNA strings (e.g., 5-mers) into a data.frame
##' with one column per nucleotide position. Each column is a factor with levels
##' A, T, C, G, suitable for downstream modeling.
##'
##' @param dna_strings Character vector of DNA strings of identical length
##'   (e.g., 5 for 5-mer). Only nucleotides A, T, C, and G are expected.
##'
##' @return A data.frame where each column corresponds to a nucleotide position
##'   (`nt_pos1`, `nt_pos2`, ...) encoded as factors with levels A, T, C, G.
##'
##' @examples
##' \dontrun{
##'   # Internal helper; not exported.
##'   dna_encoding(c("ATCGA", "TGGCA"))
##' }
##'
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


##' Predict m6A status for multiple observations
##'
##' Takes a fitted machine learning model and a feature data.frame, encodes the
##' DNA 5-mer features, and returns predicted probabilities and statuses for m6A.
##'
##' @param ml_fit A fitted classification model that supports
##'   \code{predict(..., type = "prob")} and returns a column named \code{"Positive"}.
##'   Typically a \code{randomForest} model or a caret-trained model wrapping it.
##' @param feature_df A data.frame containing at least the columns:
##'   \code{gc_content}, \code{RNA_type}, \code{RNA_region}, \code{exon_length},
##'   \code{distance_to_junction}, \code{evolutionary_conservation}, \code{DNA_5mer}.
##'   \code{RNA_type} will be coerced to factor with levels:
##'   \code{c("mRNA","lincRNA","lncRNA","pseudogene")}.
##'   \code{RNA_region} will be coerced to factor with levels:
##'   \code{c("CDS","intron","3'UTR","5'UTR")}.
##' @param positive_threshold Numeric threshold (0-1) above which a sample is
##'   labeled \code{"Positive"}; otherwise \code{"Negative"}. Default is 0.5.
##'
##' @return The input \code{feature_df} with two additional columns:
##'   \itemize{
##'     \item \code{predicted_m6A_prob}: numeric probability for class Positive
##'     \item \code{predicted_m6A_status}: character label \code{"Positive"} or \code{"Negative"}
##'   }
##'
##' @examples
##' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
##' df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
##' out <- prediction_multiple(rf_fit, df, positive_threshold = 0.6)
##' head(out)
##'
##' @import randomForest
##' @importFrom stats predict
##' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  encoded_dna <- dna_encoding(feature_df$DNA_5mer)
  input_df <- cbind(
    feature_df[, c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation")],
    encoded_dna
  )
  input_df$RNA_type <- factor(input_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  input_df$RNA_region <- factor(input_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  #prob <- predict(ml_fit, input_df, type = "prob")[, "Positive"]
  prob <- predict(ml_fit, newdata = input_df, type = "prob")[, "Positive"]
  status <- ifelse(prob >= positive_threshold, "Positive", "Negative")
  feature_df$predicted_m6A_prob <- prob
  feature_df$predicted_m6A_status <- status
  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}


##' Predict m6A status for a single observation
##'
##' Convenience wrapper around \code{prediction_multiple()} for one set of
##' features, returning a named vector with probability and status.
##'
##' @param ml_fit A fitted classification model that supports
##'   \code{predict(..., type = "prob")} and returns a column named \code{"Positive"}.
##' @param gc_content Numeric GC content feature.
##' @param RNA_type Character or factor; one of \code{c("mRNA","lincRNA","lncRNA","pseudogene")}.
##' @param RNA_region Character or factor; one of \code{c("CDS","intron","3'UTR","5'UTR")}.
##' @param exon_length Numeric exon length feature.
##' @param distance_to_junction Numeric distance-to-junction feature.
##' @param evolutionary_conservation Numeric evolutionary conservation feature.
##' @param DNA_5mer Character DNA 5-mer sequence (e.g., \code{"ATCGA"}).
##' @param positive_threshold Numeric threshold (0-1) for Positive vs Negative; default 0.5.
##'
##' @return A named vector of length 2:
##'   \code{c(predicted_m6A_prob = <numeric>, predicted_m6A_status = <character>)}.
##'
##' @examples
##' rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
##' df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
##' r1 <- df[1, ]
##' prediction_single(
##'   ml_fit = rf_fit,
##'   gc_content = r1$gc_content,
##'   RNA_type = as.character(r1$RNA_type),
##'   RNA_region = as.character(r1$RNA_region),
##'   exon_length = r1$exon_length,
##'   distance_to_junction = r1$distance_to_junction,
##'   evolutionary_conservation = r1$evolutionary_conservation,
##'   DNA_5mer = r1$DNA_5mer,
##'   positive_threshold = 0.5
##' )
##'
##' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )
  result_df <- prediction_multiple(ml_fit, single_df, positive_threshold = positive_threshold)
  returned_vector <- c(
    predicted_m6A_prob = result_df$predicted_m6A_prob[1],
    predicted_m6A_status = result_df$predicted_m6A_status[1]
  )
  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}

