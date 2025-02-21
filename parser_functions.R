library(vcfR)
library(tidyverse)
library(tidyr)
library(ggplot2)

parse_snpeff_annotated_sv <- function(ann_manta_VRanges) {
  genes <- c(rep(NA, length(ann_manta_VRanges)))
  multiple_genes <- c(rep(NA, length(ann_manta_VRanges)))
  svtype <- c(rep(NA, length(ann_manta_VRanges)))
  for (i in 1:length(ann_manta_VRanges)) {
    svtype[i] <- ann_manta_VRanges[i,]$SVTYPE
    ann_text <- ann_manta_VRanges[i,]$ANN[[1]]
    if (length(ann_text) == 0) {
      genes[i] <- NA
    } else if (length(ann_text) > 1) {
      gene <- c(rep(NA, length(ann_text)))
      for (j in 1:length(ann_text)) {
        gene[j] <- strsplit(ann_text, "|", fixed = TRUE)[[j]][4]
      }
      multiple_genes[i] <- paste(gene, collapse = "|")
    } else {
      gene <- strsplit(ann_text, "|", fixed = TRUE)[[1]][4]
      genes[i] <- gene
    }
  }
  return(list(genes, multiple_genes, svtype))
}

genes_list_post_processing <- function(sv_ann_df, filter_cutoff = 15) {
  for (i in 1:length(sv_ann_df$multi)) {
    if (!is.na(sv_ann_df$multi[i]) && sv_ann_df$multi[i] == "") {
      sv_ann_df$multi[i] <- NA
    } else {
      uniq_values <- unique(strsplit(sv_ann_df$multi[i], "|", fixed = TRUE)[[1]])
      uniq_values <- gsub("&", "|", uniq_values)
      uniq_multi <- paste(uniq_values,
                          sep = "|",
                          collapse = "|")
      sv_ann_df$multi[i] <- uniq_multi
    }
    if (sv_ann_df$multi[i] == "NA") {
      sv_ann_df$multi[i] <- NA
    }
  }
  
  for (i in 1:length(sv_ann_df$single)) {
    if (!is.na(sv_ann_df$single[i]) && sv_ann_df$single[i] == "") {
      sv_ann_df$single[i] <- NA
    } else {
      sv_ann_df$single[i] <- gsub("&", "|", sv_ann_df$single[i])
    }
  }
  
  # turn it to count/observation data (can have repeats)
  len_gene_column <- length(sv_ann_df$single)
  single_genes <- as.vector(sv_ann_df$single)
  events <- as.vector(sv_ann_df$event)
  
  for (i in 1:len_gene_column) {
    if (!is.na(single_genes[i])) {
      split_str <- strsplit(single_genes[i], "|", fixed = TRUE)[[1]]
      if (length(split_str) > 1) {
        for (j in 2:length(split_str)) {
          single_genes <- append(single_genes, split_str[j])
          events <- append(events, sv_ann_df$event[i])
        }
        single_genes[i] <- split_str[1]
      }
    }
  }
  
  len_gene_column <- length(sv_ann_df$multi)
  multi_genes <- as.vector(sv_ann_df$multi)
  for (i in 1:len_gene_column) {
    if (!is.na(multi_genes[i])) {
      split_str <- strsplit(multi_genes[i], "|", fixed = TRUE)[[1]]
      if (length(split_str) > 1) {
        for (j in 2:length(split_str)) {
          single_genes <- append(single_genes, split_str[j])
          events <- append(events, sv_ann_df$event[i])
        }
        single_genes[i] <- split_str[1]
      }
    }
  }
  
  len_gene_column <- length(single_genes)
  for (i in 1:len_gene_column) {
    if (!is.na(single_genes[i])) {
      split_str <- strsplit(single_genes[i], "-", fixed = TRUE)[[1]]
      if (length(split_str) > 1) {
        for (j in 2:length(split_str)) {
          single_genes <- append(single_genes, split_str[j])
          events <- append(events, events[i])
        }
        single_genes[i] <- split_str[1]
      }
    }
  }

  proccessed_sv_annotations <- data.frame(x=single_genes, y=events)
  colnames(proccessed_sv_annotations) <- c("gene", "event")
  
  proccessed_sv_annotations <- proccessed_sv_annotations[complete.cases(proccessed_sv_annotations), ]
  proccessed_sv_annotations <- proccessed_sv_annotations[grep("^[LOC]{3}", proccessed_sv_annotations$gene, invert = TRUE), ]
  proccessed_sv_annotations$count <- rep(1, nrow(proccessed_sv_annotations))
  
  low_gene_count <- proccessed_sv_annotations %>% 
    group_by(gene) %>%
    count %>%
    filter(n < filter_cutoff) %>% 
    pull(gene)
  
  low_gene_count <- append(low_gene_count, "")
  
  proccessed_sv_annotations <- proccessed_sv_annotations[ ! proccessed_sv_annotations$gene %in% low_gene_count, ]
  
  return(proccessed_sv_annotations)
}

parse_sv <- function(ann_manta_file) {
  ann_manta_VRanges <- VariantAnnotation::readVcfAsVRanges(ann_manta_file)
  print("Read file")
  genes_lists <- parse_snpeff_annotated_sv(ann_manta_VRanges)
  print("Parsed file")
  
  # construct df from genes list
  sv_genes_df <- data.frame(genes_lists)
  colnames(sv_genes_df) <- c("single", "multi", "event")
  
  processed_sv_genes_df <- genes_list_post_processing(sv_genes_df)
  print("Proccessed file")
  
  return(processed_sv_genes_df)
}

batch_parse_sv <- function(list_ann_manta_files) {
  joined_df = NULL
  for (i in 1:length(list_ann_manta_files)) {
    parsed_sv <- parse_sv(list_ann_manta_files[i])
    if (i == 1) {
      joined_df <- parsed_sv
    } else if (i > 1) {
      joined_df <- rbind(joined_df, parsed_sv)
    }
  }
  return(joined_df)
}
