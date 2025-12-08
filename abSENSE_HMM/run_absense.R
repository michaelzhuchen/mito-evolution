suppressMessages(library(tidyverse))
suppressMessages(library(ape))
suppressMessages(library(castor))
suppressMessages(library(ggplot2))
# Set ggplot theme
theme_set(theme_classic())

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
OG_id <- args[1]
hmm_LOO_results_prefix <- args[2]
alerax_dir <- args[3]
posterior_clades_dir <- args[4]
outdir <- args[5]
BOOL_verbose <- FALSE

# Settings
clade_support_threshold <- 50
species_overlap_threshold <- 0.1
clade_purity_threshold <- 0.9
species_overlap_support_threshold <- 0.5
duplications_rec_support_threshold <- 0
n_euk_species_in_dataset <- 203
n_archaea_species_in_dataset <- 337
n_bacteria_species_in_dataset <- 1737
fraction_euk_species_threshold <- round(2 * (n_euk_species_in_dataset) / (n_euk_species_in_dataset + n_archaea_species_in_dataset + n_bacteria_species_in_dataset), digits=2)

bitscore_at_expect_threshold <- 22.95 # corresponding to e-value = 0.001

if (BOOL_verbose) {
  options(warn = 1)
}

## Read in data
if (BOOL_verbose) {
  print("Reading in data")
}


hmmsearch_out_raw <- read.table(paste0(hmm_LOO_results_prefix, "/results/hmmsearch_out/", OG_id, "_LOO_hmmsearch.out"))
colnames(hmmsearch_out_raw)[1] <- "accession"
colnames(hmmsearch_out_raw)[3] <- "target_length"
colnames(hmmsearch_out_raw)[6] <- "hmm_length"
colnames(hmmsearch_out_raw)[7] <- "expect"
colnames(hmmsearch_out_raw)[8] <- "bitscore"
colnames(hmmsearch_out_raw)[13] <- "i_evalue"
hmmsearch_out_raw$OG_id <- gsub("_LOO.*", "", hmmsearch_out_raw$V4)
hmmsearch_out_raw <- hmmsearch_out_raw[order(hmmsearch_out_raw$bitscore, decreasing=TRUE),]

# Read in taxonomic data
uniprot_proteomes_tax <- read.table(here("data/taxonomy", "uniprot_new.eukaryota_prokgroups_other.opisthokonta_parasitic.plants_BaSk_CRuMs_downsample_combined_ncbi_taxonomy.tsv"), sep="\t", header=TRUE)

species_tree <- read.tree(here("data/species_phylogeny/processed_species_tree", "concat_cytosolic_ribosomal_proteins_97.5pct.spp_muscle5_clipkit.gappy.msa_constrained.ncbi.tree.manual.changes.v7_prokspp.collapsed_nodelabels_rooted_downsample_v2.contree"))
species_tree_prok <- drop.tip(species_tree, uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Eukaryota"])
species_tree_prok_labels <- c(species_tree_prok$tip.label, species_tree_prok$node.label)
species_tree_labels <- c(species_tree$tip.label, species_tree$node.label)
species_tree_euks_subtree <- get_subtree_at_node(species_tree, "Node34_Eukaryota")$subtree
species_tree_euks_labels <- c(species_tree_euks_subtree$tip.label, species_tree_euks_subtree$node.label)
prok_euk_parent_labels <- c("Node1_cellular_organisms", "Node30_Archaea", "Node31_Archaea", "Node32_Archaea", "Node33_Asgardgroup")

# Calculate node depths
depths_df <- data.frame(label = c(species_tree$tip.label, species_tree$node.label), distance=get_all_distances_to_root(species_tree, as_edge_count=TRUE))

## Run absense
if (BOOL_verbose) {
  print("Running absense")
}


homology_power_agg <- c()
r_squared_agg <- c()
residuals_agg <- c()
homology_power_per_species <- c()

query_OG_list <- OG_id

for (curr_index in 1:length(query_OG_list)) {
  query_OG <- query_OG_list[curr_index]
  if (BOOL_verbose) {
    print(paste0(query_OG, " ", curr_index, "/", length(query_OG_list)))
  }
  
  hmmsearch_out <- hmmsearch_out_raw[which(hmmsearch_out_raw$OG_id == query_OG),]
  hmmsearch_out <- hmmsearch_out %>% filter(!duplicated(accession))
  hmmsearch_out$tree_id <- gsub("_.*", "", hmmsearch_out$accession)
  hmm_length <- max(hmmsearch_out$hmm_length)
  
  # Select hmmsearch hits that are OG members
  hmmsearch_out$OG_member <- TRUE
  hmmsearch_out_OG_members <- hmmsearch_out[hmmsearch_out$OG_member,]
  hmmsearch_out <- hmmsearch_out[hmmsearch_out$OG_member,]
  
  curr_OG_taxids <- unique(hmmsearch_out$tree_id)
  human_accessions_collapse <- paste0(grep("9606_", hmmsearch_out$accession, value=TRUE), collapse=",")
  
  # Get top hit by bitscore
  hmmsearch_out <- hmmsearch_out %>% group_by(tree_id) %>% filter(bitscore == max(bitscore)) %>% filter(!duplicated(tree_id))
  
  inferred_species_tree_rooted <- species_tree
  
  # Exclude prokaryotes
  inferred_species_tree_rooted <- keep.tip(inferred_species_tree_rooted, uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Eukaryota"])
  hmmsearch_out <- hmmsearch_out[hmmsearch_out$tree_id %in% uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Eukaryota"],]
  
  # Find OG and outgroup MRCAs
  mrca_index <- get_mrca_of_set(inferred_species_tree_rooted, curr_OG_taxids)
  if (mrca_index > Ntip(inferred_species_tree_rooted)) {
    mrca_label <- inferred_species_tree_rooted$node.label[mrca_index - Ntip(inferred_species_tree_rooted)]
    OG_members_subtree <- get_subtree_at_node(inferred_species_tree_rooted, mrca_index - Ntip(inferred_species_tree_rooted))$subtree
  } else {
    mrca_label <- inferred_species_tree_rooted$tip.label[mrca_index]
    OG_members_subtree <- keep.tip(inferred_species_tree_rooted, mrca_label)
  }
  
  # If no hits to self OG
  if (nrow(hmmsearch_out) == 0) {
    homology_power <- data.frame(OG_id = query_OG, human_accessions = human_accessions_collapse, L = NA, R = NA, OG_mrca_label = mrca_label, OG_outgroup_mrca_label = "", outgroup_kingdom_species = "", fraction_of_detectable_outgroup_kingdoms = NA, probability_of_detection_in_any_outgroup_species = NA, hmm_length = hmm_length, r_squared = NA, comments = "No hmm hits to self OG")
    homology_power_agg <- rbind(homology_power_agg, homology_power)
    write.table(homology_power_agg, paste0(outdir, "/results/", OG_id, "_homology_power_result.tsv"), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
    quit(save = "no")
  }
  
  dist_mat <- get_all_pairwise_distances(inferred_species_tree_rooted)
  
  # # Using distance to OG MRCA
  distance_to_mrca <- data.frame(tiplabel = c(OG_members_subtree$tip.label, OG_members_subtree$node.label), distance = get_all_distances_to_root(OG_members_subtree))
  hmmsearch_out$distance_to_mrca <- distance_to_mrca$distance[match(hmmsearch_out$tree_id, distance_to_mrca$tiplabel)]
  
  # Edge case: If fewer than 3 OG member species proteins found by hmmsearch, skip since cannot fit function
  if (length(unique(hmmsearch_out_OG_members$tree_id)) < 3) {
    homology_power <- data.frame(OG_id = query_OG, human_accessions = human_accessions_collapse, L = NA, R = NA, OG_mrca_label = mrca_label, OG_outgroup_mrca_label = "", outgroup_kingdom_species = "", fraction_of_detectable_outgroup_kingdoms = NA, probability_of_detection_in_any_outgroup_species = NA, hmm_length = hmm_length, r_squared = NA, comments = "fewer than 3 OG member species hit")
    homology_power_agg <- rbind(homology_power_agg, homology_power)
    write.table(homology_power_agg, paste0(outdir, "/results/", OG_id, "_homology_power_result.tsv"), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
    quit(save = "no")
  }
  
  OG_members_subtree_only <- keep.tip(inferred_species_tree_rooted, hmmsearch_out$tree_id)
  
  dist_mat <- get_all_pairwise_distances(OG_members_subtree_only)
  
  # Match rows and columns of distance matrices to proteins
  dist_mat_tips <- matrix(data = dist_mat[,match(hmmsearch_out$tree_id, OG_members_subtree_only$tip.label)], nrow = nrow(dist_mat), ncol = length(hmmsearch_out$tree_id))
  dist_mat_tips <- matrix(data = dist_mat_tips[match(hmmsearch_out$tree_id, OG_members_subtree_only$tip.label),], nrow = length(hmmsearch_out$tree_id), ncol = length(hmmsearch_out$tree_id))
  
  hmmsearch_out$pairwise_distance <- rowMeans(dist_mat_tips, na.rm=TRUE)
  
  # Edge case: If fewer than 3 unique OG member species distances, skip since cannot fit function reliably
  if (length(unique(hmmsearch_out$pairwise_distance)) < 3) {
    homology_power <- data.frame(OG_id = query_OG, human_accessions = human_accessions_collapse, L = NA, R = NA, OG_mrca_label = mrca_label, OG_outgroup_mrca_label = "", outgroup_kingdom_species = "", fraction_of_detectable_outgroup_kingdoms = NA, probability_of_detection_in_any_outgroup_species = NA, hmm_length = hmm_length, r_squared = NA, comments = "fewer than 3 unique weighted distances")
    homology_power_agg <- rbind(homology_power_agg, homology_power)
    write.table(homology_power_agg, paste0(outdir, "/results/", OG_id, "_homology_power_result.tsv"), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
    quit(save = "no")
  }
  
  # Get distance for all species
  core_species <- c("9606", "559292", "1257118", "5741", "185431", "3702", "36329")
  species_distance_df <- data.frame(tree_id = uniprot_proteomes_tax$tree_id[uniprot_proteomes_tax$domain == "Eukaryota"], species_name = uniprot_proteomes_tax$ScientificName[uniprot_proteomes_tax$domain == "Eukaryota"], weighted_distance = NA)
  
  representative_species <- c("2769", "280699", "4787", "2903", "3055", "5762", "44689", "37360", "412133", "2049356", "461836", "905079", "1746090", "2654837", "81824", "595528", "667725", "691883", "2983909") # kingdom level
  
  selected_species <- c(core_species, representative_species)
  
  # Specify desired distance metric
  hmmsearch_out$distance_metric <- hmmsearch_out$pairwise_distance
  
  # Pairwise distance.
  species_distance_df$distance_metric <- NA
  for (i in 1:nrow(species_distance_df)) {
    curr_tree_id <- species_distance_df$tree_id[i]
    curr_subtree <- keep.tip(inferred_species_tree_rooted, c(curr_tree_id, unique(hmmsearch_out$tree_id)))
    distances_to_tip <- get_all_distances_to_tip(curr_subtree, curr_tree_id)
    
    distances_to_tip_by_protein <- distances_to_tip[match(hmmsearch_out$tree_id, curr_subtree$tip.label)]
    species_distance_df$distance_metric[i] <- mean(distances_to_tip_by_protein)
  }
  species_distance_df$distance_metric[species_distance_df$tree_id %in% hmmsearch_out$tree_id] <- hmmsearch_out$distance_metric[match(species_distance_df$tree_id[species_distance_df$tree_id %in% hmmsearch_out$tree_id], hmmsearch_out$tree_id)] # replace OG members with their distances
  species_distance_df_select <- species_distance_df[which(species_distance_df$tree_id %in% core_species | species_distance_df$distance_metric == max(species_distance_df$distance_metric)),] # Get the set of selected core species and the maximal distant species
  species_distance_df_select$OG_member <- species_distance_df_select$tree_id %in% hmmsearch_out_OG_members$tree_id
  species_distance_df_select$OG_member_color <- "red"
  species_distance_df_select$OG_member_color[species_distance_df_select$OG_member] <- "black"
  
  # Fit to score function
  nlc <- nls.control(maxiter = 1000)
  starting_params_list <- list(
    list(L = 500, R = 1),
    list(L = 100, R = 0),
    list(L =  20, R = 0),
    list(L   = 10, R = 0),
    list(L   = 0, R = 0),
    list(L   = 0, R = 1)
  )
  fit_with_starting_params <- function(formula, data, starting_params_list, control) {
    for (st in starting_params_list) {
      fit <- try(nls(formula, data = data, start = st, control = control), silent = TRUE)
      if (!inherits(fit, "try-error")) return(fit)
    }
    stop(paste0("Model failed for all starting guesses for ", OG_id, "_", query_OG))
  }
  custom_model <- NULL
  custom_model <- fit_with_starting_params(
    bitscore ~ L * exp(-R * distance_metric),
    data       = hmmsearch_out,
    starting_params_list = starting_params_list,
    control    = nlc
  )
  
  # Extract parameters
  curr_params <- custom_model$m$getPars()
  L <- curr_params["L"]
  R <- curr_params["R"]
  
  # Model fit
  r_squared <- cor(hmmsearch_out$distance_metric, log(hmmsearch_out$bitscore), method="pearson")^2
  
  if (BOOL_verbose) {
    print(paste0("L: ", L))
    print(paste0("R: ", R))
    print(paste0("R^2: ", r_squared))
  }
  
  if (L < 0 | R < 0) {
    print("Negative parameter values, poor fit to model")
    homology_power <- data.frame(OG_id = query_OG, human_accessions = human_accessions_collapse, L = L, R = R, OG_mrca_label = mrca_label, OG_outgroup_mrca_label = "", outgroup_kingdom_species = NA, fraction_of_detectable_outgroup_kingdoms = NA, probability_of_detection_in_any_outgroup_species = NA, hmm_length = hmm_length, r_squared = r_squared, comments = "negative model parameters")
    homology_power_agg <- rbind(homology_power_agg, homology_power)
    write.table(homology_power_agg, paste0(outdir, "/results/", OG_id, "_homology_power_result.tsv"), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
    quit(save = "no")
  }
  
  curr_params_df <- data.frame(OG_id = query_OG, L, R)
  
  # Generate prediction curve
  x_seq <- seq(min(hmmsearch_out$distance_metric), max(species_distance_df$distance_metric), length.out = 100)
  predicted <- data.frame(x = x_seq, y = predict(custom_model, newdata = data.frame(distance_metric = x_seq)))
  
  # Compute prediction interval based on parameter and residual variance
  # Define prediction function
  predict_function <- function(L, R, x) {
    L * exp(-R * x)
  }
  predict_interval <- function(x_val) {
    # Prediction
    y_pred <- predict_function(params["L"], params["R"], x_val)
    
    # Gradient of the function
    grad <- c(
      exp(-params["R"] * x_val),
      -params["L"] * x_val * exp(-params["R"] * x_val)
    )
    
    # Variance of prediction
    var_pred <- t(grad) %*% vcov_matrix %*% grad
    
    # Full prediction interval
    var_pred + residual_variance
    
    interval_width <- t_value * sqrt(var_pred + residual_variance)
    c(variance = var_pred + residual_variance, lower = y_pred - interval_width, upper = y_pred + interval_width)
  }
  # Residual variance (sigma^2)
  residual_variance <- sum(resid(custom_model)^2) / df.residual(custom_model)
  sigma <- sqrt(residual_variance)  # Standard deviation of residuals
  # Get parameter estimates and covariance matrix
  params <- coef(custom_model)
  vcov_matrix <- vcov(custom_model)
  # Calculate prediction intervals
  t_value <- qt(0.975, df = df.residual(custom_model))  # 95% prediction interval
  # Apply the function across a range of x values
  x_seq <- seq(min(hmmsearch_out$distance_metric), max(species_distance_df$distance_metric), length.out = 100)
  prediction_intervals <- t(sapply(x_seq, predict_interval))
  predicted$variance <- prediction_intervals[,"variance"]
  predicted$y_lower_99_prediction_interval <- prediction_intervals[,"lower"]
  predicted$y_upper_99_prediction_interval <- prediction_intervals[,"upper"]
  predicted_variance <- t(sapply(species_distance_df$distance_metric, predict_interval))[,"variance"]
  
  ## Compute probability of detecting a homolog in target species
  # Predict score and variance for target species
  predicted_score_per_species <- data.frame(x = species_distance_df$distance_metric, y = predict(custom_model, newdata = data.frame(distance_metric = species_distance_df$distance_metric)))
  species_predicted_score <- cbind(species_distance_df, predicted_score=predicted_score_per_species$y)
  species_predicted_score$variance <- predicted_variance
  # Assuming new sample is normally distributed, compute the CDF at score threshold
  species_predicted_score$probability_of_detection <- 1 - pnorm(bitscore_at_expect_threshold, mean = species_predicted_score$predicted_score, sd = sqrt(species_predicted_score$variance), lower.tail = TRUE)
  
  species_predicted_score_tax <- merge(species_predicted_score, uniprot_proteomes_tax, by.x="tree_id", by.y="tree_id")
  species_predicted_score_tax <- species_predicted_score_tax[order(species_predicted_score_tax$tree_id),]
  species_predicted_score_tax_nonOGmember <- species_predicted_score_tax[which(!species_predicted_score_tax$tree_id %in% hmmsearch_out_OG_members$tree_id),]
  
  # Power to detect homologs
  species_predicted_score_tax_nonOGmember_outgroup <- species_predicted_score_tax_nonOGmember %>% filter(!tree_id %in% OG_members_subtree$tip.label)
  species_predicted_score_tax_nonOGmember_outgroup <- species_predicted_score_tax_nonOGmember_outgroup %>% filter(tree_id %in% selected_species) %>% mutate(probability_undetected = 1 - probability_of_detection)
  
  if (nrow(species_predicted_score_tax_nonOGmember_outgroup) > 0) {
    fraction_of_detectable_outgroup_kingdoms <- sum(species_predicted_score_tax_nonOGmember_outgroup$probability_of_detection > 0.95) / nrow(species_predicted_score_tax_nonOGmember_outgroup)
    probability_of_detection_in_any_outgroup_species <- 1 - prod(species_predicted_score_tax_nonOGmember_outgroup$probability_undetected)
    
    outgroup_mrca_index <- get_mrca_of_set(inferred_species_tree_rooted, species_predicted_score_tax_nonOGmember_outgroup$tree_id)
    if (outgroup_mrca_index > Ntip(inferred_species_tree_rooted)) {
      outgroup_mrca_label <- inferred_species_tree_rooted$node.label[outgroup_mrca_index - Ntip(inferred_species_tree_rooted)]
      outgroup_mrca_label_shorten <- gsub("Node[0-9]_", "", outgroup_mrca_label)
    } else {
      outgroup_mrca_label <- inferred_species_tree_rooted$tip.label[outgroup_mrca_index]
    }
    
  } else {
    # If there are no eukaryotic outgroups, have full coverage of Eukaryota
    fraction_of_detectable_outgroup_kingdoms <- NA
    probability_of_detection_in_any_outgroup_species <- NA
    outgroup_mrca_label <- ""
  }
  
  homology_power_per_species_curr <- data.frame(OG_id = OG_id, tree_id = species_predicted_score_tax$tree_id, probability_of_detection = round(species_predicted_score_tax$probability_of_detection, 4))
  homology_power_per_species <- rbind(homology_power_per_species, homology_power_per_species_curr)
  
  homology_power <- data.frame(OG_id = query_OG, human_accessions = human_accessions_collapse, L = L, R = R, OG_mrca_label = mrca_label, OG_outgroup_mrca_label = outgroup_mrca_label, outgroup_kingdom_species = paste0(species_predicted_score_tax_nonOGmember_outgroup$tree_id, collapse=","), fraction_of_detectable_outgroup_kingdoms = fraction_of_detectable_outgroup_kingdoms, probability_of_detection_in_any_outgroup_species = probability_of_detection_in_any_outgroup_species, hmm_length = hmm_length, r_squared = r_squared, comments = "")
  homology_power_agg <- rbind(homology_power_agg, homology_power)
}

# ## Write out
# write.table(homology_power_agg, paste0(outdir, "/results/", OG_id, "_homology_power_result.tsv"), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
# write.table(homology_power_per_species, paste0(outdir, "/results_per_species/", OG_id, "_homology_power_per_species.tsv"), sep="\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
# 
# ## Plot power per OG
# plot_filename <- paste0(outdir, "/figures/", query_OG)
# 
# # Add lineage information
# hmmsearch_out$Lineage <- uniprot_proteomes_tax$superfamily[match(hmmsearch_out$tree_id, uniprot_proteomes_tax$tree_id)]
# 
# pdf(paste0(plot_filename, ".pdf"), width=8, height=6)
# p <- ggplot(data=hmmsearch_out, aes(x=distance_metric, y=bitscore, color=Lineage)) + geom_point(alpha=0.5) + geom_ribbon(data = predicted, aes(x=x, y=y, ymin=y_lower_99_prediction_interval, ymax=y_upper_99_prediction_interval), colour="#45b6fe", fill="lightblue", alpha=0.5) + geom_line(data = predicted, aes(x = x, y = y), color = "red", linewidth = 1) + geom_hline(aes(yintercept = bitscore_at_expect_threshold), linetype="dashed")
# p <- p + scale_x_continuous(
#   breaks = species_distance_df_select$distance_metric,
#   labels = paste0(gsub(" .*", "", species_distance_df_select$species_name), " (", round(species_distance_df_select$distance_metric,1), ")"),
# ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color=species_distance_df_select$OG_member_color)) + xlab("Weighted evolutionary distance (substitutions/site)") + ylab("HMM score") #+ ylim(min(c(0,predicted$y_lower_99_prediction_interval)),max(predicted$y_upper_99_prediction_interval))
# p + ggtitle(label = paste0(query_OG, ": L=", round(L,1), ", R=", round(R,2)), subtitle=paste0("Fraction of detectable outgroup kingdoms = ", round(fraction_of_detectable_outgroup_kingdoms,2), ", R^2 = ", round(r_squared,2)))
# dev.off()


