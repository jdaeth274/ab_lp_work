## Functions for Bootstrapping analysis 

recombination_gff_cleaner <- function(in_gff_pmen3){
  #browser()
  ## Takes delim file read in and cleans it up 
  pmen3_reccy_data <- data.frame(data = (matrix(data = NA, nrow = nrow(in_gff_pmen3),
                                                ncol = 9)))
  colnames(pmen3_reccy_data) <- c("isolate","number_isolates", "snp_count",
                                  "length","density", "start", "end", "start_node",
                                  "end_node")
  
  for( k in 1:nrow(in_gff_pmen3)){
    #browser()
    current_row <- in_gff_pmen3[k,]
    current_atty <- current_row$attribute
    current_atty_taxa_snp <- str_split_fixed(current_atty,"taxa=",2)[,2]
    taxa <- str_split_fixed(current_atty_taxa_snp,";snp_count=",2)[,1]
    nodes <- str_split_fixed(sub("node=","",str_split_fixed(current_atty, ";neg",2)[,1]), "->", 2)
    taxa_split <- str_split_fixed(taxa, "\\s",10)
    taxa_first <- taxa_split[1,which(taxa_split[1,] != "")[1]]
    #taxa_first <- sub("\\s*","",taxa_first)
    snp_count <- str_split_fixed(current_atty, ";snp_count=",2)[,2]
    snp_count <- as.numeric(sub(";","",snp_count))
    num_taxa <- str_count(taxa, pattern = "GCA") # This has to change wrt to dataset (can't think of a way round atm)
    length_of_event <- current_row$end - current_row$start + 1
    snp_density <- snp_count / length_of_event
    
    pmen3_reccy_data$isolate[k] <- taxa_first
    pmen3_reccy_data$number_isolates[k] <- num_taxa
    pmen3_reccy_data$snp_count[k] <- snp_count
    pmen3_reccy_data$length[k] <- length_of_event
    pmen3_reccy_data$density[k] <- snp_density
    pmen3_reccy_data$start[k] <- current_row$start
    pmen3_reccy_data$end[k] <- current_row$end 
    pmen3_reccy_data$start_node[k] <- nodes[1,1]
    pmen3_reccy_data$end_node[k] <- nodes[1,2]
    
  }
  return(pmen3_reccy_data)
}

delim_reader <- function(delim_file){
  #browser()
  ## read in the gff recombination events file from gubbins
  reccy_csv <- "No"
  try({
    reccy_csv <- read.delim(delim_file,header = FALSE, comment.char = "#")
    colnames(reccy_csv) <- c("type","prog","class","start","end","trent","alexander","arnold","attribute")
  }, silent = TRUE)
  return(reccy_csv)
}

range_func <- function(x, ranges){
  return(any((apply(ranges, 1, spatstat.utils::inside.range, x = x))))
}
range_func2 <- function(x, ranges){
  return(spatstat.utils::inside.range(x = x, r = ranges))
}

reccy_sum_up <- function(micro_csv, gubb_base, phage_starts = NULL, phage_ends = NULL,
                         clade_name = 2, phage_exclude = TRUE, branch_recons = NULL,
                         max_rec_length = NULL, classified_snps = NULL, MGE_state = "No"){
  
  rocrp_csv <- read.delim(micro_csv,
                        stringsAsFactors = FALSE) %>% 
    rename(id = 1, state = 2)
  cat("Loaded MGE data", "\n")
  ## Get the total split of complete and broken isos, see the number of branches we're working with 
  
  rocrp_isos_clade <- rocrp_csv %>% filter(state == MGE_state) %>% dplyr::select(id)
  
  ## Get the gubbins files
  gubb_pattern <- basename(gubb_base)
  dir_loc <- sub("\\/(?!.*\\/).*$","",gubb_base, perl = TRUE)
  total_gub_files <- list.files(path = dir_loc, pattern = gubb_pattern, full.names = TRUE)
  recombination_events <- total_gub_files[grep("recombination_predictions.gff", total_gub_files)]
  tree_file <- total_gub_files[grep("node_labelled.final_tree.tre", total_gub_files)]
  gub_tree <- read.tree(tree_file)
  cat("Loaded Gubbins data", "\n")
  
  # Get the mutation recon
  #mutty_recon <- read.csv(mutation_recon, stringsAsFactors = FALSE)
  
  clade_gff <- delim_reader(recombination_events)
  clade_gff_table <- recombination_gff_cleaner(clade_gff) %>% mutate(Clade = clade_name)
  cat("Cleaned Recombination predictions", "\n")
  #classified_snps <- typing_gubbins_rec(clade_gff_table, mutty_recon)
  
  if(is.null(branch_recons)){
    clade_terminal <- clade_gff_table %>% filter(number_isolates == 1) %>% rowwise() %>%
      mutate(rocrp = ifelse(isolate %in% rocrp_isos_clade$id, "rocRp", "none")) %>% as.data.frame() %>%
      mutate(Clade = clade_name)
    
    num_isos <- length(gub_tree$tip.label)
    terminal_lengths <- cbind.data.frame(gub_tree$edge.length[sapply(1:num_isos,function(x,y)   which(y==x),y=gub_tree$edge[,2])],gub_tree$tip.label)
    colnames(terminal_lengths) <- c("length","id")
    terminal_lengths <- terminal_lengths %>% mutate(rocrp = ifelse(id %in% rocrp_isos_clade$id, "rocRp","none")) %>%
      as.data.frame() %>%
      rename(end_name = id)
    num_branches <- dplyr::count(terminal_lengths, rocrp)
    
    terminal_grouped <- terminal_lengths %>% group_by(rocrp) %>% summarise(tot_length = sum(length)) %>%
      as.data.frame()
    
  }else{
    recon_table <- read.table(file = branch_recons, sep = "\t", header = TRUE) %>%
      mutate(state = ifelse(state == MGE_state, "broken","wt"))
    clade_terminal <- clade_gff_table %>% left_join(recon_table, by = c("end_node" = "node")) %>%
      rename(comM = state) %>% as.data.frame() %>%
      mutate(Clade = clade_name) 
    snp_terminal <- classified_snps %>% left_join(recon_table, by = c("end_node" = "node")) %>%
      rename(comM = state) %>% as.data.frame() %>%
      mutate(Clade = clade_name) 
    tip_n_num <- c(gub_tree$tip.label,gub_tree$node.label )
    lengths_df <- gub_tree$edge %>% as.data.frame() %>% mutate(length = gub_tree$edge.length) %>% 
      rename(start_node = 1, end_node = 2) %>%
      rowwise() %>%
      mutate(end_name = tip_n_num[end_node])
    terminal_lengths <- lengths_df %>% left_join(recon_table, by = c("end_name" = "node")) %>%
      as.data.frame() %>% rename(comM = state)
    terminal_grouped <- terminal_lengths %>% group_by(comM) %>% summarise(tot_length = sum(length)) %>%
      as.data.frame()
    num_branches <- dplyr::count(terminal_lengths, comM)
    cat("Loading up branch reconstruction data", "\n")
  }
  
  cat("Removing phage regions", "\n")
  if(length(phage_starts) > 0){
    insert_regions <- matrix(data = 0, nrow = length(phage_starts), ncol = 2)
    insert_regions[,1] <- phage_starts
    insert_regions[,2] <- phage_ends
  }
  
  if(phage_exclude){
    for(k in 1:nrow(insert_regions)){
      if(k == 1){
        clade_sans_mge <- clade_terminal
      }
      clade_sans_mge <- clade_sans_mge %>% rowwise() %>% filter(!(range_func2(insert_regions[k,1], c(start, end))) &
                                                                  !range_func2(insert_regions[k,2], c(start, end)) &
                                                                  !range_func2(start, insert_regions[k,]))
    }
    snp_terminal <- snp_terminal %>% rowwise() %>% filter(!(range_func(base_number, insert_regions))) 
  }else{
    for(k in 1:nrow(insert_regions)){
      if(k == 1){
        clade_sans_mge <- clade_terminal
      }
      clade_sans_mge <- clade_sans_mge %>% rowwise() %>% filter((range_func2(insert_regions[k,1], c(start, end))) |
                                                                  range_func2(insert_regions[k,2], c(start, end)))
    }
    snp_terminal <- snp_terminal %>% rowwise() %>% filter((range_func(base_number, insert_regions)))
  }
  
  # clade_terminal_sans <- clade_sans_mge %>% filter(number_isolates == 1) %>% rowwise() %>%
  #   mutate(comM = ifelse(isolate %in% comM_isos_clade$id, "Complete", "Broken")) %>% as.data.frame() %>%
  #   mutate(Clade = clade_name)
  
  if(!is.null(max_rec_length)){
    clade_sans_mge <- clade_sans_mge %>% filter(length <= max_rec_length)
  }
  cat("Creating plots", "\n")
  snail_plot <- ggplot(data = clade_sans_mge) +
    geom_point(aes(x = density, y = length, group = comM, colour = comM)) +
    labs(x = "SNP density", y = "Length", colour = "comM") + scale_x_log10() +
    scale_y_log10()
  
  density_box <- ggplot(data = clade_sans_mge) + geom_boxplot(aes(x = comM, y = density,
                                                                  group = comM, colour = comM)) +
    labs(x = "comM", y = "SNP density", colour = "comM") + scale_y_log10()
  length_box <- ggplot(data = clade_sans_mge) + geom_boxplot(aes(x = comM, y = length,
                                                                 group = comM, colour = comM)) +
    labs(x = "comM", y = "Length", colour = "comM") + scale_y_log10()
  
  df_2 <- clade_sans_mge %>% dplyr::select(-comM)
  length_histos <- ggplot(data = clade_sans_mge, aes(x = length)) + 
    geom_histogram(data = df_2, aes(x = length), fill = "grey70", colour = "blue") +
    geom_histogram(aes(fill = comM),colour = "black") +
    facet_wrap( ~ comM, ncol = 1) + scale_y_log10() + scale_x_log10() +
    labs( x = "Length (bp)", y = "count", title = clade_name)
  
  length_histos <- ggplot(data = clade_sans_mge, aes(x = length)) + 
    geom_histogram(data = df_2, aes(x = length), fill = "grey70", colour = "blue") +
    geom_histogram(aes(fill = comM),colour = "black") +
    facet_wrap( ~ comM, ncol = 1) + scale_y_log10() + scale_x_log10() +
    labs( x = "Length (bp)", y = "count", title = clade_name)
  
  dens_histos <- ggplot(data = clade_sans_mge, aes(x = density)) + 
    geom_histogram(data = df_2, aes(x = density), fill = "grey70", colour = "blue") +
    geom_histogram(colour = "black", aes(fill = comM)) +
    facet_wrap( ~ comM, ncol = 1) + scale_y_log10() + scale_x_log10() +
    labs( x = "Density (SNP/bp)", y = "count", title = clade_name)
  
  
  samp_szes <- dplyr::count(clade_sans_mge, comM)
  out_df <- data.frame(matrix(data = NA, nrow = 1, ncol = 21)) %>% 
    rename(comM_length = 1, wt_length = 2, p_val_length = 3,
           comM_dens = 4, wt_dens = 5, p_val_dens = 6,
           sample_comM = 7, sample_wt = 8, clade = 9, 
           branch_comM = 10, branch_wt = 11, comM_rec_per_length = 12,
           wt_rec_per_length = 13, comM_rm = 14, wt_rm = 15, 
           comM_m = 16, comM_r = 17,
           wt_m = 18, wt_r = 19,
           comM_rho_m = 20, wt_rho_m = 21)
  
  
  clade_terminal_sans <- clade_sans_mge %>% 
    left_join(terminal_lengths %>% rename(branch_length = length) %>%
                dplyr::select(end_name, branch_length), by = c("end_node" = "end_name")) %>%
    as.data.frame()
  
  cat("Creating SNP databases", "\n")
  
  #rownames(out_df) <- c("Length", "Density")
  
  complete_median_dens <- NULL
  snp_branches <- snp_terminal %>% 
    mutate(s_num = ifelse(Type == "S", 1, 0)) %>%
    mutate(r_num = ifelse(Type == "r",1,0)) %>%
    mutate(branch = paste(start_node, end_node, sep = "-")) %>%
    group_by(branch) %>%
    summarise(S_tot = sum(s_num), r_tot = sum(r_num), comM = first(comM)) %>%
    as.data.frame() %>% 
    mutate(start_node = sub("-.*$", "", branch, perl = TRUE)) %>% 
    mutate(end_node = sub("^.*-", "", branch, perl = TRUE)) %>%
    select(-branch)
  
  rec_branches <- clade_terminal_sans %>% 
    mutate(rec_num = 1) %>%
    mutate(branch = paste(start_node, end_node, sep = "-")) %>%
    group_by(branch) %>%
    summarise(rec_tot = sum(rec_num), rec_snps = sum(snp_count), rec_length = sum(length),
              comM = first(comM)) %>%
    as.data.frame() %>%
    mutate(start_node = sub("-.*$", "", branch, perl = TRUE)) %>% 
    mutate(end_node = sub("^.*-", "", branch, perl = TRUE)) %>%
    select(-branch)
  
  by_branch_snps <- snp_branches %>% 
    left_join(rec_branches, by = c("start_node", "end_node", "comM")) %>%
    as.data.frame() %>%
    select(-start_node)
  by_branch_df <- terminal_lengths %>%
    left_join(by_branch_snps , by = c("end_name" = "end_node", "comM" = "comM")) %>%
    mutate(start_node = tip_n_num[start_node]) %>%
    select(-end_node) %>%
    rename(end_node = end_name) %>%
    mutate(clade = clade_name)
  by_branch_df[is.na(by_branch_df)] <- 0
  cat("Running statistical tests", "\n")
  try({
    wc_leng <- wilcox.test(clade_terminal_sans[clade_terminal_sans$comM == "broken", "length"], clade_terminal_sans[clade_terminal_sans$comM == "wt", "length"])
    broke_median_leng <- median(clade_terminal_sans[clade_terminal_sans$comM == "broken", "length"])
    complete_median_leng <- median(clade_terminal_sans[clade_terminal_sans$comM == "wt", "length"])
    comM_r_num <- by_branch_df %>% filter(comM == "broken")
    comM_r_num <- sum(comM_r_num$rec_snps)
    comM_m_num <- snp_terminal %>% filter(comM == "broken") %>% 
      filter(Type == "S") %>% nrow()
    comM_rm <- comM_r_num / comM_m_num
    comM_rho_m <- samp_szes[1,2] / comM_m_num
    wt_r_num <- by_branch_df %>% filter(comM == "wt")
    wt_r_num <- sum(wt_r_num$rec_snps)
    wt_m_num <- snp_terminal %>% filter(comM == "wt") %>% 
      filter(Type == "S") %>% nrow()
    wt_rm <- wt_r_num / wt_m_num
    wt_rho_m <- samp_szes[2,2] / wt_m_num
    
    wc_dens <- wilcox.test(clade_terminal_sans[clade_terminal_sans$comM == "broken", "density"], clade_terminal_sans[clade_terminal_sans$comM == "wt", "density"])
    broke_median_dens <- median(clade_terminal_sans[clade_terminal_sans$comM == "broken", "density"])
    complete_median_dens <- median(clade_terminal_sans[clade_terminal_sans$comM == "wt", "density"])
  }, silent = TRUE)
  if(!is.null(wc_leng)){
    out_df[1,] <- c(broke_median_leng, complete_median_leng, wc_leng$p.value,
                    broke_median_dens, complete_median_dens, wc_dens$p.value,
                    samp_szes[1,2], samp_szes[2,2], clade_name, num_branches[1,2], num_branches[2,2],
                    samp_szes[1,2] / terminal_grouped[1,2], samp_szes[2,2] / terminal_grouped[2,2],
                    comM_rm, wt_rm, comM_m_num, comM_r_num, wt_m_num, wt_r_num, 
                    comM_rho_m, wt_rho_m)
    # out_df[2,] <- c(broke_median_dens, complete_median_dens, wc_dens$p.value, samp_szes[1,2],
    #                 samp_szes[2,2], clade_name, num_branches[1,2], num_branches[2,2],
    #                 samp_szes[1,2] / terminal_grouped[1,2], samp_szes[2,2] / terminal_grouped[2,2])
  }else{
    print("Not enough data points for wilcox test")
  }
  return(list(snail = snail_plot, density = density_box, length = length_box, histo_l = length_histos, histo_d = dens_histos,
              df = out_df, recombinations = clade_terminal_sans, branch_lengths = terminal_lengths, snp_data = snp_terminal,
              branch_df = by_branch_df))
}

typing_gubbins_rec <- function(reccy_preds, branch_base){
  ## Function to take in branch base and classify each snp,
  ## either r for in a putative recombination event or S for a clonal frame mutation
  
  branch_base$Type <- "S"
  cat("\n")
  num_zeros <- nchar(nrow(branch_base))
  num_rows <- nrow(branch_base)
  for(k in 1:nrow(branch_base)){
    nchar_k <- nchar(k)
    nchar_0 <- num_zeros - nchar_k
    cat("\r", "Completed ", rep(0, nchar_0), k, " of ", num_rows, " SNPs", sep = "")
    current_snp <- branch_base[k,]
    potential_rec <- reccy_preds %>%
      filter((start_node == current_snp$start_node) &  (end_node == current_snp$end_node))
    if(nrow(potential_rec) > 0){
      
      if(any(data.table::between(current_snp$base_number, potential_rec$start, potential_rec$end)))
        branch_base$Type[k] <- "r"
    }
    
  }
  
  return(branch_base)
}

typing_lapply_func <- function(row, reccy_preds, num_zeros, num_rows){
  nchar_k <- nchar(row[length(row)])
  nchar_0 <- num_zeros - nchar_k
  cat("\r", "Completed ", rep(0, nchar_0), row[length(row)], " of ", num_rows, " SNPs", sep = "")
  snp_type <- "S"
  potential_rec <- reccy_preds %>%
    filter((start_node == row[1]) &  (end_node == row[2]))
  if(nrow(potential_rec) > 0){
    
    if(any(data.table::between(as.integer(row[5]), potential_rec$start, potential_rec$end)))
      snp_type <- "r"
  }
  
  return(snp_type)
}

typing_gubbins_rec_apply <- function(reccy_preds, branch_base){
  cat("\n")
  num_zeros <- nchar(nrow(branch_base))
  num_rows <- nrow(branch_base)
  branch_base <- branch_base %>% mutate(index = row_number())
  snp_types <- apply(X = branch_base, MARGIN = 1, FUN = typing_lapply_func, 
                     reccy_preds = reccy_preds, num_zeros = num_zeros,
                     num_rows = num_rows)
  branch_base$Type <- snp_types
  cat("\n")
  return(branch_base)
}

total_func <- function(data, indices){
  d <- data[indices,]
  r_m <- sum(d$rec_snps) / sum(d$S_tot)
  rho_m <- sum(d$rec_tot) / sum(d$S_tot)
  rec_branch <-  sum(d$rec_tot) / sum(d$length)
  res_list <- c(r_m, rho_m, rec_branch)
  names(res_list) <- c("r/m","rho/m","branch_rec")
  return(res_list)
  
}
rm_func <- function(data, indices){
  d <- data[indices,]
  r_m <- sum(d$rec_snps) / sum(d$S_tot)
  return(r_m)
}
rho_m_func <- function(data, indices){
  d <- data[indices,]
  rho_m <- sum(d$rec_tot) / sum(d$S_tot)
  return(rho_m)
}
branch_length_func <- function(data, indices){
  d <- data[indices,]
  branch_rec <- sum(d$length) / sum(d$rec_tot)
  return(branch_rec)
}

bootstrapping_function <- function(branch_db, boot_num = 100000, num_threads = 1,
                                   clade_name = "cluster 2", statistic = "total", xmax = 100){
  
  if (is.null(statistic) || length(statistic) != 1L || !statistic %in% c("r/m", "rho/m","branch_length", "total")) {
    stop(sprintf("'statistic' argument must be one of: 'r/m', 'rho/m', 'branch_length' not %s \n", statistic))
  }
  
  # total_db = recombinations df from reccy_sum_up
  # colloid = test to focus on (r/m, rho/m)
  # snp_df = snp_df from reccy_sum_up
  ## Function to bootstrap for r/m or other features of a dataset 
  comM_data <- branch_db %>% filter(comM == "broken")
  wt_data <- branch_db %>% filter(comM == "wt")
  smp_sizes <- dplyr::count(branch_db, comM)
  
  if(statistic == "total"){
    comM_df_boot <- boot(data = comM_data, statistic = total_func, R = boot_num, ncpus = num_threads,
                         parallel = "multicore")
    wt_df_boot <- boot(data = wt_data, statistic = total_func, R = boot_num, ncpus = num_threads,
                       parallel = "multicore")
  }else if(statistic == "r/m"){
    comM_df_boot <- boot(data = comM_data, statistic = rm_func, R = boot_num, ncpus = num_threads)
    wt_df_boot <- boot(data = wt_data, statistic = rm_func, R = boot_num, ncpus = num_threads)
  }else if(statistic == "rho/m"){
    comM_df_boot <- boot(data = comM_data, statistic = rho_m_func, R = boot_num, ncpus = num_threads)
    wt_df_boot <- boot(data = wt_data, statistic = rho_m_func, R = boot_num, ncpus = num_threads)
  }else{
    comM_df_boot <- boot(data = comM_data, statistic = branch_length_func, R = boot_num, ncpus = num_threads)
    wt_df_boot <- boot(data = wt_data, statistic = branch_length_func, R = boot_num, ncpus = num_threads)
  }
  
  if(statistic == "total"){
    
    rm_diff <- as.data.frame(wt_df_boot$t[,1] - comM_df_boot$t[,1]) %>%
      rename(diff = 1) %>%
      mutate(cluster = clade_name)
    rhom_diff <- as.data.frame(wt_df_boot$t[,2] - comM_df_boot$t[,2])  %>%
      rename(diff = 1) %>%
      mutate(cluster = clade_name)
    branch_diff <- as.data.frame(wt_df_boot$t[,3] - comM_df_boot$t[,3])  %>%
      rename(diff = 1) %>%
      mutate(cluster = clade_name)
    
    
    
    wilc_test_rm <- wilcox.test(comM_df_boot$t[,1], wt_df_boot$t[,1])
    wilc_test_rhom <- wilcox.test(comM_df_boot$t[,2], wt_df_boot$t[,2])
    wilc_test_branch <- wilcox.test(comM_df_boot$t[,3], wt_df_boot$t[,3])
    comM_hist_data <- comM_df_boot$t %>% as.data.frame() %>%
      rename(r_m = 1, rho_m = 2, branch_length_per_rec = 3) %>% mutate(comM = "broken")
    wt_hist_data <- wt_df_boot$t %>% as.data.frame() %>%
      rename(r_m = 1, rho_m = 2, branch_length_per_rec = 3) %>% mutate(comM = "wt")
    
  }else{
    wilc_test <- wilcox.test(comM_df_boot$t[,1], wt_df_boot$t[,1])
    comM_hist_data <- comM_df_boot$t %>% as.data.frame() %>%
      rename(r_m = 1) %>% mutate(comM = "broken")
    wt_hist_data <- wt_df_boot$t %>% as.data.frame() %>%
      rename(r_m = 1) %>% mutate(comM = "wt")
  }
  
  
  tot_hist_data <- bind_rows(comM_hist_data, wt_hist_data)
  
  if(statistic == "total"){
    plot_title_rm <- sprintf("r/m for %s bootstrap samples in clade %s, p-val = %s", boot_num, clade_name, wilc_test_rm$p.value)
    
    histo_rm <- ggplot(data = tot_hist_data) + geom_histogram(data = filter(tot_hist_data, comM == "broken"),
                                                              aes(x = r_m, fill = comM), alpha = 0.5,
                                                              binwidth = 0.1) +
      geom_histogram(data = filter(tot_hist_data, comM == "wt"),aes(x = r_m, fill = comM),
                     alpha = 0.5,
                     binwidth = 0.1) +
      labs( y = "Count", x = "r/m", title = plot_title_rm) + 
      theme_bw() +
      geom_vline(xintercept = median(filter(tot_hist_data, comM == "broken")$r_m), linetype = "dotted") +
      geom_vline(xintercept = median(filter(tot_hist_data, comM == "wt")$r_m), linetype = "dashed") +
      xlim(c(0,(max(tot_hist_data$r_m) + mean(tot_hist_data$r_m))))
    
    plot_title_rhom <- sprintf("rho/m for %s bootstrap samples in clade %s, p-val = %s", boot_num, clade_name, wilc_test_rhom$p.value)
    
    rhom_bindwith <- (max(tot_hist_data$rho_m) - min(tot_hist_data$rho_m)) / (boot_num/100)
    
    histo_rhom <- ggplot(data = tot_hist_data) + geom_histogram(data = filter(tot_hist_data, comM == "broken"),
                                                                aes(x = rho_m, fill = comM), alpha = 0.5,
                                                                binwidth = rhom_bindwith) +
      geom_histogram(data = filter(tot_hist_data, comM == "wt"),aes(x = rho_m, fill = comM),
                     alpha = 0.5,
                     binwidth = rhom_bindwith) +
      labs( y = "Count", x = "rho/m", title = plot_title_rhom) + 
      theme_bw() +
      geom_vline(xintercept = median(filter(tot_hist_data, comM == "broken")$rho_m), linetype = "dotted") +
      geom_vline(xintercept = median(filter(tot_hist_data, comM == "wt")$rho_m), linetype = "dashed") +
      xlim(c(0,(max(tot_hist_data$rho_m) + mean(tot_hist_data$rho_m))))
    
    plot_title_branch <- sprintf("branch_length/rec for %s bootstrap samples in clade %s, p-val = %s", boot_num, clade_name, wilc_test_branch$p.value)
    
    branch_bindwith <- (max(tot_hist_data$branch_length_per_rec) - min(tot_hist_data$branch_length_per_rec)) / (boot_num/100)
    histo_branch <- ggplot(data = tot_hist_data) + geom_histogram(data = filter(tot_hist_data, comM == "broken"),
                                                                  aes(x = branch_length_per_rec, fill = comM), alpha = 0.5,
                                                                  binwidth = branch_bindwith) +
      geom_histogram(data = filter(tot_hist_data, comM == "wt"),aes(x = branch_length_per_rec, fill = comM),
                     alpha = 0.5,
                     binwidth = branch_bindwith) +
      labs( y = "Count", x = "branch_length/rec", title = plot_title_branch) + 
      theme_bw() +
      geom_vline(xintercept = median(filter(tot_hist_data, comM == "broken")$branch_length_per_rec), linetype = "dotted") +
      geom_vline(xintercept = median(filter(tot_hist_data, comM == "wt")$branch_length_per_rec), linetype = "dashed") +
      xlim(c(0,(max(tot_hist_data$branch_length_per_rec) + mean(tot_hist_data$branch_length_per_rec))))
    
    
    
    rm_diff_bindwith <- (max(rm_diff$diff) - min(rm_diff$diff)) / (boot_num/100)
    bounds <- quantile(rm_diff$diff, probs = c(0.025, 0.975))
    H0 <- !((bounds["2.5%"] > 0) || (bounds["97.5%"] < 0))
    plot_title_rmdiff <- sprintf("r/m diff (WT - comMdelta) for %s bootstrap samples in clade %s, H0 = %s", boot_num, clade_name, H0)
    rm_diff_plot <- ggplot(rm_diff, aes(x = diff)) +
      geom_histogram(binwidth = rm_diff_bindwith, fill = "midnightblue") +
      geom_vline(xintercept = bounds[1], colour = "red") +
      geom_vline(xintercept = bounds[2], colour = "red") +
      labs(x = "r/m difference", y = "Count", title = plot_title_rmdiff)
    
    rhom_diff_bindwith <- (max(rhom_diff$diff) - min(rhom_diff$diff)) / (boot_num/100)
    bounds <- quantile(rhom_diff$diff, probs = c(0.025, 0.975))
    H0 <- !((bounds["2.5%"] > 0) || (bounds["97.5%"] < 0))
    plot_title_rhomdiff <- sprintf("rho/m diff (WT - comMdelta) for %s bootstrap samples in clade %s, H0 = %s", boot_num, clade_name, H0)
    rhom_diff_plot <- ggplot(rhom_diff, aes(x = diff)) +
      geom_histogram(binwidth = rhom_diff_bindwith, fill = "midnightblue") +
      geom_vline(xintercept = bounds[1], colour = "red") +
      geom_vline(xintercept = bounds[2], colour = "red") +
      labs(x = "rho/m difference", y = "Count", title = plot_title_rhomdiff)
    
    branch_diff_bindwith <- (max(branch_diff$diff) - min(branch_diff$diff)) / (boot_num/100)
    bounds <- quantile(branch_diff$diff, probs = c(0.025, 0.975))
    H0 <- !((bounds["2.5%"] > 0) || (bounds["97.5%"] < 0))
    plot_title_branchdiff <- sprintf("r/m diff (WT - comMdelta) for %s bootstrap samples in clade %s, H0 = %s", boot_num, clade_name, H0)
    branch_diff_plot <- ggplot(branch_diff, aes(x = diff)) +
      geom_histogram(binwidth = branch_diff_bindwith, fill = "midnightblue") +
      geom_vline(xintercept = bounds[1], colour = "red") +
      geom_vline(xintercept = bounds[2], colour = "red") +
      labs(x = "branch length / recombination difference", y = "Count", title = plot_title_branchdiff) #+
    #xlim(c(bounds["2.5%"] - 20,(max(branch_diff$diff) + mean(branch_diff$diff))))
    
    
    histo_out <- list(r_m = histo_rm, rho_m =  histo_rhom, length_rec = histo_branch,
                      r_m_diff = rm_diff_plot, 
                      rhom_diff = rhom_diff_plot,
                      length_diff = branch_diff_plot)
    diff_data <- list(rm = rm_diff,
                      rhom = rhom_diff, 
                      branch = branch_diff)
    
  }else{
    plot_title <- sprintf("%s for %s bootstrap samples in clade %s, p-val = %s", statistic, boot_num, clade_name, wilc_test$p.value)
    
    histo_out <- ggplot(data = tot_hist_data) + geom_histogram(data = filter(tot_hist_data, comM == "broken"),
                                                               aes(x = r_m, fill = comM), alpha = 0.5,
                                                               binwidth = 0.1) +
      geom_histogram(data = filter(tot_hist_data, comM == "wt"),aes(x = r_m, fill = comM),
                     alpha = 0.5,
                     binwidth = 0.1) +
      labs( y = "Count", x = statistic, title = plot_title) + 
      theme_bw() +
      geom_vline(xintercept = median(filter(tot_hist_data, comM == "broken")$r_m), linetype = "dotted") +
      geom_vline(xintercept = median(filter(tot_hist_data, comM == "wt")$r_m), linetype = "dashed") +
      xlim(c(0,xmax))
    diff_data <- NULL
  }
  
  
  
  return(list(rorcp = comM_df_boot, wt = wt_df_boot, gg_hist = histo_out, hist_data = tot_hist_data,
              diff_data = diff_data))
  
}

stat_2_5_quant <- function(x){
  return(quantile(x, probs = 0.025))
}
stat_97_5_quant <- function(x){
  return(quantile(x, probs = 0.975))
}
UniquePanelCoords <- ggplot2::ggproto(
  "UniquePanelCoords", ggplot2::CoordCartesian,
  
  num_of_panels = 1,
  panel_counter = 1,
  panel_ranges = NULL,
  
  setup_layout = function(self, layout, params) {
    self$num_of_panels <- length(unique(layout$PANEL))
    self$panel_counter <- 1
    layout
  },
  
  setup_panel_params =  function(self, scale_x, scale_y, params = list()) {
    if (!is.null(self$panel_ranges) & length(self$panel_ranges) != self$num_of_panels)
      stop("Number of panel ranges does not equal the number supplied")
    
    train_cartesian <- function(scale, limits, name, given_range = NULL) {
      if (is.null(given_range)) {
        expansion <- ggplot2:::default_expansion(scale, expand = self$expand)
        range <- ggplot2:::expand_limits_scale(scale, expansion,
                                               coord_limits = self$limits[[name]])
      } else {
        range <- given_range
      }
      
      out <- list(
        ggplot2:::view_scale_primary(scale, limits, range),
        sec = ggplot2:::view_scale_secondary(scale, limits, range),
        arrange = scale$axis_order(),
        range = range
      )
      names(out) <- c(name, paste0(name, ".", names(out)[-1]))
      out
    }
    
    cur_panel_ranges <- self$panel_ranges[[self$panel_counter]]
    if (self$panel_counter < self$num_of_panels)
      self$panel_counter <- self$panel_counter + 1
    else
      self$panel_counter <- 1
    
    c(train_cartesian(scale_x, self$limits$x, "x", cur_panel_ranges$x),
      train_cartesian(scale_y, self$limits$y, "y", cur_panel_ranges$y))
  }
)

coord_panel_ranges <- function(panel_ranges, expand = TRUE, default = FALSE, clip = "on") {
  ggplot2::ggproto(NULL, UniquePanelCoords, panel_ranges = panel_ranges, 
                   expand = expand, default = default, clip = clip)
}

total_clade_runs_acba <- function(list_of_results = NULL,regions_to_exclude = "total", collect_phage = TRUE,
                             max_rec_length = NULL){
  
  #browser()
  ## Phage regions PHASTER:
  ### Clade 3:
  # Starts: 1735281, 3221404
  # Ends: 1782780, 3273622
  # Extended cps: 526707 .. 668530
  ### Clade 4:
  # Starts: 838178, 837745, 2187208, 2654953, 3107119
  # Ends: 847561, 874758, 2227712, 2694310, 3153545
  # Extended cps: 0 .. 250361 and 3852894 .. 3875775
  ### Clade 5:
  # Starts: 2034457, 2255791, 2294007, 2372403, 3164458, 3237503, 3350995
  # Ends: 2079596, 2273593, 2320972, 2386265, 3201623, 3254141, 3390318
  # Extended cps: 955900 .. 1036661
  ### Clade 6:
  # Starts: 161099, 1154539, 1572773
  # Ends: 209406, 1194090, 1591513
  # Extended cps: 2293602 .. 2682714
  ### GC2:
  # Starts: 1527936, 2596115, 2599929, 2698724, 2961618, 786209
  # Ends: 1582498, 2624498, 2636231, 2743967, 2996206, 801880
  ### GC1:
  # Starts: 2461407, 2682256, 811817
  # Ends: 2480441, 2733602, 837679
  # Extended cps: 3452029 .. 3688965 & 0 .. 192365
  
  if(regions_to_exclude == "total"){
    gc1_starts <- c(2461407, 2682256, 811817, 3452029, 1)
    gc1_ends <- c(2480441, 2733602, 837679, 3688965, 192365)
    gc2_starts <- c(1527936, 2596115, 2599929, 2698724, 2961618, 786209, 1, 3685790)
    gc2_ends <- c(1582498, 2624498, 2636231, 2743967, 2996206, 801880, 397167, 3908399)
    clade3_starts <- c(1735281, 3221404, 526707)
    clade3_ends <- c(1782780, 3273622, 668530)
    clade4_starts <- c(838178, 837745, 2187208, 2654953, 3107119, 0, 3852894)
    clade4_ends <- c(847561, 874758, 2227712, 2694310, 3153545, 250362, 3875775)
    clade5_starts <- c(2034457, 2255791, 2294007, 2372403, 3164458, 3237503, 3350995, 955900)
    clade5_ends <- c(2079596, 2273593, 2320972, 2386265, 3201623, 3254141, 3390318, 1036661)
    clade6_starts <- c(161099, 1154539, 1572773, 2293602)
    clade6_ends <- c(209406, 1194090, 1591513, 2682714)
    
  }else if(regions_to_exclude == "phage"){
    gc1_starts <- c(2461407, 2682256, 811817)
    gc1_ends <- c(2480441, 2733602, 837679)
    gc2_starts <- c(1527936, 2596115, 2599929, 2698724, 2961618)#, 786209)
    gc2_ends <- c(1582498, 2624498, 2636231, 2743967, 2996206)#, 801880)
    clade3_starts <- c(1735281, 3221404)
    clade3_ends <- c(1782780, 3273622)
    clade4_starts <- c(838178, 837745, 2187208, 2654953, 3107119)
    clade4_ends <- c(847561, 874758, 2227712, 2694310, 3153545)
    clade5_starts <- c(2034457, 2255791, 2294007, 2372403, 3164458, 3237503, 3350995)
    clade5_ends <- c(2079596, 2273593, 2320972, 2386265, 3201623, 3254141, 3390318)
    clade6_starts <- c(161099, 1154539, 1572773)
    clade6_ends <- c(209406, 1194090, 1591513)
  }else if(regions_to_exclude == "cps"){
    gc1_starts <- c( 3452029, 1)
    gc1_ends <- c(3688965, 192365)
    gc2_starts <- c(1, 3685790)
    gc2_ends <- c(397167, 3908399)
    clade3_starts <- c(526707)
    clade3_ends <- c(668530)
    clade4_starts <- c(0, 3852894)
    clade4_ends <- c(250362, 3875775)
    clade5_starts <- c(955900)
    clade5_ends <- c(1036661)
    clade6_starts <- c(2293602)
    clade6_ends <- c(2682714)
    
  }
  
  if (is.null(list_of_results)){
    print("On Clade GC1")
    snps_gc1 <- read.csv("~/Dropbox/phd/acinetobacter_baumannii_work/acba_gubbins/thesis_res/gc1_gubbins_3.2.0_ft_rax/gc1_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)
    tot_gc1 <- reccy_sum_up(micro_csv = "gc1_gubbins_3.2.0_ft_rax/gc1_comM_states.csv",
                            gubb_base = "gc1_gubbins_3.2.0_ft_rax/gc1_gubbins.",
                            phage_starts = gc1_starts,
                            phage_ends = gc1_ends,
                            clade_name = "GC1",
                            branch_recons = "gc1_gubbins_3.2.0_ft_rax/gc1_comM_res.tsv",
                            phage_exclude = collect_phage,
                            max_rec_length = max_rec_length,
                            classified_snps = snps_gc1)
    
    print("On Clade GC2")
    snps_gc2 <- read.csv("gc2_gubbins_3.2.0_ft_raxml_128/gc2_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)
    tot_gc2 <- reccy_sum_up(micro_csv = "gc2_gubbins_3.2.0_ft_raxml_128/gc2_comM_states.csv",
                            gubb_base = "gc2_gubbins_3.2.0_ft_raxml_128/gc2_gubbins.",
                            phage_starts = gc2_starts,
                            phage_ends = gc2_ends,
                            clade_name = "GC2",
                            branch_recons = "gc2_gubbins_3.2.0_ft_raxml_128/gc2_comM_res.tsv",
                            phage_exclude = collect_phage,
                            max_rec_length = max_rec_length,
                            classified_snps = snps_gc2)
    
    print("On Clade 3")
    snps_cl3 <- read.csv("clade_3_gubbins_res_jar/clade_3_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)
    clade_3_res <- reccy_sum_up(micro_csv = "clade_3_gubbins_res_jar/clad3_3_comM_states.csv",
                                gubb_base = "clade_3_gubbins_res_jar/clade_3_gubbins_res_jar.",
                                phage_starts = clade3_starts,
                                phage_ends = clade3_ends,
                                clade_name = "clade3",
                                branch_recons = "clade_3_gubbins_res_jar/clade_3_comM_res.tsv",
                                phage_exclude = collect_phage,
                                max_rec_length = max_rec_length,
                                classified_snps = snps_cl3)
    
    print("On Clade 4")
    snps_cl4 <- read.csv("clade_4_gubbins_res_jar/clade_4_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)
    clade_4_res <- reccy_sum_up(micro_csv = "clade_4_gubbins_res_jar/clade_4_comM_states.csv",
                                gubb_base =  "clade_4_gubbins_res_jar/clade_4_gubbins_res_jar.",
                                phage_starts = clade4_starts,
                                phage_ends = clade4_ends,
                                clade_name = "clade4",
                                branch_recons = "clade_4_gubbins_res_jar/clade_4_comM_res.tsv",
                                phage_exclude = collect_phage,
                                max_rec_length = max_rec_length,
                                classified_snps = snps_cl4)
    
    print("On Clade 5")
    snps_cl5 <- read.csv("clade_5_gubbins_res_jar/clade_5_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)
    clade_5_res <- reccy_sum_up(micro_csv = "clade_5_gubbins_res_jar/clade_5_comM_states.csv",
                                gubb_base = "clade_5_gubbins_res_jar/clade_5_gubbins_res_jar.",
                                phage_starts = clade5_starts,
                                phage_ends = clade5_ends,
                                clade_name = "clade5",
                                branch_recons = "clade_5_gubbins_res_jar/clade_5_comM_res.tsv",
                                phage_exclude = collect_phage,
                                max_rec_length = max_rec_length,
                                classified_snps = snps_cl5)
    
    print("On Clade 6")
    snps_cl6 <- read.csv("clade_6_gubbins_res_jar/clade_6_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)
    clade_6_res <- reccy_sum_up(micro_csv = "clade_6_gubbins_res_jar/clade_6_comM_states.csv",
                                gubb_base = "clade_6_gubbins_res_jar/clade_6_gubbins_res_jar.",
                                phage_starts = clade6_starts,
                                phage_ends = clade6_ends,
                                clade_name = "clade6",
                                branch_recons = "clade_6_gubbins_res_jar/clade_6_comM_res.tsv",
                                phage_exclude = collect_phage,
                                max_rec_length = max_rec_length,
                                classified_snps = snps_cl6)
  }else{
    tot_gc1 <- list_of_results[[1]]
    tot_gc2 <- list_of_results[[2]]
    clade_3_res <- list_of_results[[3]]
    clade_4_res <- list_of_results[[4]]
    clade_5_res <- list_of_results[[5]]
    clade_6_res <- list_of_results[[6]]
  }
  
  total_lengths <- ggarrange(tot_gc1$histo_l, tot_gc2$histo_l,
                             clade_3_res$histo_l, clade_4_res$histo_l,
                             clade_5_res$histo_l, clade_6_res$histo_l,
                             ncol = 3, nrow = 2, common.legend = TRUE) 
  
  total_branches <- bind_rows(tot_gc1$branch_lengths, tot_gc2$branch_lengths) %>%
    bind_rows(clade_3_res$branch_lengths) %>%
    bind_rows(clade_4_res$branch_lengths) %>%
    bind_rows(clade_5_res$branch_lengths) %>%
    bind_rows(clade_6_res$branch_lengths) 
  
  grouped_branches <- total_branches %>% group_by(comM) %>%
    summarise(tot_length = sum(length))
  
  total_reccies <-  bind_rows(tot_gc1$df, tot_gc2$df) %>%
    bind_rows(clade_3_res$df) %>%
    bind_rows(clade_4_res$df) %>%
    bind_rows(clade_5_res$df) %>%
    bind_rows(clade_6_res$df) 
  
  total_reccies_events <- c("broken" = sum(as.integer(total_reccies$sample_comM)) ,
                            "complete" = sum(as.integer(total_reccies$sample_wt)))
  
  reccies_per_length <- grouped_branches %>% mutate(total_reccies = total_reccies_events) %>%
    mutate(per_length = total_reccies / tot_length)
  
  ## Need some plots now!
  ## Get the total recombination and density boxplot via combining reccy dfs and 
  ## replotting with a gg_facet 
  total_branches <- bind_rows(tot_gc1$recombinations, tot_gc2$recombinations) %>%
    bind_rows(clade_3_res$recombinations) %>%
    bind_rows(clade_4_res$recombinations) %>%
    bind_rows(clade_5_res$recombinations) %>%
    bind_rows(clade_6_res$recombinations) %>%
    mutate(comM = ifelse(comM == "broken", "Disrupted", "WT"),
           comM = factor(comM, levels = c("WT","Disrupted")))
  total_branches$Clade <- factor(total_branches$Clade, levels = 
                                   c("GC1","GC2", "Clade3","Clade4",
                                     "Clade5","Clade6"))
  
  
  length_names <- list(
    "GC1" = bquote(bold("GC1")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[1])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[1])),
    "GC2" = bquote(bold("GC2")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[2])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[2])),
    "Clade3" = bquote(bold("Strain 3")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[3])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[3])),
    "Clade4" = bquote(bold("Strain 4")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[4])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[4])),
    "Clade5" = bquote(bold("Strain 5")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[5])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[5])),
    "Clade6" = bquote(bold("Strain 6")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[6])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[6]))
  )
  
  
  
  length_labeller <- function(variable,value){
    return(length_names[value])
  }
  
  density_names <- list(
    "GC1" = bquote(bold("GC1")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[1])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[1])),
    "GC2" = bquote(bold("GC2")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[2])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[2])),
    "Clade3" = bquote(bold("Strain 3")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[3])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[3])),
    "Clade4" = bquote(bold("Strain 4")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[4])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[4])),
    "Clade5" = bquote(bold("Strain 5")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[5])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[5])),
    "Clade6" = bquote(bold("Strain 6")*"," ~ WT ~ N ~ '=' ~ .(total_reccies$sample_wt[6])*',' ~ Disrupted ~ N ~ '=' ~ .(total_reccies$sample_comM[6]))
  )
  
  
  
  density_labeller <- function(variable,value){
    return(density_names[value])
  }
  
  total_branches <- total_branches %>% 
    mutate(clade_comm = paste(comM, Clade)) %>%
    group_by(clade_comm) %>%
    mutate(log_length = log(length, base = 10),
           log_dens = log(density, base = 10)) %>%
    mutate(outlier_log_length = (log_length > (quantile(log_length, probs = 0.75) + 
                                                 IQR(log_length)*1.5)) | (log_length < (quantile(log_length, probs = 0.25) -
                                                                                          IQR(log_length)*1.5))) %>%
    mutate(outlier_log_dens = (log_dens > (quantile(log_dens, probs = 0.75) + 
                                             IQR(log_dens)*1.5)) | 
             (log_dens < (quantile(log_dens, probs = 0.25) -
                            IQR(log_dens)*1.5))) %>%
    ungroup() %>%
    as.data.frame()
  
  
  length_box <- ggplot(data = total_branches,
                       aes(x = comM, y = length, group = comM,
                           colour = comM)) +
    geom_point(data = function(x) dplyr::filter(x, outlier_log_length),
               position = "jitter") +
    facet_wrap(~ Clade, labeller = length_labeller) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "comM", y = "Length", colour = "comM") + 
    geom_signif(comparisons = list(c("WT", "Disrupted")),
                test = "wilcox.test", step_increase = 0.075,
                map_signif_level = TRUE, tip_length = 0,
                y_position = log(400000, base = 10),
                colour = "black") +
    scale_y_log10(limits = c(1e+02, 1e+06))
  
  
  density_box <- ggplot(data = total_branches,
                        aes(x = comM, y = density, group = comM,
                            colour = comM)) +
    facet_wrap(~ Clade, labeller = length_labeller, scales = "fixed") +
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = function(x) dplyr::filter(x, outlier_log_dens),
               position = "jitter") +
    # stat_compare_means(method = "wilcox.test",
    #                    comparisons = list(c("WT", "Disrupted")),
    #                    step.increase = 0.075,
    #                    label = "p.signif" ) +
    geom_signif(comparisons = list(c("WT", "Disrupted")),
                test = "wilcox.test", step_increase = 0.075,
                map_signif_level = TRUE, tip_length = 0,
                y_position = log(0.4, base = 10),
                colour = "black") +
    labs(x = "comM", y = "SNP density", colour = "comM") +
    # scale_y_continuous(trans = scales::log10_trans(), 
    #                    breaks=10^c(-4,-3,-2,-1,0), 
    #                    labels=c("1e-04","1e-03","1e-02",
    #                             "1e-01", "1e+00"))
    scale_y_log10(limits = c(1e-04, 1e+00))
  
  
  ## Get the per branch recombination plots too
  
  branch_plot_df <- total_reccies %>%
    mutate(broken_col = "broken") %>%
    mutate(complete_col = "complete")
  
  branch_plot_df$clade <- factor(branch_plot_df$clade, levels = 
                                   c("GC1","GC2", "clade3","clade4",
                                     "clade5","clade6"))
  
  
  branch_plot <- ggplot(data = branch_plot_df) +
    geom_point(aes(x = clade, y = broken_rec_per_length, colour = broken_col), size = 4) +
    geom_point(aes(x = clade, y = complete_rec_per_length, colour = complete_col), size = 4) +
    geom_segment(aes(x = clade, xend = clade, y = broken_rec_per_length, yend = complete_rec_per_length), color = "grey") +
    coord_flip() + labs(x = "Clade", y = "Recombinations per branch length", colour = "comM") +
    scale_y_log10()
  
  rm_plot <- ggplot(data = branch_plot_df)  +
    geom_point(aes(x = clade, y = broken_rm, colour = broken_col), size = 4) +
    geom_point(aes(x = clade, y = wt_rm, colour = complete_col), size = 4) +
    geom_segment(aes(x = clade, xend = clade, y = broken_rm, yend = wt_rm), color = "grey") +
    coord_flip() + labs(x = "Clade", y = "r/m values", colour = "comM") #+
  #scale_y_log10()
  
  indiv_results <- list(gc1 = tot_gc1,
                        gc2 = tot_gc2,
                        clade3 = clade_3_res,
                        clade4 = clade_4_res, 
                        clade5 = clade_5_res, 
                        clade6 = clade_6_res)
  
  return(list(indiv_results = indiv_results, length_distos = total_lengths,
              branch_summary = reccies_per_length, branch_plot = branch_plot,
              density_box = density_box, length_box = length_box,
              summary_df = total_reccies, rm_plot = rm_plot,
              total_data = total_branches))
  
  
}

plotting_dens_length <- function(lp1_recombinations,
                                 lp2_recombinations,
                                 lp3_recombinations,
                                 combined_summary){
  #browser()
  ## Need some plots now!
  ## Get the total recombination and density boxplot via combining reccy dfs and 
  ## replotting with a gg_facet 
  total_branches <- bind_rows(lp1_recombinations, lp2_recombinations) %>%
    bind_rows(lp3_recombinations) %>%
    rename(RocRp = comM) %>%
    mutate(RocRp = ifelse(RocRp == "wt", "WT", "RocRp"),
           RocRp = factor(RocRp, levels = c("WT", "RocRp")))
  total_branches$Clade <- factor(total_branches$Clade, levels = 
                                   c("Cluster 1","Cluster 2",
                                     "Cluster 3"))
  
  
  length_names <- list(
    "Cluster 1" = bquote(bold("Strain 1")*"," ~ WT ~ N ~ '=' ~ .(combined_summary$sample_wt[1])*',' ~ RocRp ~ N ~ '=' ~ .(combined_summary$sample_comM[1])),
    "Cluster 2" = bquote(bold("Strain 2")*"," ~ WT ~ N ~ '=' ~ .(combined_summary$sample_wt[2])*',' ~ RocRp ~ N ~ '=' ~ .(combined_summary$sample_comM[2])),
    "Cluster 3" = bquote(bold("Strain 3")*"," ~ WT ~ N ~ '=' ~ .(combined_summary$sample_wt[3])*',' ~ RocRp ~ N ~ '=' ~ .(combined_summary$sample_comM[3]))
  )
  
  
  
  length_labeller <- function(variable,value){
    return(length_names[value])
  }
  
  total_branches <- total_branches %>% 
    mutate(clade_comm = paste(RocRp, Clade)) %>%
    group_by(clade_comm) %>%
    mutate(log_length = log(length, base = 10),
           log_dens = log(density, base = 10)) %>%
    mutate(outlier_log_length = (log_length > (quantile(log_length, probs = 0.75) + 
                                                 IQR(log_length)*1.5)) | (log_length < (quantile(log_length, probs = 0.25) -
                                                                                          IQR(log_length)*1.5))) %>%
    mutate(outlier_log_dens = (log_dens > (quantile(log_dens, probs = 0.75) + 
                                             IQR(log_dens)*1.5)) | 
             (log_dens < (quantile(log_dens, probs = 0.25) -
                            IQR(log_dens)*1.5))) %>%
    ungroup() %>%
    as.data.frame()
  
  
  length_box <- ggplot(data = total_branches,
                       aes(x = RocRp, y = length, group = RocRp,
                           colour = RocRp)) +
    geom_point(data = function(x) dplyr::filter(x, outlier_log_length),
               position = "jitter") +
    facet_wrap(~ Clade, labeller = length_labeller) +
    geom_boxplot(outlier.shape = NA) +
    labs(x = "RocRp", y = "Length", colour = "RocRp") + 
    geom_signif(comparisons = list(c("WT", "RocRp")),
                test = "wilcox.test", step_increase = 0.075,
                map_signif_level = TRUE, tip_length = 0,
                y_position = log(400000, base = 10),
                colour = "black") +
    scale_y_log10(limits = c(1e+02, 1e+06))
  
  
  density_box <- ggplot(data = total_branches,
                        aes(x = RocRp, y = density, group = RocRp,
                            colour = RocRp)) +
    facet_wrap(~ Clade, labeller = length_labeller, scales = "fixed") +
    geom_boxplot(outlier.shape = NA) +
    geom_point(data = function(x) dplyr::filter(x, outlier_log_dens),
               position = "jitter") +
    # stat_compare_means(method = "wilcox.test",
    #                    comparisons = list(c("WT", "Disrupted")),
    #                    step.increase = 0.075,
    #                    label = "p.signif" ) +
    geom_signif(comparisons = list(c("WT", "RocRp")),
                test = "wilcox.test", step_increase = 0.075,
                map_signif_level = TRUE, tip_length = 0,
                y_position = log(0.05, base = 10),
                colour = "black") +
    labs(x = "RocRp", y = "SNP density", colour = "RocRp") +
    # scale_y_continuous(trans = scales::log10_trans(), 
    #                    breaks=10^c(-4,-3,-2,-1,0), 
    #                    labels=c("1e-04","1e-03","1e-02",
    #                             "1e-01", "1e+00"))
    scale_y_log10(limits = c(1e-04, 1e-01))
  
  
  
  
  return(list(density = density_box,
              length = length_box))
  
}


