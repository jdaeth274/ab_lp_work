#' --- 
#' title: Bootstrapping good time
#' author: Josh D'Aeth
#' ---

#' Script to run the bootstrapping analysis on the __Legionella__ and \
#' __Acinetobacter__ collections. Functions loaded in from other script

#+ setup, include=FALSE
require(dplyr)
require(ape)
require(stringr)
require(ggplot2)
require(spatstat.utils)
require(boot)
require(ggpubr)
source("./src/bootstrapping_functions.R")

#'## Acinetobacter runs 
#' Run through the 6 acinetobacter clades. \
#' First let's run through the **GC2** clade - the biggest!
#'

classified_snps_gc2 <- read.csv("~/Dropbox/phd/acba_legion_paper/gubbins_results/gc2_gubbins_ft_raxml_128/gc2_gubbins_classified_snps.csv",
                                stringsAsFactors = FALSE)

tot_gc2 <- reccy_sum_up(micro_csv = "~/Dropbox/phd/acba_legion_paper/gubbins_results/gc2_gubbins_ft_raxml_128/gc2_comM_states.tsv",
                        gubb_base = "~/Dropbox/phd/acba_legion_paper/gubbins_results/gc2_gubbins_ft_raxml_128/gc2_gubbins.",
                        phage_starts = c(1527936, 2596115, 2599929, 2698724, 2961618),# 786209),#, 1, 3685790),
                        phage_ends = c(1582498, 2624498, 2636231, 2743967, 2996206), #801880),#, 397167, 3908399),
                        clade_name = "GC2",
                        branch_recons = "~/Dropbox/phd/acba_legion_paper/gubbins_results/gc2_gubbins_ft_raxml_128/gc2_comM_res.tsv",
                        phage_exclude = TRUE,
                        max_rec_length = NULL,
                        classified_snps =  classified_snps_gc2)
tot_gc2$df
tot_gc2$snail
gc2_tot_rm <- sum(tot_gc2$branch_df$rec_snps) / sum(tot_gc2$branch_df$S_tot)
gc2_tot_rhom <- sum(tot_gc2$branch_df$rec_tot) / sum(tot_gc2$branch_df$S_tot)
system.time(gc2_boot <- bootstrapping_function(tot_gc2$branch_df, num_threads = 8, clade_name = "GC2")) ## Uses about 13GB of mem
gc2_boot$gg_hist$r_m
gc2_boot$gg_hist$rho_m
gc2_boot$gg_hist$r_m_diff
gc2_boot$gg_hist$rhom_diff
gc2_boot$gg_hist$length_diff #+ scale_x_continuous(limits = c(0,10000))

#' ### GC1 
#' Now for the GC1 runs. 

gc1_snps <- read.csv("~/Dropbox/phd/acba_legion_paper/gubbins_results/gc1_gubbins_ft_rax/gc1_gubbins_classified_snps.csv",
                     stringsAsFactors = FALSE)

tot_gc1 <- reccy_sum_up(micro_csv = "~/Dropbox/phd/acba_legion_paper/gubbins_results/gc1_gubbins_ft_rax/gc1_comM_states.tsv",
                        gubb_base = "~/Dropbox/phd/acba_legion_paper/gubbins_results/gc1_gubbins_ft_rax/gc1_gubbins.",
                        phage_starts = c(2461407, 2682256, 811817),#, 3452029, 1),
                        phage_ends = c(2480441, 2733602, 837679),#, 3688965, 192365),
                        clade_name = "GC1",
                        branch_recons = "~/Dropbox/phd/acba_legion_paper/gubbins_results/gc1_gubbins_ft_rax/gc1_comM_res.tsv",
                        phage_exclude = TRUE,
                        max_rec_length = NULL,
                        classified_snps =  gc1_snps)
tot_gc1$df
gc1_tot_rm <- sum(tot_gc1$branch_df$rec_snps) / sum(tot_gc1$branch_df$S_tot)
gc1_tot_rhom <- sum(tot_gc1$branch_df$rec_tot) / sum(tot_gc1$branch_df$S_tot)
system.time(gc1_boot <- bootstrapping_function(tot_gc1$branch_df, num_threads = 10, clade_name = "GC1"))
gc1_boot$gg_hist$r_m
gc1_boot$gg_hist$rho_m
gc1_boot$gg_hist$r_m_diff
gc1_boot$gg_hist$rhom_diff
gc1_boot$gg_hist$length_diff #+ scale_x_continuous(limits = c(-100, 100))


#'### Clade 3
#' Now for the clade 3 results 

clade3_snps  <- read.csv("~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_3_gubbins_res/clade_3_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)

tot_clade3 <- reccy_sum_up(micro_csv = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_3_gubbins_res/clade_3_comM_states.tsv",
                           gubb_base = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_3_gubbins_res/clade_3_gubbins_res.",
                           phage_starts = c(1735281, 3221404),# 526707),
                           phage_ends = c(1782780, 3273622),#, 668530),
                           clade_name = "Clade3",
                           branch_recons = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_3_gubbins_res/clade_3_comM_res.tsv",
                           phage_exclude = TRUE,
                           max_rec_length = NULL,
                           classified_snps =  clade3_snps)
tot_clade3$df
clade3_tot_rm <- sum(tot_clade3$branch_df$rec_snps) / sum(tot_clade3$branch_df$S_tot)
clade3_tot_rhom <- sum(tot_clade3$branch_df$rec_tot) / sum(tot_clade3$branch_df$S_tot)
clade3_boot <- bootstrapping_function(tot_clade3$branch_df, num_threads = 10, clade_name = "Strain 3")
clade3_boot$gg_hist$r_m
clade3_boot$gg_hist$r_m_diff
clade3_boot$gg_hist$rhom_diff
clade3_boot$gg_hist$length_diff #+ scale_x_continuous(limits = c(-50,15))

#'### Clade 4 
#' Now for clade 4 results 

clade4_snps  <- read.csv("~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_4_gubbins_res/clade_4_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)

tot_clade4 <- reccy_sum_up(micro_csv = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_4_gubbins_res/clade_4_comM_states.tsv",
                           gubb_base = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_4_gubbins_res/clade_4_gubbins_res.",
                           phage_starts = c(838178, 837745, 2187208, 2654953, 3107119),# 0, 3852894),
                           phage_ends = c(847561, 874758, 2227712, 2694310, 3153545),# 250362, 3875775),
                           clade_name = "Clade4",
                           branch_recons = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_4_gubbins_res/clade_4_comM_res.tsv",
                           phage_exclude = TRUE,
                           max_rec_length = NULL,
                           classified_snps =  clade4_snps)
tot_clade4$df
clade4_tot_rm <- sum(tot_clade4$branch_df$rec_snps) / sum(tot_clade4$branch_df$S_tot)
clade4_tot_rhom <- sum(tot_clade4$branch_df$rec_tot) / sum(tot_clade4$branch_df$S_tot)
clade4_boot <- bootstrapping_function(tot_clade4$branch_df, num_threads = 10, clade_name = "Strain 4")
clade4_boot$gg_hist$r_m
clade4_boot$gg_hist$r_m_diff
clade4_boot$gg_hist$rhom_diff
clade4_boot$gg_hist$length_diff

#'### Clade 5 
#'Now for clade 5, 5 alive!

clade5_snps  <- read.csv("~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_5_gubbins_res/clade_5_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)

tot_clade5 <- reccy_sum_up(micro_csv = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_5_gubbins_res/clade_5_comM_states.tsv",
                           gubb_base = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_5_gubbins_res/clade_5_gubbins_res.",
                           phage_starts = c(2034457, 2255791, 2294007, 2372403, 3164458, 3237503, 3350995),# 955900),
                           phage_ends = c(2079596, 2273593, 2320972, 2386265, 3201623, 3254141, 3390318),# 1036661),
                           clade_name = "Clade5",
                           branch_recons = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_5_gubbins_res/clade_5_comM_res.tsv",
                           phage_exclude = TRUE,
                           max_rec_length = NULL,
                           classified_snps =  clade5_snps)
tot_clade5$df
clade5_tot_rm <- sum(tot_clade5$branch_df$rec_snps) / sum(tot_clade5$branch_df$S_tot)
clade5_tot_rhom <- sum(tot_clade5$branch_df$rec_tot) / sum(tot_clade5$branch_df$S_tot)
clade5_boot <- bootstrapping_function(tot_clade5$branch_df, num_threads = 10, clade_name = "Strain 5")
clade5_boot$gg_hist$r_m_diff
clade5_boot$gg_hist$rhom_diff
clade5_boot$gg_hist$length_diff# + scale_x_continuous(limits = c(-500,0))

#'### CLade 6 
#'Now for the clade 6 results 

clade6_snps  <- read.csv("~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_6_gubbins_res/clade_6_gubbins_classified_snps.csv",
                         stringsAsFactors = FALSE)

tot_clade6 <- reccy_sum_up(micro_csv = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_6_gubbins_res/clade_6_comM_states.tsv",
                           gubb_base = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_6_gubbins_res/clade_6_gubbins_res.",
                           phage_starts = c(161099, 1154539, 1572773),# 2293602),
                           phage_ends = c(209406, 1194090, 1591513), #2682714),
                           clade_name = "Clade6",
                           branch_recons = "~/Dropbox/phd/acba_legion_paper/gubbins_results/clade_6_gubbins_res/clade_6_comM_res.tsv",
                           phage_exclude = TRUE,
                           max_rec_length = NULL,
                           classified_snps =  clade6_snps)
tot_clade6$df
clade6_tot_rm <- sum(tot_clade6$branch_df$rec_snps) / sum(tot_clade6$branch_df$S_tot)
clade6_tot_rhom <- sum(tot_clade6$branch_df$rec_tot) / sum(tot_clade6$branch_df$S_tot)
clade6_boot <- bootstrapping_function(tot_clade6$branch_df, num_threads = 10, clade_name = "Strain 6")
clade6_boot$gg_hist$r_m_diff
clade6_boot$gg_hist$rhom_diff
clade6_boot$gg_hist$length_diff

#' ## Plot out the acinetobacter runs 
#' Create the combined plots for the acinetobacter runs 

rm_diff_boot <- bind_rows(gc2_boot$diff_data$rm, gc1_boot$diff_data$rm) %>%
  bind_rows(clade3_boot$diff_data$rm) %>%
  bind_rows(clade4_boot$diff_data$rm) %>%
  bind_rows(clade5_boot$diff_data$rm) %>%
  bind_rows(clade6_boot$diff_data$rm) %>%
  mutate(cluster = factor(cluster, levels = c("GC1","GC2","Strain 3","Strain 4",
                                              "Strain 5", "Strain 6")))

rhom_diff_boot <- bind_rows(gc2_boot$diff_data$rhom, gc1_boot$diff_data$rhom) %>%
  bind_rows(clade3_boot$diff_data$rhom) %>%
  bind_rows(clade4_boot$diff_data$rhom) %>%
  bind_rows(clade5_boot$diff_data$rhom) %>%
  bind_rows(clade6_boot$diff_data$rhom) %>%
  mutate(cluster = factor(cluster, levels = c("GC1","GC2","Strain 3","Strain 4",
                                              "Strain 5", "Strain 6")))

dummy_dat_1 <- data.frame(count = c(10,10), diff = range(rhom_diff_boot %>%
                                                           filter(cluster == "GC2") %>%
                                                           pull(diff)))
## Plot the r/m
ggplot(data = rm_diff_boot, aes(x = diff)) +
  geom_histogram(fill = "royalblue2", colour = "royalblue2", binwidth = 0.01) +
  #stat_summary(geom = "vline", fun.y = quantile, probs = 0.025) +
  stat_summary(aes(x = 0.1, y = diff, xintercept = after_stat(y), group = cluster), 
               fun = stat_2_5_quant, geom = "vline", colour = "red") + 
  stat_summary(aes(x = 0.1, y = diff, xintercept = after_stat(y), group = cluster), 
               fun = stat_97_5_quant, geom = "vline", colour = "red" ) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  facet_wrap(~cluster, nrow = 2, ncol = 3, scales = "free") +
  theme_bw() +
  labs(x = expression(paste(italic("r"),"/",italic("m"), " Difference")), y = "Count")

## Plot the rho/m

ggplot(data = rhom_diff_boot, aes(x = diff)) +
  geom_histogram(fill = "royalblue2",colour = "royalblue2", binwidth = 0.0005) +
  #stat_summary(geom = "vline", fun.y = quantile, probs = 0.025) +
  stat_summary(aes(x = 0.1, y = diff, xintercept = stat(y), group = cluster), 
               fun = stat_2_5_quant, geom = "vline", colour = "red") + 
  stat_summary(aes(x = 0.1, y = diff, xintercept = stat(y), group = cluster), 
               fun = stat_97_5_quant, geom = "vline", colour = "red" ) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  #coord_cartesian(xlim = c(-0.04,0.08)) +
  facet_wrap(vars(cluster), ncol = 3, nrow = 2,  scales = "free", shrink = TRUE) +
  coord_panel_ranges(panel_ranges = list(
    list(x = c(-0.15, 0.10)), # GC1 
    list(x = c(-0.15, 0.01)), # GC2
    list(NULL),
    list(x = c(-0.075, 0.075)),
    list(x = c(-0.1,0.25)),
    list(x = c(-0.05, 0.05))
    
  )) +
  theme_bw() +
  geom_blank(data = dummy_dat_1) +
  labs(x = expression(paste(rho, "/",italic(m), " Difference")), y = "Count") 
#geom_point(data = facter_bounds, x = NA, y = NA)

#'### Plot the boxplots 
#'Plot the length and SNP density boxplots

total_res_list <- list(tot_gc1,
                       tot_gc2,
                       tot_clade3,
                       tot_clade4,
                       tot_clade5,
                       tot_clade6)
boxies <- total_clade_runs_acba(list_of_results = total_res_list)

boxies$density_box +theme_bw()
boxies$length_box + theme_bw()

#' # Legionella analysis 
#' Now time to plot out the legionella results! We'll start with Clade 2
#' 

classified_snps <- read.csv("~/Dropbox/phd/acba_legion_paper/gubbins_results/lp2_gubbins/lp2_gubbins_classified_snps.csv",
                            stringsAsFactors = FALSE)

tot_lp2 <- reccy_sum_up(micro_csv = "~/Dropbox/phd/acba_legion_paper/gubbins_results/lp2_gubbins/lp2_rocrp_states.tsv",
                        gubb_base = "~/Dropbox/phd/acba_legion_paper/gubbins_results/lp2_gubbins/lp2_gubbins.",
                        phage_starts = c(1073977, 2268824),
                        phage_ends = c(1197938, 2379523),
                        clade_name = "Cluster 2",
                        branch_recons = "~/Dropbox/phd/acba_legion_paper/gubbins_results/lp2_gubbins/lp2_rocrp_res.tsv",
                        phage_exclude = TRUE,
                        max_rec_length = NULL,
                        classified_snps =  classified_snps,
                        MGE_state = "Yes")
tot_lp2$df
tot_lp2$snail
lp2_tot_rm <- sum(tot_lp2$branch_df$rec_snps) / sum(tot_lp2$branch_df$S_tot)
lp2_tot_rhom <- sum(tot_lp2$branch_df$rec_tot) / sum(tot_lp2$branch_df$S_tot)
lp2_boot <- bootstrapping_function(tot_lp2$branch_df, num_threads = 8, clade_name = "Strain 2")
lp2_boot$gg_hist$r_m_diff
lp2_boot$gg_hist$rhom_diff


#'## LP1 
#'Now for the premier clade, and it's live!

lp1_snps <- read.csv("~/Dropbox/phd/acba_legion_paper/gubbins_results/lp1_gubbins/lp1_gubbins_classified_snps.csv",
                     stringsAsFactors = FALSE)

tot_lp1 <- reccy_sum_up(micro_csv = "~/Dropbox/phd/acba_legion_paper/gubbins_results/lp1_gubbins/lp1_rocrp_states.tsv",
                        gubb_base = "~/Dropbox/phd/acba_legion_paper/gubbins_results/lp1_gubbins/lp1_gubbins.",
                        phage_starts = c(2251366,2747987),
                        phage_ends = c(2324375,2827323),
                        clade_name = "Cluster 1",
                        branch_recons = "~/Dropbox/phd/acba_legion_paper/gubbins_results/lp1_gubbins/lp1_rocrp_res.tsv",
                        phage_exclude = TRUE,
                        max_rec_length = NULL,
                        classified_snps =  lp1_snps,
                        MGE_state = "Yes")
tot_lp1$df
lp1_tot_rm <- sum(tot_lp1$branch_df$rec_snps) / sum(tot_lp1$branch_df$S_tot)
lp1_tot_rhom <- sum(tot_lp1$branch_df$rec_tot) / sum(tot_lp1$branch_df$S_tot)
system.time(lp1_boot <- bootstrapping_function(tot_lp1$branch_df, num_threads = 6, clade_name = "Strain 1"))
lp1_boot$gg_hist$r_m
lp1_boot$gg_hist$rho_m
lp1_boot$gg_hist$length_rec
lp1_boot$gg_hist$r_m_diff
lp1_boot$gg_hist$rhom_diff
lp1_boot$gg_hist$length_diff

#' ## LP3!!
#' Now for the final clade of legionella

lp3_snps <- read.csv("~/Dropbox/phd/acba_legion_paper/gubbins_results/lp3_gubbins/lp3_gubbins_classified_snps.csv",
                     stringsAsFactors = FALSE)

tot_lp3 <- reccy_sum_up(micro_csv = "~/Dropbox/phd/acba_legion_paper/gubbins_results/lp3_gubbins/lp3_rocrp_states.tsv",
                        gubb_base = "~/Dropbox/phd/acba_legion_paper/gubbins_results/lp3_gubbins/lp3_gubbins.",
                        phage_starts = c(2667588),
                        phage_ends = c(2688775),
                        clade_name = "Cluster 3",
                        branch_recons = "~/Dropbox/phd/acba_legion_paper/gubbins_results/lp3_gubbins/lp3_rocrp_res.tsv",
                        phage_exclude = TRUE,
                        max_rec_length = NULL,
                        classified_snps =  lp3_snps,
                        MGE_state = "Yes")
tot_lp3$df
lp3_tot_rm <- sum(tot_lp3$branch_df$rec_snps) / sum(tot_lp3$branch_df$S_tot)
lp3_tot_rhom <- sum(tot_lp3$branch_df$rec_tot) / sum(tot_lp3$branch_df$S_tot)

lp3_boot <- bootstrapping_function(tot_lp3$branch_df, num_threads = 6, clade_name = "Strain 3")
lp3_boot$gg_hist$rho_m
lp3_boot$gg_hist$r_m
lp3_boot$gg_hist$length_rec
lp3_boot$gg_hist$r_m_diff
lp3_boot$gg_hist$rhom_diff
lp3_boot$gg_hist$length_diff

rm_lp <- data.frame(matrix(ncol = 3, nrow = 3)) %>%
  rename(clade = 1, rm = 2, rhom = 3) %>%
  mutate(clade = c("lp1", "lp2","lp3"),
         rm = c(lp1_tot_rm, lp2_tot_rm, lp3_tot_rm),
         rhom = c(lp1_tot_rhom, lp2_tot_rhom, lp3_tot_rhom))

overall_rm_lp <- (sum(tot_lp1$branch_df$r_tot) +
                    sum(tot_lp2$branch_df$r_tot) +
                    sum(tot_lp3$branch_df$r_tot)) / 
  (sum(tot_lp1$branch_df$S_tot) +
     sum(tot_lp2$branch_df$S_tot) +
     sum(tot_lp3$branch_df$S_tot))

overall_rhom_lp <- (sum(tot_lp1$branch_df$rec_tot) +
                      sum(tot_lp2$branch_df$rec_tot) +
                      sum(tot_lp3$branch_df$rec_tot)) / 
  (sum(tot_lp1$branch_df$S_tot) +
     sum(tot_lp2$branch_df$S_tot) +
     sum(tot_lp3$branch_df$S_tot))

median(tot_lp2$recombinations$density)
median(tot_lp1$recombinations$density)
median(tot_lp1$recombinations$length)
median(tot_lp2$recombinations$length)


##############################################################
## Plot out the density and length boxplots ##################
##############################################################

combined_summary <- bind_rows(tot_lp1$df,
                              tot_lp2$df,
                              tot_lp3$df)

dens_length_plots <- plotting_dens_length(
  tot_lp1$recombinations, tot_lp2$recombinations, tot_lp3$recombinations,
  combined_summary)

dens_length_plots$density + theme_bw()
dens_length_plots$length + theme_bw()

###########################################################
## Plot out the total bootstraps  #########################
###########################################################
rm_diff_boot <- bind_rows(lp1_boot$diff_data$rm %>% mutate(rep = row_number()),
                          lp2_boot$diff_data$rm %>% mutate(rep = row_number())) %>%
  bind_rows(lp3_boot$diff_data$rm %>% mutate(rep = row_number())) %>%
  mutate(cluster = factor(cluster, levels = c("Strain 1","Strain 2",
                                              "Strain 3")))

rhom_diff_boot <- bind_rows(lp1_boot$diff_data$rhom %>% mutate(rep = row_number()),
                            lp2_boot$diff_data$rhom %>% mutate(rep = row_number())) %>%
  bind_rows(lp3_boot$diff_data$rhom %>% mutate(rep = row_number())) %>%
  mutate(cluster = factor(cluster, levels = c("Strain 1","Strain 2",
                                              "Strain 3")))



## r/m
ggplot(data = rm_diff_boot, aes(x = diff)) +
  geom_histogram(fill = "royalblue2", colour = "royalblue2", binwidth = 0.01) +
  #stat_summary(geom = "vline", fun.y = quantile, probs = 0.025) +
  stat_summary(aes(x = 0.1, y = diff, xintercept = stat(y), group = cluster), 
               fun = stat_2_5_quant, geom = "vline", colour = "red") + 
  stat_summary(aes(x = 0.1, y = diff, xintercept = stat(y), group = cluster), 
               fun = stat_97_5_quant, geom = "vline", colour = "red" ) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  facet_wrap(~cluster, nrow = 1, ncol = 3, scales = "free") +
  theme_bw() +
  labs(x = expression(paste(italic("r/m"), " Difference")), y = "Count") +
  coord_panel_ranges(panel_ranges = list(
    list(x = NULL, y = c(0,400)),
    list(x = NULL, y = c(0,125)),
    list(y = c(0,60))
  ))


dummy_dat_1 <- data.frame(count = c(10,10), diff = range(rhom_diff_boot %>% filter(cluster == "Cluster 1") %>% pull(diff)))

ggplot(data = rhom_diff_boot, aes(x = diff)) +
  geom_histogram(fill = "royalblue2",colour = "royalblue2", binwidth = 0.0005) +
  #stat_summary(geom = "vline", fun.y = quantile, probs = 0.025) +
  stat_summary(aes(x = 0.1, y = diff, xintercept = stat(y), group = cluster), 
               fun = stat_2_5_quant, geom = "vline", colour = "red") + 
  stat_summary(aes(x = 0.1, y = diff, xintercept = stat(y), group = cluster), 
               fun = stat_97_5_quant, geom = "vline", colour = "red" ) + 
  geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
  #coord_cartesian(xlim = c(-0.04,0.08)) +
  facet_wrap(vars(cluster), ncol = 3, nrow = 2,  scales = "free", shrink = TRUE) +
  coord_panel_ranges(panel_ranges = list(
    list(y = c(0, 1500)), # GC1 
    list(y = c(0, 600)), # GC2
    list(y = c(0, 175)))) +
  theme_bw() +
  geom_blank(data = dummy_dat_1) +
  labs(x = expression(paste(rho, italic("/m"), " Difference")), y = "Count") 
#geom_point(data = facter_bounds, x = NA, y = NA)

