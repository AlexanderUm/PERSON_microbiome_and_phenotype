#-------------------------------------------------------------------------------
# Parameters
#-------------------------------------------------------------------------------
prm.ls <- list()

prm.ls[["Paths"]] <- list("path_seq" = "data/seqdata/", 
                         "path_picrust2" = "data/seqdata/picrust/pathways_out/path_abun_unstrat.tsv")

prm.ls[["Paths"]][["meta"]] <- 
  list("phen" = "data/metadata/metadata_phenotype_only_with_FFQ.csv", 
       "with_contr" = "data/metadata/metadata_with_control.csv")

prm.ls[["Paths"]][["metabol"]] <- 
  list("plasma_FA" = "data/metabolites/plasma_FA.csv", 
       "plasma_metab_all" = "data/metabolites/plasma_metabolomics_all.csv", 
       "plasma_metab_sel" = "data/metabolites/plasma_metabolomics_select.csv",
       "plasma_GLP1" = "data/metabolites/plasma_GLP.csv", 
       "fecal_FA" = "data/metabolites/fecal_FA.csv")

prm.ls[["Paths"]][["additional"]][["metabol_batch_id"]] <- 
                                   "data/metadata/Metabol_Batch_Info.csv"

prm.ls[["General"]] <- list("Age_col" = "Age",
                            "Part_id_col" = "ID", 
                            "Sex_col" = "Sex",
                            "SeqID_col" = "SeqID", 
                            "Center_col" = "Center",
                            "Antib_col" = "Anitibiotics_use",
                            "Group_col" = "Phenotype", 
                            "Bristol_col" = "BSS", 
                            "Energy_col" = "Energy_Intake",
                            "Fiber_col" = "Fibers_Intake",
                            "Seed" = 34096, 
                            "ncores" = 4, 
                            "root_dir" = "out")

prm.ls[["Data"]] <- list("min_read_tax" = 10, 
                         "tax_lvls" = c("ASV", "Genus", "Family", "Phylum"), 
                         "z_tranform_col" = list("phen" =c("Age", 
                                                           "BSS", 
                                                           "Energy_Intake", 
                                                           "Fibers_Intake"), 
                                                 "with_contr" = c("Age")),
                         "metabolite_min_prev" = 0.75, 
                         "remove_NAs_taxa_glom" = TRUE)


# Alpha diversity
prm.ls[["Alpha"]] <- list("data_set_ps" = c("phen", "with_contr"), 
                          "Tax_lvl" = "ASV",
                          "Norm" = "Rare",
                          "measures" = c("Observed", "Shannon",
                                         "InvSimpson", "PhyloDiversity"),
                          "out_dir" = "out/alpha", 
                          "main_var" = "Phenotype", 
                          "fix_covar" = list("phen" = c("z_Age", 
                                                        "z_BSS", 
                                                        "z_Energy_Intake", 
                                                        "z_Fibers_Intake"), 
                                             "with_contr" = c("z_Age")),
                          "rand_vars" = list("phen" = c("Center", "Sex"), 
                                             "with_contr" = c("Center", "Sex")))

prm.ls[["Alpha"]][["measures_trans"]][[1]] <- list("measures" = c("InvSimpson"), 
                                                   "trans" = function(x){sqrt(x)})

prm.ls[["Alpha"]][["measures_trans"]][[2]] <- list("measures" = c("Shannon"), 
                                                   "trans" = function(x){x^2})

# Beta diversity 
prm.ls[["Beta"]] <- list("data_set_ps" = c("phen", "with_contr"), 
                         "Tax_lvl" = "ASV",
                         "Norm" = "CSS",
                         "distances" = c("Bray-Curtis" = "bray", 
                                         "Jaccard" = "jaccard",
                                         "Weighted UniFrac" = "wunifrac",
                                         "Unweighted UniFrac" = "unifrac"),
                           "n_perm" = 999, 
                           "main_var" = "Phenotype", 
                           "fix_covar" = list("phen" = c("z_Age", 
                                                         "z_BSS", 
                                                         "z_Energy_Intake", 
                                                         "z_Fibers_Intake", 
                                                         "Center", 
                                                         "Sex"), 
                                              "with_contr" = c("z_Age", 
                                                               "Center", 
                                                               "Sex")),
                            "out_path" = "out/beta")

prm.ls[["DA"]]  <- list("data_set_ps" = c("phen"), 
                        "tax_lvl" = c("ASV", "Genus"),
                        "ps_norm" = "Raw", 
                        "plot_norm" = "TMM", 
                        "summary_norm" = "Relat",
                        "feature_min_prev" = 0.25, 
                        "p_adj_method" = "BH",
                        "out_dir" = "out/DA", 
                        "max_qval" = 0.1, 
                        "maas_rand" = c("Center", "Sex"), 
                        "maas_main_effect" = "Phenotype",
                        "maas_covar" =  c("z_Age", 
                                          "z_BSS", 
                                          "z_Energy_Intake", 
                                          "z_Fibers_Intake"), 
                        "maas_method_prm" = list("set1" = c("norm" = "TMM", 
                                                            "trans" = "NONE", 
                                                            "method" = "ZINB"))
                        )


prm.ls[["RF"]] <- list("data_set_ps" = c("phen"),
                       "out_dir" = "out/RF", 
                       "tax_lvls" =  c("ASV", "Genus"), 
                       "feature_min_prev" = c(0.25), 
                       "n_permut" = 999, 
                       "max_pval_impotance" = 0.05, 
                       "n_trees" = 2501,
                       "gr_size_prop" = 0.35, 
                       "auto_tune_mtry" = FALSE, 
                       "all_feat_sets" = c("Taxa", "Metabol"),
                       "sig_feat_sets" = c("Taxa", "Metabol", "Combs"),
                       "metabol_tabs_indRF" = c("fecal_FA", 
                                                "plasma_metab_all"), 
                       "metabol_count_trans" = function(x){log10(x+1)},
                       "ASV_norm" = c("CSS"), 
                       "sets_combRF" = list("ASV & fecal FA" = 
                                              list("Taxa" = c("ASV--CSS--0.25"), 
                                                   "Metabol" = c("fecal_FA")), 
                                            "ASV & plasma" = 
                                              list("Taxa" = c("ASV--CSS--0.25"), 
                                                   "Metabol" = c("plasma_metab_all")), 
                                            "ASV & fecal FA & plasma" = 
                                              list("Taxa" = c("ASV--CSS--0.25"), 
                                                   "Metabol" = c("plasma_metab_all", 
                                                                 "fecal_FA"))), 
                       "plot_roc_size" = c("h" = 3.5, "w" = 5), 
                       "plot_signif_size" = c("h_coef" = 0.2, 
                                              "h_add" = 2, 
                                              "w" = 12), 
                       "plot_top_n_signif" = 10)



prm.ls[["metab"]] <- list("dir_out" = "out/metabolites",
                          "min_prev_metabol" = 0.75,
                          "data_set_ps" = c("phen"),
                          "tax_lvls" =  "Genus",
                          "ASV_norm" = "CSS",
                          "model_type" = c("LMM"), # options: c("LM", "LMM")
                          "da_tabs" = c("plasma_FA", 
                                        "plasma_metab_all", 
                                        "plasma_metab_sel", 
                                        "plasma_GLP1", 
                                        "fecal_FA"),
                          "da_trans_count" = function(x){log10(x+1)}, 
                          "da_rand" = c("Center", "Sex"), 
                          "da_rand_specific" = c("plasma_metab_all" = "Batch", 
                                                 "plasma_metab_sel" = "Batch"),
                          "da_main_effect" = "Phenotype",
                          "da_covar" = c("z_Energy_Intake", 
                                         "z_Age", 
                                         "z_BSS", 
                                         "z_Fibers_Intake"), 
                          "da_p_adj_method" = "BH",
                          "da_max_qval" = 0.1,
                          "da_ncol_plots" = 5, 
                          "cca_tabs" = c("plasma_metab_all", 
                                         "plasma_metab_sel", 
                                         "fecal_FA", 
                                         "plasma_FA"),
                          "cca_trans_count" = function(x){log10(x+1)},
                          "cca_main" = "Phenotype", 
                          "cca_covar" = c("Energy_Intake", 
                                          "Age", 
                                          "BSS", 
                                          "Fibers_Intake", 
                                          "Center", 
                                          "Sex"),
                          "cca_covar_specific" = c("plasma_metab_all" = "Batch", 
                                                   "plasma_metab_sel" = "Batch"), 
                          "cca_nperm" = 999,
                          "cca_max_pval" = 0.05, 
                          "cca_plot_dim" = list("w" = 7.5, "h"=7.5),
                          "cor_method" = "spearman",
                          "cor_balance_size" = FALSE,
                          "cor_tax_lvl" = c("Genus"), 
                          "cor_tax_norm" = c("CSS"), 
                          "cor_strata" = "Phenotype",
                          "cor_tax_prev" = 0.25,
                          "cor_padj_filt_est" = 0.3,
                          "cor_metab_tab" = c("plasma_metab_all", 
                                              "fecal_FA", 
                                              "plasma_FA"), 
                          "cor_p_adj_method" = "BH",
                          "cor_filt_plot_prm" = c("abs(estimate)>=0.1", 
                                                  "qval_overall<=0.1"), 
                          "cor_network_min_est" = c(0.1, 0.2, 0.25, 
                                                    0.3, 0.35, 0.4),
                          "cor_network_max_sig" = c("p.value<=0.005", 
                                                    "p.value<=0.001"),
                          "cor_network_jitter" = 0.75,
                          "cor_netheat_val" = c("estimate", "qval_overall"),
                          "cor_full_heat_min_est" = 0,
                          "cor_heat_val" = "qval_overall",
                          "cor_heat_val_cut" = 0.1,
                          "feature_min_prev" = 0.25, 
                          "max_corr_qval" = 0.1)

  
#-------------------------------------------------------------------------------
# Load libraries 
#-------------------------------------------------------------------------------
libs.list <- c("phyloseq", "tidyverse", "metagenomeSeq", "qiime2R",
               "ggsignif", "broom", "MicrobiomeStat", "vegan", "cowplot", 
               "ComplexHeatmap", "FSA", "rfPermute", "Maaslin2", 
               "network", "ggnetwork")

for (i in libs.list) {library(i, character.only = TRUE, )}

# Create directory 
dir.create("out/supp", recursive = TRUE, showWarnings = FALSE)

# Write objects 
save(list = c("prm.ls"), 
     file = "out/supp/prm.Rdata")

# Clean environment 
rm(list = ls())
gc()
