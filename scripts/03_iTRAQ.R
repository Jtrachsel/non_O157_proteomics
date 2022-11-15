library(tidyverse)
library(vegan)
library(broom)
library(NOMAD)

# function for running tests on clr normalized data
run_tests <- 
  function(matrix){
    # browser()
    tst <- 
      matrix %>%
      as.data.frame() %>% 
      rownames_to_column(var='ID') %>% 
      separate(ID, into=c('strain', 'condition'),remove = FALSE) %>% 
      mutate(condition=factor(condition, levels = c('vitro', 'vivo'))) %>% 
      pivot_longer(names_to = 'accno', values_to = 'intensity', -c(ID, strain, condition)) %>% 
      group_by(accno) %>% 
      nest() %>% 
      mutate(TTEST=map(.x=data,
                       .f=~t.test(formula=intensity ~ condition, data=.x) %>% tidy())) %>% 
      dplyr::select(-data) %>% 
      unnest(TTEST) %>% 
      ungroup() %>% 
      mutate(FDR=p.adjust(p.value, method = 'fdr')) %>% 
      arrange(estimate)
    return(tst)
  }

# function for running NOMAD normalization and ttests
NOMAD_diff_expression_tests <- 
  function(long_peptides){
    
    NOMAD_NORM <- nomadNormalization(y=long_peptides$Abundance, x=long_peptides, factors = list('Peptide', 'iTRAQ'))
    nomad_proteins <- nomadAssembleProteins(NOMAD_NORM$y, NOMAD_NORM$x)
    
    TESTS <- 
      nomad_proteins$scores %>% 
      as.data.frame() %>% 
      rownames_to_column(var='accno') %>% 
      pivot_longer(cols = -accno, values_to = 'intensity') %>% 
      mutate(reporter=sub('Run1_iTRAQ','',name),
             condition=ifelse(reporter %in% c(1,2,3), 'vitro', 'vivo'), 
             condition=factor(condition, levels = c('vitro', 'vivo')),
             strain=case_when(
               reporter %in% c(1,4) ~'O26', 
               reporter %in% c(2,5) ~'O111', 
               reporter %in% c(3,6) ~'O145'
             )) %>% 
      group_by(accno) %>%
      mutate(VAR=var(intensity)) %>% 
      filter(VAR != 0) %>%#pull(accno) %>% unique()
      nest() %>% 
      mutate(TTEST=map(.x=data,
                       .f=~t.test(formula=intensity ~ condition, data=.x) %>% tidy())) %>% 
      dplyr::select(-data) %>% 
      unnest(TTEST) %>% 
      ungroup() %>% 
      mutate(FDR=p.adjust(p.value, method = 'fdr')) %>% 
      arrange(estimate)
    return(TESTS)
  }


run_all_norm_tests <- function(MATRIX, NOMAD_TESTS){
  
  MATRIX <- MATRIX[,colSums(MATRIX) != 0]
  MATRIX_log <- apply(MATRIX, c(1,2), FUN = log)
  MATRIX_log_ratio <- t(scale(t(MATRIX_log), scale=FALSE))
  MATRIX_clr <- decostand(x = MATRIX, method = 'clr')
  MATRIX_rab <- decostand(x = MATRIX, method = 'total')
  MATRIX_log_rab <- decostand(x = MATRIX_log, method = 'total')
  
  # hist(MATRIX_log_rab)
  
  raw_tests <- 
    run_tests(MATRIX) %>% 
    mutate(order=1:n(), 
           method='raw')
  
  
  nomad_tests <- 
    NOMAD_TESTS %>% 
    mutate(order=1:n(), 
           method='nomad')
  
  clr_tests <-
    run_tests(MATRIX_clr) %>%
    mutate(order=1:n(), 
           method='clr')
  
  rab_tests <-  run_tests(MATRIX_rab) %>%
    mutate(order=1:n(), 
           method='rab')
  
  log_rab_tests <-  run_tests(MATRIX_log_rab) %>%
    mutate(order=1:n(), 
           method='log_rab')
  
  log_tests <- 
    run_tests(MATRIX_log) %>% 
    mutate(order=1:n(), 
           method='log')
  
  
  
  all_tests <- 
    bind_rows(raw_tests, nomad_tests,clr_tests, rab_tests, log_rab_tests, log_tests) %>% 
    group_by(accno) %>% 
    mutate(MEAN_ORDER=mean(order), 
           ref_order=order[method=='clr']) %>% 
    mutate(SIG=ifelse(p.value < 0.05, T, F))
  
  return(all_tests)
  
  
}


run_multivariate <- 
  function(MATRIX){
    META <- 
      tibble(ID=rownames(MATRIX)) %>% 
      separate(ID, into=c('strain', 'condition'), remove = F)
    
    
    plot(metaMDS(MATRIX, distance = 'robust.aitchison', autotransform = FALSE))
    NMDS <- metaMDS(MATRIX, distance = 'robust.aitchison', autotransform = FALSE)
    
    prot_NMDS_coords <- NMDS$species %>% as.data.frame() %>% rownames_to_column(var='accno')
    
    
    META <- 
      NMDS$points %>%
      as.data.frame() %>%
      rownames_to_column(var='ID') %>% 
      left_join(META)
    
    NMDS_plot <- 
      META %>%  
      ggplot(aes(x=MDS1, y=MDS2, fill=condition, label=strain)) +
      geom_point(size=3, shape=21)+
      geom_label_repel(force = 10,box.padding = unit(1,'lines'), min.segment.length = unit(0, 'lines')) +
      theme_bw() + 
      lims(x=c(-8,8), y=c(-6,6)) + 
      ggtitle('NMDS of Aitchison distances')
    NMDS_plot
    
    adonis_results <- adonis2(method = 'robust.aitchison',data = META, formula = iTRAQ_L_mat ~ strain * condition)
    return(list(NMDS_plot, adonis_results))
  }



### strainwise comparisons ###


strain_comparisons <- function(MATRIX){
  # browser()
  
  MATRIX <- MATRIX[,colSums(MATRIX) != 0]
  MATRIX_clr <- decostand(x = MATRIX, method = 'clr')
  
  RESULTS <- 
    MATRIX_clr %>%
    as.data.frame() %>% 
    rownames_to_column(var='ID') %>% 
    pivot_longer(cols=-ID, names_to = 'accno', values_to = 'clr') %>% 
    mutate(strain=sub('(.*)_(.*)','\\1',ID), 
           condition=sub('(.*)_(.*)','\\2',ID)) %>% 
    group_by(strain, accno) %>% 
    summarise(condition_dif=clr[2] - clr[1],  
              enriched_in=ifelse(condition_dif > 0 , condition[2], condition[1]),
              .groups='drop') #%>% 
    # pull(strain_ratio) %>% hist()
  return(RESULTS)
    
  
}


# annotations from uniprot

uniprot_annotations <- read_tsv('output/pan_0.05_annotations.tsv') %>% 
  mutate(accno=From) %>% 
  dplyr::select(accno, everything(), -From, -`pH dependence`)

prot_descipts <- read_tsv('output/protein_descriptions.tsv')

# output of maxquant, protein groups
# filtered to only consider proteins that make the cutoffs

PG <- 
  read_tsv('output/protein_groups_cleaned.tsv') %>% 
  filter(`Q-value` < 0.05) %>%
  filter(Peptides > 1) %>%
  filter(!grepl('CON__', accno))

valid_accnos <- PG %>% pull(accno)
##########

# Peptide level intensities for NOMAD 
# 
# ALL_PEPTIDES <- read_tsv('maxquant_results/combined/txt/peptides.txt')
ALL_PEPTIDES <- read_tsv('maxquant_results_pan_pro_0.05_pep_0.01/combined/txt/peptides.txt')

ALL_PEPTIDES <-
  ALL_PEPTIDES %>% 
  dplyr::select(Sequence,`Leading razor protein`,  matches('Reporter intensity . [a-z]+')) %>% 
  mutate(Protein=`Leading razor protein`, 
         Peptide=Sequence) %>% 
  dplyr::select(Protein, Peptide,everything(), -`Leading razor protein`, -Sequence) %>% 
  mutate(
    tmp_id=str_split(gsub('REV__tr','',Protein), pattern = '\\|'), 
    first_underscore=map_int(tmp_id, ~min(grep('_', .x)))) %>%  
  # dplyr::select(first_underscore, tmp_id) %>% 
  mutate(TMP=map2_chr(tmp_id,first_underscore, ~pluck(.x,.y) ), 
         TMP=sub(';$','',TMP), 
         TMP=sub(';[trsp]+$','',TMP),
         TMP=sub('REV__', '', TMP), 
         Protein=TMP) %>% 
  dplyr::select(-first_underscore, -tmp_id, -TMP) %>%
  filter(Protein %in% valid_accnos)



LONG_PEPTIDES <- 
  ALL_PEPTIDES %>%
  pivot_longer(cols = -c(Protein, Peptide), names_to = 'Sample', values_to = 'Abundance') %>% 
  mutate(iTRAQ=sub('Reporter intensity ([0-8]+) [a-z]+','\\1',Sample), 
         Run=sub('Reporter intensity ([0-8]+) ([a-z]+)','\\2',Sample),) %>% 
  dplyr::select(-Sample) %>% 
  filter(iTRAQ %in% c(1:6))%>%
  filter(!grepl('CON', Protein))
# summary_stats


protein_summary <- 
  LONG_PEPTIDES %>%
  filter(Abundance > 0) %>% 
  group_by(Run,Protein) %>%
  summarise(num_peptides=sum(length(unique(Peptide))), .groups = 'drop')

run_summary <- 
  protein_summary %>% 
  filter(num_peptides >1) %>% 
  group_by(Run) %>%
  summarise(tot_proteins=sum(length(unique(Protein))))



unique_proteins <- 
  protein_summary %>% 
  filter(num_peptides >1) %>% 
  summarise(lactation=sum(!(Protein[Run == 'lactation'] %in% Protein[Run == 'maintenance'])), 
            maintenance=sum(!(Protein[Run == 'maintenance'] %in% Protein[Run == 'lactation'])), 
            both=sum((Protein[Run == 'lactation'] %in% Protein[Run == 'maintenance']))) %>% 
  pivot_longer(cols=everything(),names_to = 'Run', values_to = 'unique_proteins')

### USE THIS TABLE
# unique_proteins %>% add_row(Run='total', unique_proteins=sum(unique_proteins$unique_proteins))


lact_peptides <- 
  LONG_PEPTIDES %>%
  filter(Run == 'lactation' & Abundance > 0) %>%
  # filter(Run == 'lactation') %>% 
  mutate(Run=1)


maint_peptides <-
  LONG_PEPTIDES %>%
  filter(Run == 'maintenance'& Abundance > 0) %>%
  # filter(Run == 'lactation') %>% 
  mutate(Run=1)



NOMAD_LACT_TESTS <- NOMAD_diff_expression_tests(long_peptides = lact_peptides)
NOMAD_MAINT_TESTS <- NOMAD_diff_expression_tests(long_peptides = maint_peptides)


p_NOMAD_LACT_HIST <- 
  NOMAD_LACT_TESTS %>%
  # filter(p.value < 0.05) %>%
  ggplot(aes(x=estimate, fill=ifelse(p.value < 0.05, T, F))) +
  annotate(x=2, y=50, geom='label', label='Enriched in vivo')+
  annotate(x=-2, y=50, geom='label', label='Enriched in vitro')+
  geom_histogram() + xlim(-3,3)+
  geom_vline(xintercept = 0) + 
  labs(y='count', 
       fill='P < 0.05')+
  theme(legend.position = 'bottom')

p_NOMAD_MAINT_HIST <- 
  NOMAD_MAINT_TESTS %>% 
  # filter(p.value < 0.05) %>%
  ggplot(aes(x=estimate, fill=ifelse(p.value < 0.05, T, F))) +
  annotate(x=2, y=50, geom='label', label='Enriched in vivo')+
  annotate(x=-2, y=50, geom='label', label='Enriched in vitro')+
  geom_histogram() + xlim(-3,3)+
  geom_vline(xintercept = 0) + 
  labs(y='count', 
       fill='P < 0.05')+
  theme(legend.position = 'bottom')
# All 'sigs' are more abundant in the vivo condition in both the lact and maint diets



# This section for t.tests with other normalization techniques
# 
# tibble(strain=paste0('Strain', rep(c(1:3),4)), 
#        condition=rep(c(rep('vivo', 3), rep('vitro', 3)),2), 
#        diet=c(rep('Lact', 6), rep('maint', 6)), 
#        LC_MSMS_run=c(rep('Run1',6), rep('Run2',6)), 
#        iTRAQ_label=c(1:6, 1:6))
#

iTRAQ <- 
  PG %>% 
  dplyr::select(accno, contains('Reporter intensity corrected')) %>% 
  pivot_longer(-accno, names_to = c('reporter', 'treatment'),
               names_prefix = 'Reporter intensity corrected ', 
               names_sep = ' ', values_to = 'intensity')

#  itraq reporter intensities
iTRAQ %>% 
  ggplot(aes(x=reporter, y=intensity)) + 
  geom_col(color='black') + 
  facet_wrap(~treatment)



iTRAQ <- 
  iTRAQ %>%
  filter(reporter %in% c(1:6)) %>% # only 6 of the 8 reporters used
  group_by(reporter, treatment) %>% 
  mutate(rel_intensity=intensity/sum(intensity),
         log_intensity=log(intensity),
         # log_ratio_intensity=log_intensity/sum(log_intensity),
         condition=ifelse(reporter %in% c(1,2,3), 'vitro', 'vivo'), 
         strain=case_when(
           reporter %in% c(1,4) ~'O26', 
           reporter %in% c(2,5) ~'O111', 
           reporter %in% c(3,6) ~'O145'
         )) %>% 
  ungroup() 
  
# this shows us that within each MSMS run we only dectect proteins that
# were detectable in all strains across both in-vivo and in-vitro conditions.

# iTRAQ %>%
#   group_by(treatment, accno) %>% 
#   summarise(ANY=any(intensity == 0), 
#             ALL=all(intensity == 0)) %>% 
#   ungroup() %>% 
#   mutate(ANY_NOT_ALL=ANY & !ALL) %>% 
#   filter(ANY_NOT_ALL)

###

iTRAQ_L <-
  iTRAQ %>% 
  filter(treatment == 'lactation' & is.finite(log_intensity)) %>% 
  group_by(reporter) %>% 
  mutate(norm_log_intensity=log_intensity - mean(log_intensity)) %>% 
  ungroup()



iTRAQ_L %>%
  ggplot(aes(x=reporter, y=norm_log_intensity, fill=condition)) + 
  geom_col(color='black') + 
  facet_wrap(~treatment)



iTRAQ_L_mat <- 
  iTRAQ_L %>% 
  dplyr::select(accno, strain, condition, intensity) %>% 
  pivot_wider(names_from = c(strain, condition), values_from = intensity ) %>% 
  column_to_rownames('accno') %>% as.matrix() %>% t()
  
iTRAQ_L_mat <- iTRAQ_L_mat[,colSums(iTRAQ_L_mat) != 0]

iTRAQ_M <- iTRAQ %>% filter(treatment == 'maintenance')


iTRAQ_M_mat <- 
  iTRAQ_M %>% 
  dplyr::select(accno, strain, condition, intensity) %>% 
  pivot_wider(names_from = c(strain, condition), values_from = intensity ) %>% 
  column_to_rownames('accno') %>% as.matrix() %>% t()

iTRAQ_M_mat <- iTRAQ_M_mat[,colSums(iTRAQ_M_mat) != 0]



# Run all normalizations and tests:
# need to supply nomad tests separately


lact_all_tests <- run_all_norm_tests(MATRIX =iTRAQ_L_mat, NOMAD_TESTS = NOMAD_LACT_TESTS )
lact_all_tests <- lact_all_tests %>% left_join(prot_descipts)

lact_all_tests %>% filter(method == 'clr') %>% filter(FDR < 0.05)
lact_all_tests %>% filter(method == 'nomad') %>% filter(FDR < 0.05)

p_lact_tests <- 
  lact_all_tests %>%
  filter(method %in% c('clr', 'nomad')) %>% 
  ggplot(aes(x=estimate, fill=SIG))+
  geom_histogram()+
  facet_wrap(~method, scales = 'free') + 
  geom_vline(xintercept = 0) +
  labs(fill='uncorrected P < 0.05') + 
  theme(legend.position = 'bottom') + 
  ggtitle('Histograms of log2(fold change) for all proteins', 
          'Negative estimates indicate greater expression in-vitro')


lact_all_tests <- 
  lact_all_tests%>%
  ungroup() %>% 
  transmute(gene_name = accno, 
            method,
            estimate, 
            p.value, 
            FDR,description, Localization ) %>% 
  mutate(across(where(is.numeric), ~signif(.x, digits = 2))) %>% 
  mutate(description=sub('([^=]) [A-Z][A-Z]=.*','\\1',description)) 




CLR_LACT_VIVO <-
  lact_all_tests %>% 
  filter(method == 'clr') %>% 
  filter(p.value < 0.05) %>%
  arrange(FDR) %>%
  filter(estimate > 0)

CLR_LACT_VIVO %>%
  left_join(uniprot_annotations, by = c('gene_name' = 'accno')) %>%
  write_tsv('output/CLR_LACT_VIVO.tsv')


CLR_LACT_VITRO <-
  lact_all_tests %>% 
  filter(method == 'clr') %>% 
  filter(p.value < 0.05) %>%
  arrange(FDR) %>%
  filter(estimate < 0) 

CLR_LACT_VITRO %>%
  left_join(uniprot_annotations, by = c('gene_name' = 'accno')) %>%
  write_tsv('output/CLR_LACT_VITRO.tsv')



NOMAD_LACT_VIVO <-
  lact_all_tests %>% 
  filter(method == 'nomad') %>% 
  filter(p.value < 0.05) %>%
  arrange(FDR) %>%
  filter(estimate > 0) 

NOMAD_LACT_VIVO %>%
  left_join(uniprot_annotations, by = c('gene_name' = 'accno')) %>%
  write_tsv('output/NOMAD_LACT_VIVO.tsv')


NOMAD_LACT_VITRO <-
  lact_all_tests %>% 
  filter(method == 'nomad') %>% 
  filter(p.value < 0.05) %>%
  arrange(FDR) %>%
  filter(estimate < 0)

NOMAD_LACT_VITRO %>%
  left_join(uniprot_annotations, by = c('gene_name' = 'accno')) %>%
  write_tsv('output/NOMAD_LACT_VITRO.tsv')




maint_all_tests <- run_all_norm_tests(MATRIX =iTRAQ_M_mat, NOMAD_TESTS = NOMAD_MAINT_TESTS )
maint_all_tests <- maint_all_tests %>% left_join(prot_descipts)

# only one sig in maint (nomad)
maint_all_tests %>% filter(method == 'clr') %>% filter(FDR < 0.05)
maint_all_tests %>% filter(method == 'nomad') %>% filter(FDR < 0.05)


p_maint_tests <- 
  maint_all_tests %>%
  filter(method %in% c('clr', 'nomad')) %>% 
  ggplot(aes(x=estimate, fill=SIG))+
  geom_histogram()+
  facet_wrap(~method, scales = 'free') + 
  geom_vline(xintercept = 0) +
  labs(fill='uncorrected P < 0.05') + 
  theme(legend.position = 'bottom') + 
  ggtitle('Histograms of log2(fold change) for all proteins', 
          'Negative estimates indicate greater expression in-vitro')




maint_all_tests <- 
  maint_all_tests %>% ungroup() %>% 
  transmute(gene_name = accno, 
            method,
            estimate, 
            p.value, 
            FDR,description, Localization ) %>% 
  mutate(across(where(is.numeric), ~signif(.x, digits = 2))) %>% 
  mutate(description=sub('([^=]) [A-Z][A-Z]=.*','\\1',description)) 

CLR_MAINT_VIVO <-
  maint_all_tests %>% 
  filter(method == 'clr') %>% 
  filter(p.value < 0.05) %>%
  arrange(FDR) %>%
  filter(estimate > 0) 

CLR_MAINT_VIVO %>%
  left_join(uniprot_annotations, by = c('gene_name' = 'accno')) %>% 
  write_tsv('output/CLR_MAINT_VIVO.tsv')


CLR_MAINT_VITRO <-
  maint_all_tests %>% 
  filter(method == 'clr') %>% 
  filter(p.value < 0.05) %>%
  arrange(FDR) %>%
  filter(estimate < 0) 

CLR_MAINT_VITRO %>%
  left_join(uniprot_annotations, by = c('gene_name' = 'accno')) %>%
  write_tsv('output/CLR_MAINT_VITRO.tsv')



NOMAD_MAINT_VIVO <-
  maint_all_tests %>% 
  filter(method == 'nomad') %>% 
  filter(p.value < 0.05) %>%
  arrange(FDR) %>%
  filter(estimate > 0) 

NOMAD_MAINT_VIVO %>%
  left_join(uniprot_annotations, by = c('gene_name' = 'accno')) %>% 
  write_tsv('output/NOMAD_MAINT_VIVO.tsv')


NOMAD_MAINT_VITRO <-
  maint_all_tests %>% 
  filter(method == 'nomad') %>% 
  filter(p.value < 0.05) %>%
  arrange(FDR) %>%
  filter(estimate < 0) 

NOMAD_MAINT_VIVO %>%
  left_join(uniprot_annotations, by = c('gene_name' = 'accno')) %>%
  write_tsv('output/NOMAD_MAINT_VITRO.tsv')





# pull out nomad and clr data for each 




### end differential abundance  

# multivariate similarity  

lact_multivariate <- run_multivariate(MATRIX =iTRAQ_L_mat)

maint_multivariate <- run_multivariate(MATRIX = iTRAQ_M_mat)


### strain ratio analysis

maint_within_strain_comps <- 
  strain_comparisons(MATRIX = iTRAQ_M_mat)



maint_prot_SDs <- 
  maint_within_strain_comps %>%
  dplyr::select(-enriched_in) %>% group_by(accno) %>% 
  summarise(SD=sd(condition_dif)) %>% 
  arrange(desc(SD))

### NEED TO RE-ORG THIS AND JOIN IN ANNOTATIONS FOR OUTPUT TABLE
library(ComplexHeatmap)

p_maint_within_strain <- 
  maint_within_strain_comps %>%
  pivot_wider(-enriched_in, names_from = strain, values_from = condition_dif) %>% 
  left_join(maint_prot_SDs) %>% 
  arrange(desc(SD)) %>% 
  left_join(uniprot_annotations) %>% 
  write_tsv('output/maint_within_strains.tsv') %>% 
  slice_max(SD, n=50) %>% 
  dplyr::select(accno, O111, O145,O26 ) %>% 
  column_to_rownames(var='accno') %>% 
  as.matrix() %>%# hist()
  ComplexHeatmap::Heatmap(name = 'Log2(Vivo / Vitro)',
                          column_title = 'Maintenance: within strain ratios between conditions.
                          top 50 proteins with most variation between strains, 
                          (proteins where strains behave differently between conditions)', 
                          row_names_gp = gpar(fontsize = 8))


# pivot_wider(names_from = 'strain', values_from = 'condition_dif')


lact_within_strain_comps <- 
  strain_comparisons(MATRIX = iTRAQ_L_mat) 

p_lact_within_strain <- 
  lact_within_strain_comps %>% 
  pivot_wider(-enriched_in, names_from = strain, values_from = condition_dif) %>% 
  left_join(maint_prot_SDs) %>% 
  arrange(desc(SD)) %>% 
  left_join(uniprot_annotations) %>% 
  write_tsv('output/lact_within_strains.tsv') %>% 
  slice_max(SD, n=50) %>% 
  dplyr::select(accno, O111, O145,O26 ) %>% 
  column_to_rownames(var='accno') %>% 
  as.matrix() %>%# hist()
  ComplexHeatmap::Heatmap(name = 'Log2(Vivo / Vitro)',
                          column_title = 'Lactation: within strain ratios between conditions.
                          top 50 proteins with most variation between strains, 
                          (proteins where strains behave differently between conditions)', 
                          row_names_gp = gpar(fontsize = 8))




