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
      select(-data) %>% 
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
             condition=ifelse(reporter %in% c(1,2,3), 'vivo', 'vitro'), 
             condition=factor(condition, levels = c('vitro', 'vivo')),
             strain=case_when(
               reporter %in% c(1,4) ~'strain1', 
               reporter %in% c(2,5) ~'strain2', 
               reporter %in% c(3,6) ~'strain3'
             )) %>% 
      group_by(accno) %>%
      mutate(VAR=var(intensity)) %>% 
      filter(VAR != 0) %>%#pull(accno) %>% unique()
      nest() %>% 
      mutate(TTEST=map(.x=data,
                       .f=~t.test(formula=intensity ~ condition, data=.x) %>% tidy())) %>% 
      select(-data) %>% 
      unnest(TTEST) %>% 
      ungroup() %>% 
      mutate(FDR=p.adjust(p.value, method = 'fdr')) %>% 
      arrange(estimate)
    return(TESTS)
  }


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
ALL_PEPTIDES <- read_tsv('maxquant_results/combined/txt/peptides.txt')

ALL_PEPTIDES <-
  ALL_PEPTIDES %>% 
  select(Sequence,`Leading razor protein`,  matches('Reporter intensity . [a-z]+')) %>% 
  mutate(Protein=`Leading razor protein`, 
         Peptide=Sequence) %>% 
  select(Protein, Peptide,everything(), -`Leading razor protein`, -Sequence) %>% 
  mutate(
    tmp_id=str_split(gsub('REV__tr','',Protein), pattern = '\\|'), 
    first_underscore=map_int(tmp_id, ~min(grep('_', .x)))) %>%  
  # select(first_underscore, tmp_id) %>% 
  mutate(TMP=map2_chr(tmp_id,first_underscore, ~pluck(.x,.y) ), 
         TMP=sub(';$','',TMP), 
         TMP=sub(';[trsp]+$','',TMP),
         TMP=sub('REV__', '', TMP), 
         Protein=TMP) %>% 
  select(-first_underscore, -tmp_id, -TMP) %>%
  filter(Protein %in% valid_accnos)



LONG_PEPTIDES <- 
  ALL_PEPTIDES %>%
  pivot_longer(cols = -c(Protein, Peptide), names_to = 'Sample', values_to = 'Abundance') %>% 
  mutate(iTRAQ=sub('Reporter intensity ([0-8]+) [a-z]+','\\1',Sample), 
         Run=sub('Reporter intensity ([0-8]+) ([a-z]+)','\\2',Sample),) %>% 
  select(-Sample) %>% 
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



LACT_TESTS <- NOMAD_diff_expression_tests(long_peptides = lact_peptides)
MAINT_TESTS <- NOMAD_diff_expression_tests(long_peptides = maint_peptides)


p_NOMAD_LACT_HIST <- 
  LACT_TESTS %>%
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
  MAINT_TESTS %>% 
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

tibble(strain=paste0('Strain', rep(c(1:3),4)), 
       condition=rep(c(rep('vivo', 3), rep('vitro', 3)),2), 
       diet=c(rep('Lact', 6), rep('maint', 6)), 
       LC_MSMS_run=c(rep('Run1',6), rep('Run2',6)), 
       iTRAQ_label=c(1:6, 1:6))
#

iTRAQ <- 
  PG %>% 
  select(accno, contains('Reporter intensity corrected')) %>% 
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
         condition=ifelse(reporter %in% c(1,2,3), 'vivo', 'vitro'), 
         strain=case_when(
           reporter %in% c(1,4) ~'strain1', 
           reporter %in% c(2,5) ~'strain2', 
           reporter %in% c(3,6) ~'strain3'
         )) %>% 
  ungroup() 
  
# iTRAQ %>% 
#   ggplot(aes(x=reporter, y=log_intensity, fill=condition)) + 
#   geom_col(color='black') + 
#   facet_wrap(~treatment)

# iTRAQ %>% 
#   ggplot(aes(x=reporter, y=rel_intensity, fill=condition)) + 
#   geom_col(color='black') + 
#   facet_wrap(~treatment)

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
  select(accno, strain, condition, intensity) %>% 
  pivot_wider(names_from = c(strain, condition), values_from = intensity ) %>% 
  column_to_rownames('accno') %>% as.matrix() %>% t()
  

iTRAQ_L_mat <- iTRAQ_L_mat[,colSums(iTRAQ_L_mat) != 0]
iTRAQ_L_mat_log <- apply(iTRAQ_L_mat, c(1,2), FUN = log)
iTRAQ_L_mat_log_ratio <- t(scale(t(iTRAQ_L_mat_log), scale=FALSE))
iTRAQ_L_mat_clr <- decostand(x = iTRAQ_L_mat, method = 'clr')
iTRAQ_L_mat_rab <- decostand(x = iTRAQ_L_mat, method = 'total')
iTRAQ_L_mat_log_rab <- decostand(x = iTRAQ_L_mat_log, method = 'total')

hist(iTRAQ_L_mat_log_rab)

raw_tests <- 
  run_tests(iTRAQ_L_mat) %>% 
  mutate(order=1:n(), 
         method='raw')


nomad_tests <- 
  LACT_TESTS %>% 
  mutate(order=1:n(), 
         method='nomad')

clr_tests <-
  run_tests(iTRAQ_L_mat_clr) %>%
  mutate(order=1:n(), 
         method='clr')

rab_tests <-  run_tests(iTRAQ_L_mat_rab) %>%
  mutate(order=1:n(), 
         method='rab')

log_rab_tests <-  run_tests(iTRAQ_L_mat_log_rab) %>%
  mutate(order=1:n(), 
         method='log_rab')

log_tests <- 
  run_tests(iTRAQ_L_mat_log) %>% 
  mutate(order=1:n(), 
         method='log')

#


all_tests <- 
  bind_rows(raw_tests, nomad_tests,clr_tests, rab_tests, log_rab_tests, log_tests) %>% 
  group_by(accno) %>% 
  mutate(MEAN_ORDER=mean(order), 
         ref_order=order[method=='clr']) %>% 
  mutate(SIG=ifelse(p.value < 0.05, T, F))

# 
# 
# all_tests <- 
#   bind_rows(raw_tests, nomad_tests,clr_tests, rab_tests, log_rab_tests, log_tests) %>% 
#   group_by(accno) %>% 
#   mutate(MEAN_ORDER=mean(order), 
#          ref_order=order[method=='nomad'])

#
all_tests %>% 
  ggplot(aes(x=order, y=ref_order, fill=SIG,)) + 
  geom_point(shape=21) + 
  facet_wrap(~method+SIG)

all_tests %>% 
  ggplot(aes(x=estimate, y=p.value, fill=SIG,)) + 
  geom_point(shape=21) + 
  facet_wrap(~method+SIG, scales = 'free')

all_tests %>% ggplot(aes(x=estimate, fill=SIG))+geom_histogram()+
  facet_wrap(~method, scales = 'free') + 
  geom_vline(xintercept = 0)


META <- 
  tibble(ID=rownames(iTRAQ_L_mat_clr)) %>% 
  separate(ID, into=c('strain', 'condition'), remove = F)


plot(metaMDS(iTRAQ_L_mat_clr, distance = 'euclidean', autotransform = FALSE))
NMDS <- metaMDS(iTRAQ_L_mat_clr, distance = 'euclidean', autotransform = FALSE)

prot_NMDS_coords <- NMDS$species %>% as.data.frame() %>% rownames_to_column(var='accno')


META <- 
  NMDS$points %>%
  as.data.frame() %>%
  rownames_to_column(var='ID') %>% 
  left_join(META)
library(ggrepel)

NMDS_plot <- 
  META %>%  
  ggplot(aes(x=MDS1, y=MDS2, fill=condition, label=strain)) +
  geom_point()+
  geom_label_repel() +
  theme_bw() + 
  lims(x=c(-8,8), y=c(-6,6))

adonis_results <- adonis2(method = 'euclidean',data = META, formula = iTRAQ_L_mat_clr ~ strain * condition)

# tst_clr <- 
#   iTRAQ_L_mat_clr %>%
#   as.data.frame() %>% 
#   rownames_to_column(var='ID') %>% 
#   separate(ID, into=c('strain', 'condition'),remove = FALSE) %>% 
#   pivot_longer(names_to = 'accno', values_to = 'intensity', -c(ID, strain, condition)) %>% 
#   group_by(accno) %>% 
#   nest() %>% 
#   mutate(TTEST=map(.x=data,
#                    .f=~t.test(formula=intensity ~ condition, data=.x) %>% tidy())) %>% 
#   select(-data) %>% 
#   unnest(TTEST) %>% 
#   ungroup() %>% 
#   mutate(FDR=p.adjust(p.value, method = 'fdr')) %>% 
#   arrange(p.value)
# 
# tst_log_ratios <- 
#   iTRAQ_L_mat_log_ratio %>%
#   as.data.frame() %>% 
#   rownames_to_column(var='ID') %>% 
#   separate(ID, into=c('strain', 'condition'),remove = FALSE) %>% 
#   pivot_longer(names_to = 'accno', values_to = 'intensity', -c(ID, strain, condition)) %>% 
#   group_by(accno) %>% 
#   nest() %>% 
#   mutate(TTEST=map(.x=data,
#                    .f=~t.test(formula=intensity ~ condition, data=.x) %>% tidy())) %>% 
#   select(-data) %>% 
#   unnest(TTEST) %>% 
#   ungroup() %>% 
#   mutate(FDR=p.adjust(p.value, method = 'fdr')) %>% 
#   arrange(p.value)
# 
# tst_raw_log <- 
#   iTRAQ_L_mat_log %>%
#   as.data.frame() %>% 
#   rownames_to_column(var='ID') %>% 
#   separate(ID, into=c('strain', 'condition'),remove = FALSE) %>% 
#   pivot_longer(names_to = 'accno', values_to = 'intensity', -c(ID, strain, condition)) %>% 
#   group_by(accno) %>% 
#   nest() %>% 
#   mutate(TTEST=map(.x=data,
#                    .f=~t.test(formula=intensity ~ condition, data=.x) %>% tidy())) %>% 
#   select(-data) %>% 
#   unnest(TTEST) %>% 
#   ungroup() %>% 
#   mutate(FDR=p.adjust(p.value, method = 'fdr')) %>% 
#   arrange(p.value)
# tst_raw_log$estimate %>% hist()


#########


iTRAQ_M <- iTRAQ %>% filter(treatment == 'maintenance')


iTRAQ_M_mat <- 
  iTRAQ_M %>% 
  select(accno, strain, condition, intensity) %>% 
  pivot_wider(names_from = c(strain, condition), values_from = intensity ) %>% 
  column_to_rownames('accno') %>% as.matrix() %>% t()



iTRAQ_M_mat <- iTRAQ_M_mat[,colSums(iTRAQ_M_mat) != 0]
library(vegan)
iTRAQ_M_mat_clr <- decostand(x = iTRAQ_M_mat, method = 'clr')


any(colSums(iTRAQ_L_mat) == 0)


adonis2(method = 'euclidean',data = META, formula = iTRAQ_M_mat_clr ~ strain + condition)





tst_clr %>% 
  mutate(FDR=p.adjust(p.value, 'fdr')) %>% 
  mutate(l2fc=log(estimate1/estimate2)) %>% 
  arrange(l2fc)

### 

# tests on raw log2 transformed intensities




####
# Relative abundance methods
# paired or unpaired T.test, doesnt really matter, very few sig diffs btw proteins in L
# none in M

tst <-
  iTRAQ_L %>% 
  select(-intensity, -reporter) %>% 
  pivot_wider(names_from = condition, values_from = rel_intensity ) %>% 
  group_by(accno) %>%
  nest() %>% 
  mutate(ttest=map(.x = data, ~t.test(.x$vitro, .x$vivo) %>% tidy())) %>% 
  select(accno, ttest) %>% 
  unnest(ttest) %>% 
  ungroup() %>% 
  filter(!is.nan(p.value)) %>% 
  mutate(FDR=p.adjust(p.value, 'fdr')) %>% 
  arrange(FDR)
  tst$method
# t.test(x=tst$data[[1]]$vivo, y=tst$data[[1]]$vitro, paired = T) %>% tidy()

iTRAQ_L %>% group_by(accno) %>% nest() %>% 
mutate(lms=map(.x = data, ~lm(data = .x, formula=rel_intensity~condition)), 
         res=map(.x=lms, .f=~tidy(.x, conf.int=T) %>% 
                   filter(term=='conditionvivo'))) %>% 
  select(accno, res) %>% 
  unnest(res) %>% 
  ungroup() %>% 
  filter(!is.nan(p.value)) %>% 
  mutate(FDR=p.adjust(p =p.value, method = 'fdr')) %>% 
  arrange((FDR))

iTRAQ_M_res <- 
  iTRAQ_M %>% 
  select(-intensity, -reporter) %>% 
  pivot_wider(names_from = condition, values_from = rel_intensity ) %>% 
  group_by(accno) %>%
  nest() %>% 
  mutate(ttest=map(.x = data, ~t.test(.x$vitro, .x$vivo, paired=T) %>% tidy())) %>% 
  select(accno, ttest) %>% 
  unnest(ttest) %>% 
  ungroup() %>% 
  filter(!is.nan(p.value)) %>% 
  mutate(FDR=p.adjust(p.value, 'fdr')) %>% 
  arrange(FDR)

######## multivariate analyses #########
# relative abundance, bray curtis dissimilarities
tst <- iTRAQ_L %>%
  select(accno, strain, condition, rel_intensity) %>% 
  pivot_wider(names_from = c(strain, condition), values_from = rel_intensity) %>% 
  column_to_rownames(var='accno') %>% as.matrix() %>% t()

tst[1:6, 1:3]

tst <- tst[,colSums(tst) != 0]

library(vegan)
library(ggrepel)
MDS_tst <- metaMDS(tst)
META <- MDS_tst$points %>% as.data.frame() %>% rownames_to_column(var="ID") %>% 
  separate(ID, into = c('strain', 'condition'))
META %>% 
  ggplot(aes(x=MDS1, y=MDS2)) + geom_point(aes(fill=condition), shape=21)+geom_text_repel(aes(label=strain))
vegan::vegdist(tst)

adonis2(data = META, formula = tst ~ condition + strain)



###

