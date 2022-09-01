
#####
# 
# PG %>%
#   filter(!grepl('CON', x = `Protein IDs`)) %>% 
#   ggplot(aes(x=log(`Reporter intensity 1 lactation`),
#              y=log(`Reporter intensity corrected 1 lactation`))) + 
#   geom_point()



# iBAQ
PG <- read_tsv('output/protein_groups_cleaned.tsv')

iBAQ_df <- 
  PG %>% 
  filter(!is.nan(iBAQ)) %>% 
  filter(iBAQ != 0) %>% 
  filter(!(grepl('CON', accno))) %>% 
  select(accno, contains('iBAQ')) %>%
  transmute(accno, 
            iBAQ_L=`iBAQ lactation`, 
            iBAQ_M=`iBAQ maintenance`) %>% 
  mutate(riBAQ_L=iBAQ_L/sum(iBAQ_L), 
         riBAQ_M=iBAQ_M/sum(iBAQ_M)) %>% 
  pivot_longer(names_to = c('type', 'condition'), names_sep = '_', values_to = 'value', cols = -c(accno))

iBAQ_df %>% 
  ggplot(aes(x=condition, y=value, fill=condition)) + 
  geom_col(color='black') + 
  facet_wrap(vars(type), scales = 'free')
  

riBAQ_df <- iBAQ_df %>% filter(type == 'riBAQ')

pseudo_val <- min(riBAQ_df$value[riBAQ_df$value !=0]) / 2

# riBAQ_df %>% 
#   filter(value < .0001) %>%
#   ggplot(aes(x=(value))) + 
#   geom_histogram(bins = 50)


riBAQ_res <- 
  riBAQ_df %>% 
  group_by(accno) %>% 
  mutate(only_L=value[condition == 'L'] > 0 & value[condition == 'M'] == 0, 
         only_M=value[condition == 'M'] > 0 & value[condition == 'L'] == 0, 
         both=only_L == only_M) %>% 
  ungroup() %>% 
  mutate(value=value + pseudo_val) %>% 
  group_by(accno) %>% 
  summarise(only_L=unique(only_L), only_M=unique(only_M),both=unique(both), 
            l2FC=log2(value[condition=='L']/value[condition== 'M']))


riBAQ_res %>%
  ggplot(aes(x=l2FC, fill=both)) + 
  geom_histogram(bins=50)+
  lims(x=c(-17,17), 
       y=c(0,100))

riBAQ_res %>% 
  filter(abs(l2FC) > 1) %>%
  arrange((l2FC)) %>% 
  left_join(prot_descripts) %>% 
  ggplot(aes(x=l2FC, fill=both)) + 
  geom_histogram(bins=50)+
  lims(x=c(-17,17), 
       y=c(0,100))


sig_L_iBAQ <- 
  riBAQ_res %>% 
  filter(l2FC > 1) %>% 
  left_join(prot_descripts)

up_L_iBAQ <- nrow(sig_L_iBAQ)
up_L_iBAQ

sig_M_iBAQ <- 
  riBAQ_res %>% 
  filter(l2FC < -1) %>% 
  left_join(prot_descripts)

up_M_iBAQ <- nrow(sig_M_iBAQ)

num_only_L <- riBAQ_res$only_L %>% sum()
num_only_M <- riBAQ_res$only_M %>% sum()

########

##########


tibble(strain=paste0('Strain', rep(c(1:3),4)), 
       condition=rep(c(rep('vivo', 3), rep('vitro', 3)),2), 
       diet=c(rep('Lact', 6), rep('maint', 6)), 
       LC_MSMS_run=c(rep('Run1',6), rep('Run2',6)), 
       iTRAQ_label=c(1:6, 1:6))


iTRAQ <- 
  PG %>% 
  select(accno, contains('Reporter intensity corrected')) %>% 
  pivot_longer(-accno, names_to = c('reporter', 'treatment'),
               names_prefix = 'Reporter intensity corrected ', 
               names_sep = ' ', values_to = 'intensity')

#  itraq reporter intensities
iTRAQ %>% 
  ggplot(aes(x=reporter, y=intensity, fill=grepl('CON', accno))) + 
  geom_col(color='black') + 
  facet_wrap(~treatment)



iTRAQ <- 
  iTRAQ %>%
  filter(!grepl('CON', accno)) %>% 
  filter(reporter %in% c(1:6)) %>% 
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
  


riBAQ_res %>% filter(both)

# this shows us that within each MSMS run we only dectect proteins that
# were detectable in all strains across both in-vivo and in-vitro conditions.
iTRAQ %>%
  group_by(treatment, accno) %>% 
  summarise(ANY=any(intensity == 0), 
            ALL=all(intensity == 0)) %>% 
  ungroup() %>% 
  mutate(ANY_NOT_ALL=ANY & !ALL) %>% 
  filter(ANY_NOT_ALL)


iTRAQ %>%
  ggplot(aes(x=reporter, y=rel_intensity, fill=condition)) + 
  geom_col(color='black') + 
  facet_wrap(~treatment)

iTRAQ %>%
  ggplot(aes(x=reporter, y=log_intensity, fill=condition)) + 
  geom_col(color='black') + 
  facet_wrap(~treatment)

iTRAQ_L <-
  iTRAQ %>% 
  filter(treatment == 'lactation' & is.finite(log_intensity)) %>% 
  group_by(reporter) %>% 
  mutate(norm_log_intensity=log_intensity - mean(log_intensity), 
         log_rel_intensity=log(rel_intensity)) %>% 
  ungroup()



iTRAQ_L %>%
  ggplot(aes(x=reporter, y=norm_log_intensity, fill=condition)) + 
  geom_col(color='black') + 
  facet_wrap(~treatment)



iTRAQ_L$log_intensity

iTRAQ_L_mat <- 
  iTRAQ_L %>% 
  select(-rel_intensity,-log_intensity, -norm_log_intensity, -log_rel_intensity, -reporter, -treatment) %>% 
  pivot_wider(names_from = c(strain, condition), values_from = intensity ) %>% 
  column_to_rownames('accno') %>% as.matrix() %>% t()
  

iTRAQ_L_mat <- iTRAQ_L_mat[,colSums(iTRAQ_L_mat) != 0]




iTRAQ_M <- iTRAQ %>% filter(treatment == 'maintenance')


iTRAQ_M_mat <- 
  iTRAQ_M %>% 
  select(-rel_intensity,-log_intensity, -reporter, -treatment) %>% 
  pivot_wider(names_from = c(strain, condition), values_from = intensity ) %>% 
  column_to_rownames('accno') %>% as.matrix() %>% t()


iTRAQ_M_mat <- iTRAQ_M_mat[,colSums(iTRAQ_M_mat) != 0]
library(vegan)
iTRAQ_L_mat_clr <- decostand(x = iTRAQ_L_mat, method = 'clr')
iTRAQ_M_mat_clr <- decostand(x = iTRAQ_M_mat, method = 'clr')


any(colSums(iTRAQ_L_mat) == 0)

META <- 
  tibble(ID=rownames(iTRAQ_L_mat_clr)) %>% 
  separate(ID, into=c('strain', 'condition'), remove = F)


any(colSums(iTRAQ_L_mat) == 0)
plot(metaMDS(iTRAQ_L_mat, distance = 'robust.aitchison', autotransform = FALSE))



adonis2(method = 'euclidean',data = META, formula = iTRAQ_L_mat_clr ~ strain + condition)
adonis2(method = 'euclidean',data = META, formula = iTRAQ_M_mat_clr ~ strain + condition)

tst_clr <- 
  iTRAQ_L_mat_clr %>%
  as.data.frame() %>% 
  rownames_to_column(var='ID') %>% 
  separate(ID, into=c('strain', 'condition'),remove = FALSE) %>% 
  pivot_longer(names_to = 'accno', values_to = 'intensity', -c(ID, strain, condition)) %>% 
  group_by(accno) %>% 
  nest() %>% 
  mutate(TTEST=map(.x=data,
                .f=~t.test(formula=intensity ~ condition, data=.x) %>% tidy())) %>% 
  select(-data) %>% 
  unnest(TTEST) %>% 
  ungroup() 



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
t.test(x=tst$data[[1]]$vivo, y=tst$data[[1]]$vitro, paired = T) %>% tidy()

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

t.test()


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

