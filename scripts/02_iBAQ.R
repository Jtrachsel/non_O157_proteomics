library(tidyverse)

# iBAQ
PG <- read_tsv('output/protein_groups_cleaned.tsv')
prot_descripts <- read_tsv('output/protein_descriptions.tsv')
# PG$Peptides

iBAQ_df <- 
  PG %>% 
  filter(!is.nan(iBAQ)) %>% 
  filter(iBAQ != 0) %>% 
  filter(Peptides > 1) %>%
  filter(!(grepl('CON', accno))) %>% 
  select(accno, contains('iBAQ')) %>%
  transmute(accno, 
            iBAQ_L=`iBAQ lactation`, 
            iBAQ_M=`iBAQ maintenance`) %>% 
  mutate(riBAQ_L=iBAQ_L/sum(iBAQ_L), 
         riBAQ_M=iBAQ_M/sum(iBAQ_M)) %>% 
  pivot_longer(names_to = c('type', 'condition'), names_sep = '_', values_to = 'value', cols = -c(accno))

p1 <- 
  iBAQ_df %>% 
  ggplot(aes(x=condition, y=value, fill=condition)) + 
  geom_col(color='black') + 
  facet_wrap(vars(type), scales = 'free') + 
  ggtitle('Figure 1: Raw iBAQ and riBAQ normalized intensities')


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


# THIS ONE
p2 <- 
  riBAQ_res %>%
  ggplot(aes(x=l2FC, fill=both)) + 
  geom_histogram(bins=50)+
  lims(x=c(-17,17), 
       y=c(0,100)) +
  annotate(x=9, y=50, geom='label', label='Enriched in Lactation')+
  annotate(x=-9, y=50, geom='label', label='Enriched in Maintenance')+
  geom_vline(xintercept = 0)+
  labs(fill='Detected in both diets', 
       y='count') +
  theme(legend.position = 'bottom')+
  ggtitle('Figure 2: Histogram of all log2FoldChange values')


# THIS ONE SIGS
p3 <- riBAQ_res %>% 
  filter(abs(l2FC) > 1) %>%
  arrange((l2FC)) %>% 
  left_join(prot_descripts) %>% 
  ggplot(aes(x=l2FC, fill=both)) + 
  annotate(x=9, y=50, geom='label', label='Enriched in Lactation')+
  annotate(x=-9, y=50, geom='label', label='Enriched in Maintenance')+
  geom_vline(xintercept = 0)+
  geom_histogram(bins=50)+
  lims(x=c(-17,17), 
       y=c(0,100))+
  labs(fill='Detected in both diets', 
       y='count') +
  theme(legend.position = 'bottom')+
  ggtitle('Figure 3: Histogram of log2FoldChange values  > 1')



sig_L_iBAQ <- 
  riBAQ_res %>% 
  filter(l2FC > 1) %>% 
  arrange(desc(l2FC)) %>% 
  left_join(prot_descripts)

up_L_iBAQ <- nrow(sig_L_iBAQ)

sig_M_iBAQ <- 
  riBAQ_res %>% 
  filter(l2FC < -1) %>% 
  arrange(l2FC) %>% 
  left_join(prot_descripts)


up_M_iBAQ <- nrow(sig_M_iBAQ)


##
num_not_not_diff <- riBAQ_res %>% 
  filter(abs(l2FC) < 1) %>% nrow()

#
num_only_L <- riBAQ_res$only_L %>% sum()
num_only_M <- riBAQ_res$only_M %>% sum()
num_both <- riBAQ_res$both %>% sum()
total <- nrow(riBAQ_res)

# num proteins detected
T1 <- tribble(~category, ~ 'number of proteins', 
        'total',total, 
        'both diets',num_both,
        'lactation only',num_only_L,
        'maintenance only',num_only_M)


# 'sig' different proteins
T2 <- tribble(~category, ~'number of proteins', 
        'not different', num_not_not_diff, 
        'up in lactation',up_L_iBAQ, 
        'up in maintenance',up_M_iBAQ)

# membrane prots up in L  
T3 <- sig_L_iBAQ %>% filter(both)%>%
  filter(grepl('Mem', Localization))%>% 
  select(accno, l2FC, description, Localization)


T4 <- sig_L_iBAQ %>% 
  select(accno, l2FC, description, Localization)

# membrane prots up in M  
T5 <- sig_M_iBAQ %>% filter(both) %>%
  filter(grepl('Mem', Localization)) %>% 
  select(accno, l2FC, description, Localization)

T6 <- sig_M_iBAQ %>% 
  select(accno, l2FC, description, Localization)





exp_design_tibble <- 
  tibble(strain=paste0('Strain', rep(c(1:3),4)), 
         condition=rep(c(rep('vivo', 3), rep('vitro', 3)),2), 
         diet=c(rep('Lact', 6), rep('maint', 6)), 
         LC_MSMS_run=c(rep('Run1',6), rep('Run2',6)), 
         iTRAQ_label=c(1:6, 1:6))
#




