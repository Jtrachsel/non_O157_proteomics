library(tidyverse)

# iBAQ
PG <- read_tsv('output/protein_groups_cleaned.tsv')
prot_descripts <- read_tsv('output/protein_descriptions.tsv')
PG$Peptides

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
  arrange(desc(l2FC)) %>% 
  left_join(prot_descripts)

up_L_iBAQ <- nrow(sig_L_iBAQ)
up_L_iBAQ

sig_M_iBAQ <- 
  riBAQ_res %>% 
  filter(l2FC < -1) %>% 
  arrange(l2FC) %>% 
  left_join(prot_descripts)

sig_M_iBAQ %>% filter(both)

up_M_iBAQ <- nrow(sig_M_iBAQ)
up_M_iBAQ

num_only_L <- riBAQ_res$only_L %>% sum()
num_only_M <- riBAQ_res$only_M %>% sum()
