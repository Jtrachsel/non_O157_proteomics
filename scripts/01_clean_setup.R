library(tidyverse)
library(broom)
library(Biostrings)
usethis::use_directory('output')

# read in maxquant file
PG <- read_tsv('data/proteinGroups.txt')


sum(PG$`Unique peptides` > 1)
# library(IHW)

# BiocManager::install('IHW')

# overview #



# Identify all detected proteins
PG <- 
  PG %>% 
  mutate(
    tmp_id=str_split(gsub('REV__tr','',`Majority protein IDs`), pattern = '\\|'), 
    first_underscore=map_int(tmp_id, ~min(grep('_', .x)))) %>%  
  # select(first_underscore, tmp_id) %>% 
  mutate(TMP=map2_chr(tmp_id,first_underscore, ~pluck(.x,.y) ), 
         TMP=sub(';$','',TMP), 
         TMP=sub(';[trsp]+$','',TMP),
         TMP=sub('REV__', '', TMP), 
         accno=TMP) %>% #pull(accno)
  select(accno, everything(),-first_underscore, -tmp_id, -TMP)

PG %>% write_tsv('output/protein_groups_cleaned.tsv')

# for independent hypothesis weighting
IHW_dat <- 
  PG %>%
  filter(!grepl('CON', accno)) %>% 
  transmute(accno, unique_peptides=`Unique peptides`) %>% 
  arrange(desc(unique_peptides))


### protein accession and annotation 

prot_descripts <- 
  read_tsv('reference/IDs.txt', col_names = c('raw_ID')) %>% 
  mutate(accno=      sub('.*\\|([[:alnum:]]+)\\|([[:alnum:]]+_[[:alnum:]]+) (.*)','\\2',raw_ID), 
         description=sub('.*\\|([[:alnum:]]+)\\|([[:alnum:]]+_[[:alnum:]]+) (.*)','\\3',raw_ID)) %>% 
  select(-raw_ID)


prot_descripts <- 
  prot_descripts %>% 
  filter(accno %in% PG$accno)



ref_fasta <- readAAStringSet('reference/UP000000625.fasta')

names(ref_fasta) <-
  sub('.*\\|([[:alnum:]]+)\\|([[:alnum:]]+_[[:alnum:]]+) (.*)',
      '\\2',
      names(ref_fasta))

ref_fasta <- ref_fasta[prot_descripts$accno]

ref_fasta %>%
  writeXStringSet('reference/detected_proteins.faa')
names(ref_fasta) %>% 
  write_lines('output/detected_accnos.txt')

### psort

psort <- 
  read_tsv('output/psort_res.tsv') %>% 
  transmute(accno=SeqID, 
            Localization)

prot_descripts <- prot_descripts %>% left_join(psort)


prot_descripts %>% write_tsv('output/protein_descriptions.tsv')
