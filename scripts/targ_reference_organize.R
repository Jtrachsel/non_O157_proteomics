library(tidyverse)

# got this by mapping uniparc_IDs from mmseqs2 clustered rep sequences
KB_map <- read_tsv('./reference/targeted_search_database/map_to_uniprot_KB.tsv')


# take the first mapping result for each rep protein
KB_map %>%
  group_by(From) %>% 
  slice_head(n=1) %>% 
  ungroup() %>% 
  write_tsv('reference/targeted_search_database/uniparc_to_uniprot.tsv') %>% 
  pull(To) %>% 
  write_lines('reference/targeted_search_database/uniprot_IDs.txt')

### RUN UNIPROT THINGS HERE ###

###

uniparc_to_uniprot <- read_tsv('reference/targeted_search_database/uniparc_to_uniprot.tsv') %>% 
  transmute(uniparc_ID=From, 
            uniprot_ID=To)

uniprot_annotations <- read_tsv('reference/targeted_search_database/targeted_annotations.tsv')%>% 
  mutate(uniprot_ID=From) %>% 
  select(-From)


mapped_annotations <- uniparc_to_uniprot %>% left_join(uniprot_annotations)


library(Biostrings)


FASTA <- readAAStringSet('reference/targeted_search_database/targ_clustered_rep_seq.fasta')

names(FASTA) <- sub('(UPI.*) status=.*','\\1', names(FASTA))
writeXStringSet(FASTA,filepath = 'reference/targeted_search_database.faa')

mapped_annotations
