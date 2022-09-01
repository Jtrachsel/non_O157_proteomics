
#!/bin/bash
set -e

echo "download reference from uniprot"
echo "This proteome is part of the Escherichia coli (strain K12) pan proteome"
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/UP000000625.fasta.gz

mkdir reference
gunzip UP000000625.fasta.gz && mv UP000000625.fasta ./reference

# extract just the contig IDs for downstream protein identification
cat reference/UP000000625.fasta |grep '>' > ./reference/IDs.txt

