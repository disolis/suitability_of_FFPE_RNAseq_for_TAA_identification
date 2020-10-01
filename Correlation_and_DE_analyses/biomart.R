library("biomaRt")
listMarts()

ensembl = useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)

ensembl = useDataset("clfamiliaris_gene_ensembl",mart=ensembl)
BM <- getBM(attributes=c('ensembl_gene_id','external_gene_name', 'gene_biotype','description'),
      mart = ensembl)

