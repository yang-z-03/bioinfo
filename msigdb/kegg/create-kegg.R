
remotes::install_github('YuLab-SMU/createKEGGdb')
setwd('~/bioinfo/msigdb/kegg')

createKEGGdb::create_kegg_db('hsa')
