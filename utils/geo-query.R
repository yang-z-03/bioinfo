
library(GEOquery)
library(tidyverse)
library(limma)
library(magrittr)
library(dplyr)

## fetching information from GEO ----------------------------------------------

# here we provide several methods to easily download information from geo
# database. geo database contains the following level of objects, both have an
# unique identifier starting with:

# platforms (GPL#)    a platform record is composed of a summary description of
#                     the **array** or **sequencer** and, for array-based
#                     Platforms, a data table defining the array template. Each
#                     Platform record is assigned a unique and stable GEO
#                     accession number (GPLxxx). A Platform may reference many
#                     Samples that have been submitted by multiple submitters.

# samples (GSM#)      a sample record describes the conditions under which an
#                     individual Sample was handled, the manipulations it
#                     underwent, and the abundance measurement of each element
#                     derived from it. Each Sample record is assigned a unique
#                     and stable GEO accession number (GSMxxx). A Sample entity
#                     must reference only one Platform and may be included in
#                     multiple Series.

# series (GSE#)       a series record links together a group of related Samples
#                     and provides a focal point and description of the whole
#                     study. Series records may also contain tables describing
#                     extracted data, summary conclusions, or analyses. Each
#                     Series record is assigned a unique and stable GEO
#                     accession number (GSExxx).

# advanced summaries of original experiments

# datasets (GDS#)     a dataset represents a curated collection of biologically
#                     and statistically comparable GEO Samples and forms the
#                     basis of GEO"s suite of data display and analysis tools.
#
#                     samples within a DataSet refer to the same Platform, that
#                     is, they share a common set of array elements. Value
#                     measurements for each Sample within a DataSet are assumed
#                     to be calculated in an equivalent manner, that is,
#                     considerations such as background processing and
#                     normalization are consistent across the DataSet.
#                     Information reflecting experimental factors is provided
#                     through DataSet subsets. Not all submitted data are
#                     suitable for DataSet assembly and we are experiencing a
#                     backlog in DataSet creation, so not all Series have
#                     corresponding DataSet record(s).

# profiles            a profile consists of the expression measurements for an
#   (gene @ GDS#)     individual gene across all Samples in a DataSet.


## downloading data from GEO --------------------------------------------------

# 联网下载 GSE 的补充文件到指定的文件夹（默认：geo）
geo_get_supplemental_files <- function(index,
                                       base_dir = "geo",
                                       make_dir = TRUE,
                                       filter_regex = NULL) {
  GEOquery::getGEOSuppFiles(index,
    baseDir = base_dir,
    makeDirectory = make_dir,
    filter_regex = filter_regex
  )
}

# 只联网获得 GSE 的补充文件名称和网络位置以供自行下载
# 返回一个两列的表：第一列，文件名；第二列，URL
geo_get_supplemental_links <- function(index) {
  url <- lapply(index, function(x) {
    GEOquery::getGEOSuppFiles(x, makeDirectory = FALSE,
                              fetch_files = FALSE)
  }) |> bind_rows()

  # 构造一个两列的表：["fname", "url"]
  return(url)
}

## grouping info and platform info --------------------------------------------

# 使用 GSE 号获得系列信息：这个函数一般会在 geo 文件夹下下载三个文件：
# - gse***_series_matrix.txt.gz
# - gpl***.annot.gz
# - gpl***.soft.gz
#
# 系列信息是一个列表，其中通常只有一个元素，名字就是 .txt.gz 的文件名，类型 ExpressionSet
# ExpressionSet 是一个很大的 S4 对象，他有如下五个属性：
#
# .@ experimentData <MIAME>
#    包含作者信息、实验室信息、地址、简介、标题、PMID 等文献信息
#
# .@ assayData @ exprs
#    原始数据，一个巨大的 double 数据表，行数不定，是芯片测试的基因数，这个 GSE 中有多少
#    个样本，就有多少列。行名称和列名称是字符串 [ngenes * nsample]
#
# .@ phenoData @ data
#    描述分组信息的数据表，一行是一个样本，列名实例如下：
#    "title"                   "geo_accession"           "status"
#    "submission_date"         "last_update_date"        "type"
#    "channel_count"           "source_name_ch1"         "organism_ch1"
#    "characteristics_ch1"     "characteristics_ch1.1"   "characteristics_ch1.2"
#    "treatment_protocol_ch1"  "growth_protocol_ch1"     "molecule_ch1"
#    "extract_protocol_ch1"    "label_ch1"               "label_protocol_ch1"
#    "taxid_ch1"               "hyb_protocol"            "scan_protocol"
#    "description"             "data_processing"         "platform_id"
#    "contact_name"            "contact_email"           "contact_phone"
#    "contact_institute"       "contact_address"         "contact_city"
#    "contact_state"           "contact_zip/postal_code" "contact_country"
#    "supplementary_file"      "data_row_count"          "cell line:ch1"
#    "cell origin:ch1"         "genotype/variation:ch1"
#
# .@ featureData @ data
#    表示这个芯片中测定的基因信息，一行一个基因，列名实例如下：
#    "ID"                               "GB_ACC"
#    "SPOT_ID"                          "Species Scientific Name"
#    "Annotation Date"                  "Sequence Type"
#    "Sequence Source"                  "Target Description"
#    "Representative Public ID"         "Gene Title"
#    "Gene Symbol"                      "ENTREZ_GENE_ID"
#    "RefSeq Transcript ID"             "Gene Ontology Biological Process"
#    "Gene Ontology Cellular Component" "Gene Ontology Molecular Function"
#
# .@ annotation
#    这个系列所用到的平台的 GPL 编号
#
# .@ protocolData @ data 一般为空
geo_get_series_annot <- function(series_id, base_dir = "geo") {
  gpl <- GEOquery::getGEO(series_id, destdir = base_dir)

  exprset_names <- names(gpl)
  if (length(exprset_names) == 1) return(gpl[[exprset_names]])
  else return()
}

# 从输入的 ExpressionSet 对象中自动提取并下载 GPL 对象
# 返回一个很大的 GPL 对象，它的常用属性如下表：
#
# .@ header
#    .$ date                 这个 GPL 登记的日期
#    .$ platform             平台的 GPL 号
#    .$ platform_organism    是给那个生物的基因制作的
#    .$ platform_title       公司给出的全名
#
# .@ dataTable <GEODataTable>
#    .@ columns <data.frame>
#       对于自己给出的基因信息列的说明
#
#       1                     ID         ID from Platform data table
#       2             Gene title                    Entrez Gene name
#       3            Gene symbol                  Entrez Gene symbol
#       4                Gene ID              Entrez Gene identifier
#       5          UniGene title                 Entrez UniGene name
#       6         UniGene symbol               Entrez UniGene symbol
#       7             UniGene ID           Entrez UniGene identifier
#       8       Nucleotide Title             Entrez Nucleotide title
#       9                     GI                  GenBank identifier
#       10     GenBank Accession                   GenBank accession
#       11      Platform_CLONEID   CLONE_ID from Platform data table
#       12          Platform_ORF        ORF from Platform data table
#       13       Platform_SPOTID    SPOT_ID from Platform data table
#       14   Chromosome location Entrez gene chromosome and location
#       15 Chromosome annotation   Entrez gene chromosome annotation
#       16           GO:Function         Gene Ontology Function term
#       17            GO:Process          Gene Ontology Process term
#       18          GO:Component        Gene Ontology Component term
#       19        GO:Function ID   Gene Ontology Function identifier
#       20         GO:Process ID    Gene Ontology Process identifier
#       21       GO:Component ID  Gene Ontology Component identifier
#
#    .@ table <data.frame>
#       一行一个基因，每个基因都有上面所示的那些列
geo_get_platform_annot <- function(expr_set, base_dir = "geo") {
  gpl <- GEOquery::getGEO(expr_set @ annotation,
                          destdir = base_dir,
                          AnnotGPL = TRUE)
  return(gpl)
}

# 从 GSE 的信息表中提取实验分组信息
geo_get_groupings <- function(expr_set,
                              sel_cols = c("geo_accession", "title",
                                           "type", "source_name_ch1",
                                           "taxid_ch1", "platform_id")) {
  return(expr_set @ phenoData @ data[, sel_cols, drop = FALSE])
}