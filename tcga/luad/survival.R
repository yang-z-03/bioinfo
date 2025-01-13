
files <- readRDS('tcga/luad/index.rds')[[1]]
clinical <- readRDS('tcga/luad/clinical.rds')
drugs <- readRDS('tcga/luad/drugs.rds')
radiations <- readRDS('tcga/luad/radiations.rds')
followups <- readRDS('tcga/luad/followups.rds')

dataset <- clinical[, c(
  'bcr_patient_barcode',
  'vital_status', 'days_to_birth',
  'days_to_death', 'days_to_last_followup',
  'person_neoplasm_cancer_status',
  'radiation_therapy', 'postoperative_rx_tx',
  'primary_therapy_outcome_success',
  'has_new_tumor_events_information',
  'has_radiations_information',
  'has_drugs_information',
  'stage_event_pathologic_stage',
  'stage_event_tnm_categories'
)]

colnames(dataset) <- c(
  'barcode', 'vital', 'age', 'death', 'last.followup',
  'status', 'rt.spon', 'chemo.spon', 'success.spon',
  'recurrence', 'rt.any', 'chemo.any', 'stage', 'tnm'
)

expr.records <- files[, c(
  'cases.submitter_id', 'sample_type', 'id', 'file_name'
)]

colnames(expr.records) <- c(
  'barcode', 'sample', 'id', 'file'
)

setdiff(dataset $ barcode, expr.records $ barcode)
expr.records $ sample |> table()
primary.records <- expr.records[expr.records $ sample == 'Primary Tumor', ]
setdiff(dataset $ barcode, primary.records $ barcode)

dataset $ barcode |> duplicated() |> sum()
primary.records $ barcode |> duplicated() |> sum()

# there are patients with multiple primary tumor samples for sequencing
# we just take the first one.

primary.records <- primary.records[!duplicated(primary.records $ barcode), ]
setdiff(dataset $ barcode, primary.records $ barcode)
setdiff(primary.records $ barcode, dataset $ barcode)

# there are 6 patients with no matching rna-seq. just ignore them.

meta.data <- merge(dataset, primary.records, by = 'barcode')

survival.t <- meta.data $ death
survival.t[is.na(survival.t)] <- meta.data $ last.followup[is.na(survival.t)]
meta.data $ t <- survival.t

meta.data $ vital <- 'alive'
meta.data $ vital[!is.na(meta.data $ death)] <- 'death'
meta.data $ vital[is.na(meta.data $ t)] <- '.'

# remove those not alive nor dead. :)

meta.data <- meta.data[meta.data $ vital != '.', ]

# we now have 437 cases.
# the next step is to construct the expression matrix.

meta.data $ expr <- paste('tcga/luad/gexp', meta.data $ id, meta.data $ file, sep = '/')
norm.counts <- list()

for (i in seq_along(meta.data $ expr)) {
  fpath <- meta.data $ expr[i]
  barcode <- meta.data $ barcode[i]
  expr.table <- read.table(fpath, header = TRUE, sep = '\t', comment.char = '#')
  
  norm.counts[['.id']] <- expr.table $ gene_id
  norm.counts[['.name']] <- expr.table $ gene_name
  norm.counts[[barcode]] <- expr.table $ fpkm_unstranded
}

fpkm <- as.data.frame(norm.counts)
fpkm <- fpkm[c(-1, -2, -3, -4), ]
rownames(fpkm) <- fpkm $ .name
saveRDS(fpkm, 'tcga/luad/fpkm.rds')

zdhhc15 <- fpkm[fpkm $ .name == 'ZDHHC15', ] |> t()
colnames(zdhhc15) <- 'ZDHHC15'
zdhhc15 <- zdhhc15[c(-1, -2), ]

length(zdhhc15)
# meta.data $ expr.zdhhc15 <- zdhhc15
meta.data $ expr.zdhhc15 <- as.numeric(zdhhc15)

saveRDS(meta.data, 'tcga/luad/metadata.rds')

# plotting survival analysis.

meta.data <- readRDS('tcga/luad/metadata.rds')
hist(meta.data $ expr.zdhhc15)
meta.data $ expr.log.zdhhc15 <- meta.data $ expr.zdhhc15 |> log1p()
hist(meta.data $ expr.log.zdhhc15)
median(meta.data $ expr.log.zdhhc15)

meta.data $ zdhhc.high <- 'ZDHHC15 High'
meta.data $ zdhhc.high[meta.data $ expr.log.zdhhc15 > 0.25] <- 'ZDHHC15 Low'

# status: censoring status 0 = censored, 1 = dead
meta.data $ censor <- 0
meta.data $ censor[meta.data $ vital == 'death'] <- 1
meta.data $ censor <- as.integer(meta.data $ censor)

# n = 437, overall, ungrouped survival data
s.overall <- ggsurvfit::survfit2(survival::Surv(t / 30, censor) ~ 1, data = meta.data)
ggsurvfit::ggsurvfit(s.overall) + ggplot2::labs(
  x = "Months",
  y = "Overall survival probability"
)

# overall, but in discrimination of expression data
s.expr <- ggsurvfit::survfit2(survival::Surv(t / 30, censor) ~ zdhhc.high, data = meta.data)
ggsurvfit::ggsurvfit(s.expr) +
  ggplot2::labs(
    x = "Months", y = "Overall survival probability"
  ) + ggplot2::theme(legend.position = c(0.7, 0.85)) +
  ggsurvfit::add_confidence_interval() +
  ggsurvfit::add_risktable()

# mixed groupings with chemotherapy and radiotherapy status.

meta.data $ state.chemo <- 'No chemotherapy'
meta.data $ state.chemo[meta.data $ chemo.any == 'YES'] <- 'Chemotherapy'
meta.data $ state.rt <- 'No radiation'
meta.data $ state.rt[meta.data $ rt.any == 'YES'] <- 'Radiation'

meta.data $ zdhhc.chemo <- paste(meta.data $ zdhhc.high, meta.data $ state.chemo, sep = ', ')
meta.data $ zdhhc.rt <- paste(meta.data $ zdhhc.high, meta.data $ state.rt, sep = ', ')

s.zdhhc.chemo <- ggsurvfit::survfit2(survival::Surv(t / 30, censor) ~ zdhhc.chemo, data = meta.data)
ggsurvfit::ggsurvfit(s.zdhhc.chemo) +
  ggplot2::labs(
    x = "Months", y = "Overall survival probability"
  ) + ggplot2::theme(legend.position = 'top') +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) + 
  ggsurvfit::add_confidence_interval()

s.zdhhc.rt <- ggsurvfit::survfit2(survival::Surv(t / 30, censor) ~ zdhhc.rt, data = meta.data)
ggsurvfit::ggsurvfit(s.zdhhc.rt) +
  ggplot2::labs(
    x = "Months", y = "Overall survival probability"
  ) + ggplot2::theme(legend.position = 'top') +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) + 
  ggsurvfit::add_confidence_interval()

s1 <- survminer::surv_fit(
  survival::Surv(t / 30, censor) ~ zdhhc.high, 
  group.by = 'state.chemo', data = meta.data
)

survminer::surv_pvalue(s1)