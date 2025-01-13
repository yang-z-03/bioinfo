
filter <- 'TCGA-LUAD'

clinical <- TCGAbiolinks::GDCquery_clinic(
  project = filter,
  type = 'clinical'
)

query.clinic <- TCGAbiolinks::GDCquery(
  project = filter, 
  data.category = "Clinical", 
  data.type = "Clinical Supplement",
  data.format = "bcr xml"
)

TCGAbiolinks::GDCdownload(
  query.clinic, method = "api", files.per.chunk = 50,
  directory = 'tcga/luad/xmls'
)

xml.patient <- TCGAbiolinks::GDCprepare_clinic(
  query.clinic, clinical.info = 'patient',
  directory = 'tcga/luad/xmls'
)

xml.drug <- TCGAbiolinks::GDCprepare_clinic(
  query.clinic, clinical.info = 'drug',
  directory = 'tcga/luad/xmls'
)

xml.admin <- TCGAbiolinks::GDCprepare_clinic(
  query.clinic, clinical.info = 'admin',
  directory = 'tcga/luad/xmls'
)

xml.followup <- TCGAbiolinks::GDCprepare_clinic(
  query.clinic, clinical.info = 'follow_up',
  directory = 'tcga/luad/xmls'
)

xml.radiation <- TCGAbiolinks::GDCprepare_clinic(
  query.clinic, clinical.info = 'radiation',
  directory = 'tcga/luad/xmls'
)

xml.stage <- TCGAbiolinks::GDCprepare_clinic(
  query.clinic, clinical.info = 'stage_event',
  directory = 'tcga/luad/xmls'
)

xml.tumor <- TCGAbiolinks::GDCprepare_clinic(
  query.clinic, clinical.info = 'new_tumor_event',
  directory = 'tcga/luad/xmls'
)

# expression profiles:
# we will use the expression profiles only to download these files
# and parse the tables ourself to construct an expression matrix.

query <- TCGAbiolinks::GDCquery(
  project = filter, 
  data.category = "Transcriptome Profiling", 
  data.type = "Gene Expression Quantification"
)

TCGAbiolinks::GDCdownload(
  query, method = "api", files.per.chunk = 50,
  directory = 'tcga/luad/gexp'
)

saveRDS(query $ results, file = 'tcga/luad/index.rds')

# merging the xml for complete record.

is.unique <- function(x) {
  duplicated(x) |> sum() == 0
}

is.unique(xml.admin $ bcr_patient_barcode)
is.unique(xml.patient $ bcr_patient_barcode)
is.unique(xml.drug $ bcr_patient_barcode) # false
is.unique(xml.radiation $ bcr_patient_barcode) # false
is.unique(xml.stage $ bcr_patient_barcode)
is.unique(xml.tumor $ bcr_patient_barcode)
is.unique(xml.followup $ bcr_patient_barcode) # false

m <- merge(xml.admin, xml.patient, by = 'bcr_patient_barcode')
m <- merge(m, xml.stage, by = 'bcr_patient_barcode')
m <- merge(m, xml.tumor, by = 'bcr_patient_barcode')

# merged data:
# 
#  [1] "bcr_patient_barcode"                            
#  [2] "bcr"                                            
#  [3] "file_uuid"                                      
#  [4] "batch_number"                                   
#  [5] "project_code"                                   
#  [6] "disease_code"                                   
#  [7] "day_of_dcc_upload"                              
#  [8] "month_of_dcc_upload"                            
#  [9] "year_of_dcc_upload"                             
# [10] "patient_withdrawal"                             
# [11] "program"                                        
# [12] "dbgap_registration_code"                        
# [13] "additional_studies"                             
# [14] "tumor_tissue_site"                              
# [15] "histological_type"                              
# [16] "other_dx"                                       
# [17] "gender"                                         
# [18] "vital_status"                                   
# [19] "days_to_birth"                                  
# [20] "days_to_last_known_alive"                       
# [21] "days_to_death"                                  
# [22] "days_to_last_followup"                          
# [23] "race_list"                                      
# [24] "tissue_source_site"                             
# [25] "patient_id"                                     
# [26] "bcr_patient_uuid"                               
# [27] "history_of_neoadjuvant_treatment"               
# [28] "informed_consent_verified"                      
# [29] "icd_o_3_site"                                   
# [30] "icd_o_3_histology"                              
# [31] "icd_10"                                         
# [32] "tissue_prospective_collection_indicator"        
# [33] "tissue_retrospective_collection_indicator"      
# [34] "days_to_initial_pathologic_diagnosis"           
# [35] "age_at_initial_pathologic_diagnosis"            
# [36] "year_of_initial_pathologic_diagnosis"           
# [37] "ethnicity"                                      
# [38] "person_neoplasm_cancer_status"                  
# [39] "performance_status_scale_timing"                
# [40] "day_of_form_completion"                         
# [41] "month_of_form_completion"                       
# [42] "year_of_form_completion"                        
# [43] "karnofsky_performance_score"                    
# [44] "eastern_cancer_oncology_group"                  
# [45] "tobacco_smoking_history"                        
# [46] "year_of_tobacco_smoking_onset"                  
# [47] "stopped_smoking_year"                           
# [48] "number_pack_years_smoked"                       
# [49] "anatomic_neoplasm_subdivision"                  
# [50] "anatomic_neoplasm_subdivision_other"            
# [51] "diagnosis"                                      
# [52] "location_in_lung_parenchyma"                    
# [53] "residual_tumor"                                 
# [54] "kras_mutation_found"                            
# [55] "kras_gene_analysis_performed"                   
# [56] "kras_mutation_result"                           
# [57] "egfr_mutation_performed"                        
# [58] "egfr_mutation_identified"                       
# [59] "egfr_mutation_result"                           
# [60] "eml4_alk_translocation_performed"               
# [61] "eml4_alk_translocation_identified"              
# [62] "eml4_alk_translocation_result"                  
# [63] "eml4_alk_translocation_method"                  
# [64] "pulmonary_function_test_performed"              
# [65] "pre_bronchodilator_fev1_percent"                
# [66] "post_bronchodilator_fev1_percent"               
# [67] "pre_bronchodilator_fev1_fvc_percent"            
# [68] "post_bronchodilator_fev1_fvc_percent"           
# [69] "dlco_predictive_percent"                        
# [70] "radiation_therapy"                              
# [71] "postoperative_rx_tx"                            
# [72] "primary_therapy_outcome_success"                
# [73] "has_new_tumor_events_information"               
# [74] "has_drugs_information"                          
# [75] "has_radiations_information"                     
# [76] "has_follow_ups_information"                     
# [77] "stage_event_system_version"                     
# [78] "stage_event_clinical_stage"                     
# [79] "stage_event_pathologic_stage"                   
# [80] "stage_event_tnm_categories"                     
# [81] "stage_event_psa"                                
# [82] "stage_event_gleason_grading"                    
# [83] "stage_event_ann_arbor"                          
# [84] "stage_event_serum_markers"                      
# [85] "stage_event_igcccg_stage"                       
# [86] "stage_event_masaoka_stage"                      
# [87] "system_version"                                 
# [88] "clinical_stage"                                 
# [89] "pathologic_stage"                               
# [90] "tnm_categories"                                 
# [91] "psa"                                            
# [92] "gleason_grading"                                
# [93] "ann_arbor"                                      
# [94] "serum_markers"                                  
# [95] "igcccg_stage"                                   
# [96] "masaoka_stage"                                  
# [97] "days_to_new_tumor_event_after_initial_treatment"
# [98] "new_neoplasm_event_types"                       
# [99] "progression_determined_by_list"                 
# [100] "locoregional_procedure"                         
# [101] "metastatic_procedure"                           
# [102] "additional_radiation_therapy"                   
# [103] "additional_pharmaceutical_therapy"              
# [104] "project" 

saveRDS(m, file = 'tcga/luad/clinical.rds')
saveRDS(xml.drug, file = 'tcga/luad/drugs.rds')
saveRDS(xml.radiation, file = 'tcga/luad/radiations.rds')
saveRDS(xml.followup, file = 'tcga/luad/followups.rds')
