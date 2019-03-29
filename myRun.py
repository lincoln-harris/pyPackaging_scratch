#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////
# script: myRun.py
# author: Lincoln
# date: 3.28.19
#
# Takes as input GOI_out_AA.csv files (from getMutationCounts_overall_and_GOI.py),
# patient metadata, seurat metadata, fusionsDF, and creates a BY CELL 
# summaryTable. The goal with this table is to provide an answer to questions like
# 'which patients have which mutations?', and 'how many cells have clinically relevant
# mutations?'. Currently i've got an ipynb that accomplishes this, but its like
# super long and unwieldy, so i though converting it to a script would be a good
# idea. Lets try and make this more modular and flowy. 
#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////
import myLib
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # want to disable this SettingWithCopyWarning

print('running...')

# READ IN ALL OF THESE BY-GENE AMINO-ACID LEVEL MUTATION COUNTS OBJECTS
mutsPATH = '/Users/lincoln.harris/code/SNP_calling_pipeline/getMutationCounts/'
egfrPATH = mutsPATH + 'egfr_germline_out_AA.csv'
brafPATH = mutsPATH + 'braf_germline_out_AA.csv'
krasPATH = mutsPATH + 'kras_germline_out_AA.csv'

egfr_df = pd.read_csv(egfrPATH, header=None, names=['cell', 'mutations'])
braf_df = pd.read_csv(brafPATH, header=None, names=['cell', 'mutations'])
kras_df = pd.read_csv(krasPATH, header=None, names=['cell', 'mutations'])

# FIRST STEP IS TO GENERATE THE mutationsDF
print('creating mutationsDF')
mutationsDF['cell'] = egfr_df['cell']
mutationsDF['egfrMut'] = egfr_df['mutations'] # fill in EGFR first -- this is ok bc
                                              #  the cell order is based on egfr_df

myLib.mutationsDF_fillIn('braf', braf_df) 
myLib.mutationsDF_fillIn('kras', kras_df)

myLib.removeExtraCharacters_mutationsDF('egfr')
myLib.removeExtraCharacters_mutationsDF('braf')
myLib.removeExtraCharacters_mutationsDF('kras')

# READ IN patientMetadata
patientMetadata = patientMetadata.drop([0,1]) # first two rows are wierd

# INIT THE SUMMARY TABLE
summaryTable['cell'] = mutationsDF['cell']

# FILL IN VARIOUS METADATA COLS
myLib.genericSummaryTableFillIn('patient_id', 'patient')
myLib.genericSummaryTableFillIn('driver_gene', 'clinical_driver_gene')
myLib.genericSummaryTableFillIn('driver_mutation', 'clinical_mutation')

# FILL IN MUTATIONS FOUND COL 
summaryTable['mutations_found_EGFR'] = mutationsDF['egfrMut']
summaryTable['mutations_found_KRAS'] = mutationsDF['krasMut']
summaryTable['mutations_found_BRAF'] = mutationsDF['brafMut']

# READ IN FUSIONS DATAFRAME, THEN FILL IN summaryTable
fusionsDF = pd.read_csv('./fusion_dataframe.csv')
myLib.fusionsFillIn(fusionsDF)

# SET UP A COL TO TRANSLATE 'RAW' MUTATION CALLS TO 'CLINICAL'
print('translating mutations')
summaryTable['mutations_found_translated'] = ""
myLib.translatedMutsFillIn_EGFR()
myLib.translatedMutsFillIn_nonEGFR('KRAS')
myLib.translatedMutsFillIn_nonEGFR('BRAF')
myLib.translatedMutsFillIn_fusions()

# CONVERT LISTS TO STRING, SO I CAN GET SET -- probably not necessary, actually 
myLib.convertToString()

# FILL IN clin_mut_found_bool
myLib.clinMutFound_fillIn()

# FILL IN  tumorCellBool
myLib.tumorCellBoolFillIn()

# GET PER-CELL ROI COVERAGE DFs
print('getting coverage to ROIs...')
braf_V600E_cov_nonZero = myLib.getNonZeroCovROI('braf', 'V600E')
egfr_L858R_cov_nonZero = myLib.getNonZeroCovROI('egfr', 'L858R')
egfr_exon19del_cov_nonZero = myLib.getNonZeroCovROI('egfr', 'exon19del')
egfr_exon20ins_cov_nonZero = myLib.getNonZeroCovROI('egfr', 'exon20ins') # this guy is totally empty...
egfr_G719X_cov_nonZero = myLib.getNonZeroCovROI('egfr', 'G719X')
egfr_L861Q_cov_nonZero = myLib.getNonZeroCovROI('egfr', 'L861Q')
egfr_S768I_cov_nonZero = myLib.getNonZeroCovROI('egfr', 'S768I')
egfr_T790M_cov_nonZero = myLib.getNonZeroCovROI('egfr', 'T790M')
kras_G12C_cov_nonZero = myLib.getNonZeroCovROI('kras', 'G12C')
kras_G13X_cov_nonZero = myLib.getNonZeroCovROI('kras', 'G13X')
kras_Q61X_cov_nonZero = myLib.getNonZeroCovROI('kras', 'Q61X')

# FIX UP SOME OF THE WEIRD ONES
kras_G13X_cov_nonZero['depth_gvcf'][4202] = 34	
kras_Q61X_cov_nonZero['depth_gvcf'][6431] = 92
egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip('[')
egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip(']')
egfr_exon19del_cov_nonZero['depth_gvcf'] = egfr_exon19del_cov_nonZero['depth_gvcf'].str.strip("'")

# FILL IN ROI COVERAGE TO SUMMARY TABLE
myLib.ROI_coverage_fillIn(braf_V600E_cov_nonZero, 'BRAF', 'V600E')
myLib.ROI_coverage_fillIn(egfr_G719X_cov_nonZero, 'EGFR', 'G719X')
myLib.ROI_coverage_fillIn(egfr_L858R_cov_nonZero, 'EGFR', 'L858R')
myLib.ROI_coverage_fillIn(egfr_L861Q_cov_nonZero, 'EGFR', 'L861Q')
myLib.ROI_coverage_fillIn(egfr_S768I_cov_nonZero, 'EGFR', 'S768I')
myLib.ROI_coverage_fillIn(egfr_T790M_cov_nonZero, 'EGFR', 'T790M')
myLib.ROI_coverage_fillIn(kras_G12C_cov_nonZero, 'KRAS', 'G12C')
myLib.ROI_coverage_fillIn(kras_G13X_cov_nonZero, 'KRAS', 'G13X')
myLib.ROI_coverage_fillIn(kras_Q61X_cov_nonZero, 'KRAS', 'Q61X')
myLib.ROI_coverage_fillIn(egfr_exon19del_cov_nonZero, 'EGFR', 'del19')
myLib.ROI_coverage_fillIn(egfr_exon20ins_cov_nonZero, 'EGFR', 'ins20')

# TRIM IT DOWN
print('trimming down...')
summaryTable_trimmed = summaryTable[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 'tumorCell_bool', 'mutations_found_translated']]
summaryTable_trimmed.columns = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool', 'mutations_found']
summaryTable_trimmed = summaryTable_trimmed[['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'mutations_found', 'coverage_to_ROI', 'clinical_mutation_found_bool', 'tumorCell_bool']]

# write this fucker
print('writing...')
summaryTable_trimmed.to_csv('/Users/lincoln.harris/Desktop/validationTable_TEST.csv', index=False)

print('done!')

#/////////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////////