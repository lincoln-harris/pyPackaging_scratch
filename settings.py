import pandas as pd 

def init():
	global mutationsDF
	global summaryTable
	global patientMetadata

	mutationsDF = pd.DataFrame(columns=['cell', 'brafMut', 'egfrMut', 'krasMut'])

	cols = ['cell', 'patient', 'clinical_driver_gene', 'clinical_mutation', 'coverage_to_ROI', 'clin_mut_found_bool', 'mutations_found_EGFR', 'mutations_found_BRAF', 'mutations_found_KRAS', 'fusions_found', 'tumorCell_bool']
	summaryTable = pd.DataFrame(columns=cols)

	patientMetadata = pd.read_csv('/Users/lincoln.harris/code/SNP_calling_pipeline/cDNA_plate_metadata.csv')