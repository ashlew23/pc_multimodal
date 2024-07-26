import os
import hail as hl
hl.init(default_reference='GRCh38', idempotent=True)

from bokeh.io import show, output_notebook
from bokeh.layouts import gridplot
output_notebook()

# Path to the Matrix Table in GCS
mt_path = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/multiMT/hail.mt'

# Read the Matrix Table
mt_wgs = hl.read_matrix_table(mt_path)

mt_wgs.count()

#Merge Phenotypes
phenotype_filename = f'{bucket}/disease_phenotype.tsv'
phenotype_filename
df = phenotype_filename

# Create binary columns for each category in 'race'
race_dummies = pd.get_dummies(df['race'])

# Concatenate these new columns with the original DataFrame
df = pd.concat([df, race_dummies], axis=1)

# Drop the original 'race' column if it's no longer needed
df = df.drop(columns=['race'])

#Convert True/False to 1/0
df = df.apply(lambda x: x.astype(int) if x.dtype == 'bool' else x)

df['person_id'] = df['person_id'].astype(str)  # Ensure person_id is a string

# Convert DataFrame to Hail Table
hl_df = hl.Table.from_pandas(df)

# Define column types
column_types = {
    'person_id': hl.tstr,
    'age': hl.tint32,
    'PC': hl.tint32,
    'is_male': hl.tint32,
    'race_Asian': hl.tint32,
    'race_Black': hl.tint32,
    'race_MiddleEastern_NorthAfrican': hl.tint32,
    'race_Multiracial': hl.tint32,
    'race_Native_Hawaiian_Pacific_Islander': hl.tint32,
    'race_Other': hl.tint32,
    'race_Skip': hl.tint32,
    'race_White': hl.tint32
}

# Annotate Hail Table columns to the specified types, ensuring person_id is a string
hl_df = hl_df.annotate(
    **{
        'person_id': hl_df['person_id'],  # person_id already string in DataFrame
        **{col: hl.int32(hl_df[col]) for col in column_types if col != 'person_id' and col in hl_df.row}
    }
)

# Key the Hail Table by 'person_id'
phenotypes = hl_df.key_by('person_id')

# Save the Table if needed
phenotypes.write('prostate_cancer_phenotype.mt', overwrite=True)

# Verify the Table
phenotypes.describe()

mt = mt_wgs.semi_join_cols(phenotypes)
mt.count()

mt = mt.annotate_cols(pheno = phenotypes[mt.s])

#View mt
mt.entries().show()

#Remove related samples
# Define the path to the related samples file
related_samples_path = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv'

# Import the table, ensuring the column types are correctly specified
related_remove = hl.import_table(
    related_samples_path,
    types={"sample_id": hl.tstr},  # Use the correct column name
    key="sample_id"  # Set the correct column as the key
)

# Verify the table structure to ensure everything is set correctly
related_remove.describe()

mt = mt.anti_join_cols(related_remove)

#Remove Flagged Samples
flagged_samples_path = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv'
flagged_samples = hl.import_table(flagged_samples_path, key='s')

mt = mt.anti_join_cols(flagged_samples)
mt.count()

#Annotate with ancestry
ancestry_pred_path = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'
ancestry_pred = hl.import_table(ancestry_pred_path, key='research_id')

mt = mt.annotate_cols(ancestry_pred = ancestry_pred[mt.s].ancestry_pred)

samples_pd = mt.cols().to_pandas()

# see counts by ancestry
samples_pd.groupby("ancestry_pred")["ancestry_pred"].count()

selected_samples_pd = samples_pd

# check that ancestry counts look good
selected_samples_pd.groupby("ancestry_pred")["ancestry_pred"].count()

#Split into training and testing
# split the downsampled_cohort randomly in half (for each ancestry) between training and testing
training_cohort_pd = selected_samples_pd.groupby("ancestry_pred").sample(frac=0.8)
test_cohort_pd = selected_samples_pd.drop(training_cohort_pd.index)

# check counts
print(f'training: {training_cohort_pd.groupby("ancestry_pred")["ancestry_pred"].count()}')
print(f'test: {test_cohort_pd.groupby("ancestry_pred")["ancestry_pred"].count()}')

training_cohort = hl.Table.from_pandas(training_cohort_pd.reset_index(), key='s')
test_cohort = hl.Table.from_pandas(test_cohort_pd.reset_index(), key='s')
full_cohort = hl.Table.from_pandas(selected_samples_pd.reset_index(), key='s')

#Extract PCs
ancestry_pred_training = ancestry_pred.semi_join(training_cohort)
ancestry_pred_training_pd = ancestry_pred_training.to_pandas()

ancestry_pred_test = ancestry_pred.semi_join(test_cohort)
ancestry_pred_test_pd = ancestry_pred_test.to_pandas()

ancestry_pred_training_pd[[f'pc_{i}' for i in range(1,17)]] = pd.DataFrame([ast.literal_eval(x) for x in ancestry_pred_training_pd.pca_features.tolist()],
                                                                    index = ancestry_pred_training_pd.index)

ancestry_pred_test_pd[[f'pc_{i}' for i in range(1,17)]] = pd.DataFrame([ast.literal_eval(x) for x in ancestry_pred_test_pd.pca_features.tolist()],
                                                                    index = ancestry_pred_test_pd.index)

#Subset to locations in weights files and write to bucket
mt = mt.filter_rows(hl.is_defined(weights_sites[mt.locus]))
mt = mt.select_entries(mt.GT)
mt = mt.select_rows()
# Repartitioning helps distribute the data more evenly across the workers, which can improve parallel processing efficiency
n_partitions = 100  # Example value; adjust based on your needs
mt = mt.repartition(n_partitions)

mt.write(f'{bucket}/cohorts/test_and_train_file.mt', overwrite=True)



