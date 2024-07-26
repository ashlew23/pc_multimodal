import os
import hail as hl
hl.init(default_reference='GRCh38', idempotent=True)

from bokeh.io import show, output_notebook
from bokeh.layouts import gridplot
import pandas as pd
output_notebook()

#Load test/train mt
bucket = os.getenv("WORKSPACE_BUCKET")
path = 'path/to/mt'
mt_path = f'{bucket}/path/test_and_train_file.mt'

mt_cohort = hl.read_matrix_table(mt_path)

#Score

name_of_file_in_bucket = 'weights_file.csv'
weights_file = pd.read_csv(name_of_file_in_bucket)
weights_file.head()

weights_files_base = [name_of_file_in_bucket]

weights_hg38= f'{bucket}/weights/hg38/data'

def load_weights(path):
    weights_ht = hl.import_table(path, types={"id":"tstr",
                                                  "effect_allele":"tstr",
                                                  "weight":"tfloat",
                                                  "chr":"tstr",
                                                  "pos":"tint32",
                                                  "a1":"tstr",
                                                  "a2":"tstr"},
                             delimiter=",")
    condition = os.path.basename(path).replace("_PRS_weights.txt","")
    weights_ht = weights_ht.annotate_globals(condition=condition)
    weights_ht = weights_ht.annotate(locus=hl.locus(weights_ht.chr, weights_ht.pos))
    weights_ht = weights_ht.key_by(weights_ht.locus)
    return weights_ht

for weight_file in weights_files_base:
    weights = load_weights(os.path.join(weights_hg38, weight_file))
    
    mt_cohort_scored = mt_cohort.annotate_rows(effect_allele = weights[mt_cohort.locus].risk, 
                            weight=weights[mt_cohort.locus].weight)
    mt_cohort_scored = mt_cohort_scored.filter_rows(hl.is_defined(mt_cohort_scored.effect_allele))
    mt_cohort_scored = mt_cohort_scored.annotate_rows(effect_allele_index = 
                                                      mt_cohort_scored.alleles.index(
                                                          mt_cohort_scored.effect_allele)
                                                     )
    mt_cohort_scored = mt_cohort_scored.annotate_cols(score=
                                              hl.agg.sum(
                                                  hl.if_else(
                                                      hl.is_defined(mt_cohort_scored.effect_allele_index),
                                                      mt_cohort_scored.GT.
                                                      one_hot_alleles(mt_cohort_scored.alleles.length())[mt_cohort_scored.effect_allele_index],
                                                      0)
                                                  * mt_cohort_scored.weight
                                              )
                                             )
    
    score_pd = mt_cohort_scored.cols().to_pandas()
    #out_path = os.path.join(bucket, "cohorts", f'{weight_file.replace("_PRS_weights.txt","")}_train_and_test_scores_new.csv')
    sample_path = os.path.join(bucket, "cohorts", 'train_and_test_scores_jul25.csv')
    score_pd.to_csv(sample_path, index = False)