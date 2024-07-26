import os
import hail as hl
hl.init(default_reference='GRCh38', idempotent=True)

weights_file_from_local = 'name_of_your_weights_file.csv'

#includes chr, risk, ref, rsid, weight, pos, a1, a1
weights_dataframe = pd.read_csv(weights_file_from_local)

#import hg19 weights from All of Us genomics bucket location

#Liftover preparation
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')  
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

#Liftover locus

weights_files = [weights_file_from_local]
path_to_file = 'your/path/here'

for file in weights_files:
    weights_path = f'{my_bucket}/path_to_file/{file}'

    weights = hl.import_table(weights_path, types={"rsid":"tstr",
                                                  "risk_allele":"tstr",
                                                  "weight":"tfloat",
                                                  "chr":"tstr",
                                                  "pos":"tint32",
                                                  "a1":"tstr",
                                                  "a2":"tstr"}, delimiter=',')
    
    weights = weights.annotate(locus_hg19 = hl.locus(weights.chr, weights.pos, reference_genome='GRCh37'))
    weights = weights.annotate(locus_hg38 = hl.liftover(weights.locus_hg19, 'GRCh38'))
    
    n_orig = weights.count()
    n_lifted = weights.aggregate(hl.agg.count_where(hl.is_defined(weights.locus_hg38)))
    
    weights = weights.annotate(chr=weights.locus_hg38.contig,
                          pos = weights.locus_hg38.position,
                          id=hl.variant_str(weights.locus_hg38, hl.sorted([weights.a1, weights.a2])))
    
    weights.filter(hl.is_defined(weights.locus_hg38)).drop('locus_hg19', 'locus_hg38').to_pandas().to_csv(f'{weights_path.replace("hg19", "hg38")}',
                                                               index=False)
    
    print(f'{file}: {n_lifted}/{n_orig} ({100*n_lifted/n_orig:.2f}%) sites successfully lifted over')
    print(f'lifted over weights file written to {weights_path.replace("hg19", "hg38")}')
















































