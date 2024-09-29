import pandas as pd
import numpy as np

def convert_to_wide(matrix, bait, prey, values):
    wide = matrix.pivot_table(index=prey, columns=bait, values=values)
    wide = wide.fillna(0)
    return wide

def standardized_rank(col, pval=False):
    #if pval, then negative log transform first. add a pseudocount to 0.
    if pval:
        nonzero = col[col > 0]
        if len(nonzero) == 0:
            pseudocount = 1e-10
        else:
            pseudocount = min(nonzero)/2
        col += pseudocount
        # neg log transform col
        col = np.log(col) * -1
        
    rank = col.rank(ascending=True, method='average')
    # standardize by dividing all ranks by max_rank
    std_rank_score = 1 - rank / max(rank)
    return std_rank_score

def to_output_matrix_pair(std_rank_score, top_k=100):
    prot_mat = []
    score_mat = []
    for x in std_rank_score.columns:
        prot_mat.append(_top_proteins(std_rank_score[x], top_k))
        score_mat.append(_top_scores(std_rank_score[x], top_k))
    prot_mat = pd.concat(prot_mat, axis=1)
    prot_mat.columns = std_rank_score.columns
    score_mat = pd.concat(score_mat, axis=1)
    return prot_mat, score_mat
    
def _top_proteins(col, top_k):
    col = col.sort_values(ascending=False)
    hits = col[0:top_k].index
    hits = hits.to_series()
    hits.index = range(1,top_k+1)
    return hits

def _top_scores(col, top_k):
    col = col.sort_values(ascending=False)
    scores = col[0:top_k]
    scores.index = range(1, top_k+1)
    return scores

def write_to_csv(matrix, filename='my_matrix.csv'):
    matrix.to_csv(f'../output/{filename}')
    
def full_workflow(in_file, bait_col, prey_col, values_col, top_k, out_file_prefix, pval=False):
    # read 
    print(in_file)
    df = pd.read_csv(in_file, sep='\t')
    df_wide = convert_to_wide(df, bait=bait_col, prey=prey_col, values=values_col)
    df_wide = df_wide.fillna(0)
    
    std_rank_score = df_wide.apply(lambda x: standardized_rank(x, pval=pval), axis=1)
    
    prot, score = to_output_matrix_pair(std_rank_score, top_k)
    
    write_to_csv(df_wide, filename=f'../output/{out_file_prefix}_ppi_matrix.csv')
    write_to_csv(prot, filename=f'../output/{out_file_prefix}_top_prot_per_bait.csv')
    write_to_csv(score, filename=f'../output/{out_file_prefix}_top_score_per_bait.csv')
    print('done')

if __name__ == '__main__':
    
    full_workflow('../output/goos2022.preprocessed_PPIs.HEK293.standardized.tsv', 'gene1', 'gene2', 'conf_score', 100, 'goos2022-HEK293')
    full_workflow('../output/huttlin2021.preprocessed_PPIs.HEK293T.standardized.tsv', 'gene1', 'gene2', 'pscore', 100, 'huttlin2021-HEK293T')
    full_workflow('../output/khoroshkin2024.preprocessed_PPIs.K562.standardized.tsv', 'gene1', 'gene2', 'conf_score', 100, 'khoroshkin2024-K562')
    full_workflow('../output/johnson2021.preprocessed_PPIs.HEK293T.standardized.tsv', 'gene1', 'gene2', 'pvalue', 100, 'johnson2021-HEK293T', pval=True)