import seaborn as sns
import pandas as pd

# Creating a DataFrame with the rnaseq coexpression data
file_name = 'small_huvec_corr_mat.csv'
rnaseq = pd.read_csv(file_name, index_col=None)
rnaseq.index = list(rnaseq.columns)


# Create the heatmap
img = sns.clustermap(rnaseq, cmap="viridis", annot=False)
img.figure.savefig('562_test_heatmap.png')