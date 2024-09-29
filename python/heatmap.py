import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

# Creating a DataFrame with the rnaseq coexpression data
file_name = 'perturb_rbp_coexp_test.csv'
rnaseq = pd.read_csv(file_name, index_col=None)
rnaseq.index = list(rnaseq.columns)
print(len(rnaseq.index))

# Create the heatmap
img = sns.clustermap(rnaseq, cmap="viridis", annot=False)#, annot=False, linewidths=0.5, linecolor='gray')
img.figure.savefig('562_test_heatmap.png')