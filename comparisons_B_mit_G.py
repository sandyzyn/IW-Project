import pandas as pd
from scipy.stats import spearmanr, wilcoxon
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv('orco gene - MIT-Broad - Geneious Comparisons.csv')

results = {
    'Exon': [],
    'Average_Benchling': [],
    'Average_Geneious': [],
    'Benchling_std': [],
    'Geneious_std': [],
    'Wilcoxon p-value (Paired)': [],
    'Spearman_Corr': [],
    'Spearman_p-value': []
}

# 8 exons
num_exons = 8  

boxplot_data = []

# Iterate through each exon
for i in range(num_exons):
    start_col = i * 4  
    exon_data = data.iloc[:, [start_col + 1, start_col + 2]].dropna()

    # Rename 
    exon_data.columns = ['Benchling', 'Geneious']

    # only want numeric data
    exon_data['Benchling'] = pd.to_numeric(exon_data['Benchling'], errors='coerce')
    exon_data['Geneious'] = pd.to_numeric(exon_data['Geneious'], errors='coerce')
    exon_data = exon_data.dropna()

    if exon_data.empty:
        continue

    exon_data['Exon'] = f'Exon {i + 1}'
    boxplot_data.append(exon_data)

    avg_benchling = exon_data['Benchling'].mean()
    avg_geneious = exon_data['Geneious'].mean()
    std_benchling = exon_data['Benchling'].std()
    std_geneious = exon_data['Geneious'].std()

    paired_wilcoxon = wilcoxon(exon_data['Benchling'], exon_data['Geneious'])

    spearman_corr, spearman_pval = spearmanr(exon_data['Benchling'], exon_data['Geneious'])

    # Append results
    results['Exon'].append(f'Exon {i + 1}')
    results['Average_Benchling'].append(avg_benchling)
    results['Average_Geneious'].append(avg_geneious)
    results['Benchling_std'].append(std_benchling)
    results['Geneious_std'].append(std_geneious)
    results['Wilcoxon p-value (Paired)'].append(paired_wilcoxon.pvalue)
    results['Spearman_Corr'].append(spearman_corr)
    results['Spearman_p-value'].append(spearman_pval)

boxplot_df = pd.concat(boxplot_data)

# Save to Excel
results_df = pd.DataFrame(results)
results_df.to_excel('comparisons_B_mit_G.xlsx', index=False)

# Boxplot visual
plt.figure(figsize=(10, 6))
sns.boxplot(
    data=boxplot_df.melt(id_vars=['Exon'], value_vars=['Benchling', 'Geneious'], 
                         var_name='Software', value_name='Specificity Score'),
    x='Exon',
    y='Specificity Score',
    hue='Software',
)
plt.title('Comparison of Specificity Scores: Benchling (MIT-Broad) vs Geneious')
plt.xlabel('Exon')
plt.ylabel('Specificity Score')
plt.legend(title='Software')
plt.tight_layout()
plt.show()

# Correlation visual
regplot_data = pd.concat(
    [data.iloc[:, [1]].rename(columns={data.columns[1]: 'Benchling'}),
     data.iloc[:, [2]].rename(columns={data.columns[2]: 'Geneious'})],
    axis=1
).dropna()

regplot_data['Benchling'] = pd.to_numeric(regplot_data['Benchling'], errors='coerce')
regplot_data['Geneious'] = pd.to_numeric(regplot_data['Geneious'], errors='coerce')

regplot_data = regplot_data.dropna()

plt.figure(figsize=(12,6))

sns.scatterplot(
    data=boxplot_df,  
    x='Benchling',
    y='Geneious',
    hue='Exon',  
)

plt.title('Correlation Between Benchling(MIT-Broad) and Geneious Scores', fontsize=14)
plt.xlabel('Benchling Specificity Score', fontsize=12)
plt.ylabel('Geneious Specificity Score', fontsize=12)
plt.legend(title='Exon', bbox_to_anchor=(1.05, 1), loc='upper left')  
plt.tight_layout()
plt.show()
