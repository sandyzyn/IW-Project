import pandas as pd
from scipy.stats import spearmanr, wilcoxon
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv('orco gene - Guidescan-Geneious Comparison.csv')

results = {
    'Exon': [],
    'Average_Guidescan2': [],
    'Average_Geneious': [],
    'Guidescan2_std': [],
    'Geneious_std': [],
    'Wilcoxon p-value (Paired)': [],
    'Spearman_Correlation': [],
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
    exon_data.columns = ['Guidescan2', 'Geneious']

    # only want numeric data
    exon_data['Guidescan2'] = pd.to_numeric(exon_data['Guidescan2'], errors='coerce')
    exon_data['Geneious'] = pd.to_numeric(exon_data['Geneious'], errors='coerce')
    exon_data = exon_data.dropna()

    if exon_data.empty:
        continue

    exon_data['Exon'] = f'Exon {i + 1}'
    boxplot_data.append(exon_data)

    avg_guidescan2 = exon_data['Guidescan2'].mean()
    avg_geneious = exon_data['Geneious'].mean()
    std_guidescan2 = exon_data['Guidescan2'].std()
    std_geneious = exon_data['Geneious'].std()

    paired_wilcoxon = wilcoxon(exon_data['Guidescan2'], exon_data['Geneious'])

    spearman_corr, spearman_pval = spearmanr(exon_data['Guidescan2'], exon_data['Geneious'])

    results['Exon'].append(f'Exon {i + 1}')
    results['Average_Guidescan2'].append(avg_guidescan2)
    results['Average_Geneious'].append(avg_geneious)
    results['Guidescan2_std'].append(std_guidescan2)
    results['Geneious_std'].append(std_geneious)
    results['Wilcoxon p-value (Paired)'].append(paired_wilcoxon.pvalue)
    results['Spearman_Correlation'].append(spearman_corr)
    results['Spearman_p-value'].append(spearman_pval)


boxplot = pd.concat(boxplot_data)

# Save to Excel
results = pd.DataFrame(results)
results.to_excel('comparisons_G_GS.xlsx', index=False)

# Boxplot visual 
plt.figure(figsize=(10, 6))
sns.boxplot(
    data=boxplot.melt(id_vars=['Exon'], value_vars=['Guidescan2', 'Geneious'], 
                         var_name='Software', value_name='Specificity Score'),
    x='Exon',
    y='Specificity Score',
    hue='Software',
)
plt.title('Comparison of Specificity Scores: Guidescan2 vs Geneious')
plt.xlabel('Exon')
plt.ylabel('Specificity Score')
plt.legend(title='Software')
plt.tight_layout()
plt.show()

# Correlation visual
correlation = pd.concat(
    [data.iloc[:, [1]].rename(columns={data.columns[1]: 'Guidescan2'}),
     data.iloc[:, [2]].rename(columns={data.columns[2]: 'Geneious'})],
    axis=1
).dropna()

correlation['Guidescan2'] = pd.to_numeric(correlation['Guidescan2'], errors='coerce')
correlation['Geneious'] = pd.to_numeric(correlation['Geneious'], errors='coerce')

correlation = correlation.dropna()

plt.figure(figsize=(12,6))

sns.scatterplot(
    data=boxplot,  
    x='Guidescan2',
    y='Geneious',
    hue='Exon',  
)

plt.title('Correlation Between Guidescan2 and Geneious Scores', fontsize=14)
plt.xlabel('Guidescan2 Specificity Score', fontsize=12)
plt.ylabel('Geneious Specificity Score', fontsize=12)
plt.legend(title='Exon', bbox_to_anchor=(1.05, 1), loc='upper left') 
plt.tight_layout()
plt.show()
