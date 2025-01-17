import pandas as pd
from scipy.stats import spearmanr, wilcoxon
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv('orco gene - CFD - Guidescan Comparison .csv')

results = {
    'Exon': [],
    'Average_Benchling': [],
    'Average_Guidescan2': [],
    'Benchling_std': [],
    'Guidescan2_std': [],
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
    exon_data.columns = ['Benchling', 'Guidescan2']

    # only want numeric data
    exon_data['Benchling'] = pd.to_numeric(exon_data['Benchling'], errors='coerce')
    exon_data['Guidescan2'] = pd.to_numeric(exon_data['Guidescan2'], errors='coerce')
    exon_data = exon_data.dropna()

    if exon_data.empty:
        continue

    exon_data['Exon'] = f'Exon {i + 1}'
    boxplot_data.append(exon_data)

    avg_benchling = exon_data['Benchling'].mean()
    avg_guidescan2 = exon_data['Guidescan2'].mean()
    std_benchling = exon_data['Benchling'].std()
    std_guidescan2 = exon_data['Guidescan2'].std()

    paired_wilcoxon = wilcoxon(exon_data['Benchling'], exon_data['Guidescan2'])

    spearman_corr, spearman_pval = spearmanr(exon_data['Benchling'], exon_data['Guidescan2'])

    results['Exon'].append(f'Exon {i + 1}')
    results['Average_Benchling'].append(avg_benchling)
    results['Average_Guidescan2'].append(avg_guidescan2)
    results['Benchling_std'].append(std_benchling)
    results['Guidescan2_std'].append(std_guidescan2)
    results['Wilcoxon p-value (Paired)'].append(paired_wilcoxon.pvalue)
    results['Spearman_Corr'].append(spearman_corr)
    results['Spearman_p-value'].append(spearman_pval)

boxplot = pd.concat(boxplot_data)

# Save to Excel
results = pd.DataFrame(results)
results.to_excel('comparisons_B_cfd_GS.xlsx', index=False)

# Boxplot visual
plt.figure(figsize=(10, 6))
sns.boxplot(
    data=boxplot.melt(id_vars=['Exon'], value_vars=['Benchling', 'Guidescan2'], 
                         var_name='Software', value_name='Specificity Score'),
    x='Exon',
    y='Specificity Score',
    hue='Software',
)
plt.title('Comparison of Specificity Scores: Benchling (CFD) vs Guidescan2')
plt.xlabel('Exon')
plt.ylabel('Specificity Score')
plt.legend(title='Software')
plt.tight_layout()
plt.show()

# Correlation visual
correlation = pd.concat(
    [data.iloc[:, [1]].rename(columns={data.columns[1]: 'Benchling'}),
     data.iloc[:, [2]].rename(columns={data.columns[2]: 'Guidescan2'})],
    axis=1
).dropna()

correlation['Benchling'] = pd.to_numeric(correlation['Benchling'], errors='coerce')
correlation['Guidescan2'] = pd.to_numeric(correlation['Guidescan2'], errors='coerce')

correlation = correlation.dropna()

plt.figure(figsize=(12,6))

sns.scatterplot(
    data=boxplot,  
    x='Benchling',
    y='Guidescan2',
    hue='Exon',  
)

plt.title('Correlation Between Benchling(CFD) and Guidescan2 Scores', fontsize=14)
plt.xlabel('Benchling Specificity Score', fontsize=12)
plt.ylabel('Guidescan2 Specificity Score', fontsize=12)
plt.legend(title='Exon', bbox_to_anchor=(1.05, 1), loc='upper left') 
plt.tight_layout()
plt.show()
