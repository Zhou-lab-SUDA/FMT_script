import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from matplotlib.backends.backend_pdf import PdfPages

def calculate_roc_auc(df, cohort_name=None):
    """Calculate ROC curve and AUC for a given cohort"""
    if cohort_name:
        cohort_df = df[df['cohort'] == cohort_name]
    else:
        cohort_df = df  # All cohorts
    
    fpr, tpr, thresholds = roc_curve(cohort_df['status_binary'], 
                                     cohort_df['proportion_healthy'])
    roc_auc = auc(fpr, tpr)
    return fpr, tpr, roc_auc, len(cohort_df)

# Load data
df = pd.read_csv('input_data.tsv', sep=',')  # Update with your filename

# Create PDF for output
with PdfPages('roc_curves.pdf') as pdf:
    plt.figure(figsize=(10, 8))
    
    # Plot ROC for all cohorts combined
    fpr_all, tpr_all, auc_all, n_all = calculate_roc_auc(df)
    plt.plot(fpr_all, tpr_all, color='navy', lw=2, 
             label=f'All Cohorts (AUC = {auc_all:.2f}, N={n_all})')
    
    # Define colors and markers for each cohort
    cohort_styles = {
        'rCDI': ('firebrick', 'o'),
        'IBS': ('darkorange', 's'),
        'LUAD': ('forestgreen', '^'),
        'MEL': ('darkviolet', 'd')
    }
    
    # Plot ROC for each cohort
    for cohort, (color, marker) in cohort_styles.items():
        fpr_cohort, tpr_cohort, auc_cohort, n_cohort = calculate_roc_auc(df, cohort)
        plt.plot(fpr_cohort, tpr_cohort, color=color, lw=1.5, 
                 marker=marker, markevery=0.1, markersize=8,
                 label=f'{cohort} (AUC = {auc_cohort:.2f}, N={n_cohort})')
    
    # Format plot
    plt.plot([0, 1], [0, 1], 'k--', lw=1)  # Diagonal line
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title('ROC Curves by Cohort', fontsize=14)
    plt.legend(loc="lower right", fontsize=10)
    plt.grid(True, alpha=0.3)
    
    # Add AUC table
    table_data = []
    for cohort in cohort_styles.keys():
        _, _, auc_val, n_val = calculate_roc_auc(df, cohort)
        table_data.append([cohort, f'{auc_val:.3f}', n_val])
    
    # Add overall AUC
    table_data.append(['All Cohorts', f'{auc_all:.3f}', n_all])
    
    # Create table
    col_labels = ['Cohort', 'AUC', 'Samples']
    plt.table(cellText=table_data,
              colLabels=col_labels,
              cellLoc='center',
              loc='lower left',
              bbox=[0.15, 0.55, 0.25, 0.35])
    
    # Save to PDF
    pdf.savefig(bbox_inches='tight')
    plt.close()
    
    # Print AUC values to console
    print("AUC Values:")
    print("-----------")
    for row in table_data:
        print(f"{row[0]}: {row[1]} (N={row[2]})")