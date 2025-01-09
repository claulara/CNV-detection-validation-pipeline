import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec

def load_data(filepath, fillna_value=0):
    data = pd.read_csv(filepath)
    data.fillna(fillna_value, inplace=True)
    data['svtype'] = data['svtype'].astype('category')
    return data


def plot_boxplots_by_algorithm_by_svtype(data, metrics, algorithm_order):
    """
    Create boxplots for metrics grouped by algorithm and SV type.
    """
    sns.set_style("whitegrid")
    custom_palette = {
        'Deletion': '#FF69B4',  
        'Duplication': '#D2691E' }
    # Create subplots
    data['Algoritmo'] = pd.Categorical(data['Algoritmo'], categories=algorithm_order, ordered=True)
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 5))
    for ax, metric in zip(axes, metrics):
        sns.boxplot(
            x="Algoritmo",
            y=metric,
            hue="svtype",
            data=data,
            palette=custom_palette,
            ax=ax)
        ax.set_title(f"{metric}")
        ax.set_ylabel(metric)
        ax.tick_params(axis='x', rotation=45)
        if ax.get_legend():
            ax.legend(title="Tipo de variante", loc='upper left')
    for ax in axes[1:]:
        ax.legend().remove()
    plt.tight_layout()
    plt.show()



def plot_lines_by_algorithm_by_svtype(data):
    """
    Generate line plots for metrics grouped by algorithm and SV type.
    """
    sns.set_style("whitegrid")
    # Transform data to long format for lineplot
    df = pd.melt(
        data,
        id_vars=['Algoritmo', 'svtype'],
        value_vars=['Precision', 'Recall', 'Fscore'],
        var_name='Métrica',
        value_name='Valor')
    df['Algoritmo'] = pd.Categorical(df['Algoritmo'], ordered=True)
    # Create plot
    plt.figure(figsize=(12, 7))
    sns.lineplot(
        data=df,
        x='Algoritmo',
        y='Valor',
        hue='svtype',
        style='Métrica',
        markers=True,
        dashes={'Precision': '', 'Recall': (2, 2), 'Fscore': (4, 2)},
        palette="bright")
    plt.xlabel('Algoritmo', fontsize=11)
    plt.ylabel('Valor', fontsize=11)
    plt.xticks(fontsize=12, rotation=0)
    plt.legend(fontsize=11, loc='upper left')
    plt.tight_layout()
    plt.show()




def plot_total_variants_by_algorithm_and_replica(data):
    """
    Generate bar plots for variant number mean grouped by sample and algorithm, spliting by CNV type.
    """
    sns.set_style("whitegrid")
    custom_palette = {
        'Deletion': '#FF69B4',  
        'Duplication': '#D2691E' 
    }
    data_filtered = data[data['gold_standard'] == 'general']
    
    # QueryTotalCount is the variant number 
    data_filtered = data_filtered.dropna(subset=['QueryTotalCount'])
    data_filtered['QueryTotalCount'] = pd.to_numeric(data_filtered['QueryTotalCount'], errors='coerce')
    data_filtered = data_filtered.dropna(subset=['QueryTotalCount'])
    data_grouped = data_filtered.groupby(['Algoritmo', 'Muestra', 'svtype'], as_index=False)['QueryTotalCount'].mean()
    counts = data_filtered.groupby(['Algoritmo', 'Muestra', 'svtype']).size().reset_index(name='count')
    # three pipelines for each sample
    expected_count = 3
    valid_groups = counts[counts['count'] == expected_count][['Algoritmo', 'Muestra', 'svtype']]
    data_grouped = pd.merge(data_grouped, valid_groups, on=['Algoritmo', 'Muestra', 'svtype'])
    algorithms = data_grouped['Algoritmo'].unique()
    fig_width_per_algo = 8 
    fig_height = 6
    fig = plt.figure(figsize=(fig_width_per_algo * len(algorithms), fig_height))
    
    # GridSpec for the subplots 
    gs = gridspec.GridSpec(1, len(algorithms), figure=fig, wspace=0.2)
    axes = []
    for i, algorithm in enumerate(algorithms):
        ax = fig.add_subplot(gs[0, i])
        axes.append(ax)
        algo_data = data_grouped[data_grouped['Algoritmo'] == algorithm]
        sns.barplot(
            data=algo_data,
            x='Muestra',
            y='QueryTotalCount',
            hue='svtype',
            palette=custom_palette,
            ax=ax,
            width=0.6)
        ax.set_title(f"{algorithm}", fontsize=12)
        ax.tick_params(axis='x', rotation=70)
        ax.tick_params(axis='y', rotation=90) 
        if i == 0:
            ax.set_ylabel("Nº promedio de variantes reportadas en las tres pipelines", fontsize=11)
        else:
            ax.set_ylabel("")
        ax.legend_.remove()  
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, fontsize=12)
    plt.tight_layout(rect=[0.03, 0.03, 0.97, 0.90])
    plt.show()



def main(filepath):
    """
    Main function to load data and generate all plots
    """
    # load data
    data = load_data(filepath)
    # define orders and metrics
    kit_order = ["IDT-V1", "IDT-V2", "ROCHE-V1"]
    replica_order = ["Réplica 1", "Réplica 2", "Réplica 3", "Réplica 4", "Réplica 5"]
    preprocessing_order = ["bowtie2-picard-gatk3", "bowtie2-samtools-gatk4", "minimap2-samtools-gatk4"]
    algorithm_order= ["cnmops", "cnvkit", "contra", "exomedepth", "manta", "laconv", "xhmm"]
    metrics = ["Precision", "Recall", "Fscore"]
    # generate plots
    plot_boxplots_by_algorithm_by_svtype(data, metrics, algorithm_order)
    plot_lines_by_algorithm_by_svtype(data)
    plot_total_variants_by_algorithm_and_replica(data)

# execute the script
if __name__ == "__main__":
    filepath = "/Users/claulara/Desktop/TFM/results_posprocess/metrics/svtype_metrics.csv"
    main(filepath)
