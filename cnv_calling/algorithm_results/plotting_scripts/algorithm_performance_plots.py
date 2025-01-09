import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec

def load_data(filepath, fillna_value=0):
    data = pd.read_csv(filepath)
    data.fillna(0, inplace=True)
    return data


def plot_boxplots_by_algorithm_by_goldstandard(data, algorithm_order, metrics):
    """
    Generate boxplot plots for metrics, grouped by algorithm and gold standard.
    """
    sns.set_style("whitegrid")
    # algorithm's order
    data['Algoritmo'] = pd.Categorical(data['Algoritmo'], categories=algorithm_order, ordered=True)
    # Create subplots
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 5))
    for ax, metric in zip(axes, metrics):
        sns.boxplot(
            x="Algoritmo",
            y=metric,
            hue="gold_standard",
            data=data,
            palette="bright",
            ax=ax)
        ax.set_title(f"{metric}")
        ax.set_ylabel(metric)
        ax.tick_params(axis='x', rotation=45)
        if ax.get_legend():
            ax.legend(title="Gold Standard", loc='upper left')
    for ax in axes[1:]:
        ax.legend().remove()
    plt.tight_layout()
    plt.show()


def plot_lines_by_algorithm_by_goldstandard(data, algorithm_order):
    """
    Generate line plots for metrics, grouped by capture kit and gold standard.
    """
    sns.set_style("whitegrid")
    custom_palette = {'general': 'blue', 'pass': 'orange'}
    # transform data to a long format for lineplot and sort by kit
    df = pd.melt(
        data,
        id_vars=['Algoritmo', 'gold_standard'],
        value_vars=['Precision', 'Recall'],
        var_name='Métrica',
        value_name='Valor')
    df['Algoritmo'] = pd.Categorical(df['Algoritmo'], categories=algorithm_order, ordered=True)
    # create plot
    plt.figure(figsize=(12, 6))
    sns.lineplot(
        data=df,
        x='Algoritmo',
        y='Valor',
        hue='gold_standard',
        style='Métrica',
        markers=True,
        dashes={'Precision': '', 'Recall': (2, 2)},
        palette=custom_palette)
    plt.xlabel('Algoritmo', fontsize=11)
    plt.ylabel('Valor', fontsize=11)
    plt.xticks(fontsize=12, rotation=0)
    plt.legend(fontsize=11, loc='upper left')
    plt.tight_layout()
    plt.show()


def plot_lines_by_algorithm_by_pipeline(data, pipeline_order):
    """
    Generate side-by-side line plots for Precision and Recall by replica.
    """
    sns.set_style("whitegrid")
    
    # Ensure 'Muestra' is categorical and ordered
    data['Preprocesado'] = pd.Categorical(data['Preprocesado'], categories=pipeline_order, ordered=True)
    
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Precision lineplot
    sns.lineplot(
        data=data,
        x="Preprocesado",
        y="Precision",
        hue="Algoritmo",
        marker="o",
        ax=axes[0]
    )
    axes[0].set_title("Precision por Pipeline Bioinformática")
    axes[0].set_xlabel("Pipeline Bioinformática")
    axes[0].set_ylabel("Precision")
    axes[0].tick_params(axis='x', rotation=0)
    axes[0].legend(loc='upper right')

    # Recall lineplot
    sns.lineplot(
        data=data,
        x="Preprocesado",
        y="Recall",
        hue="Algoritmo",
        marker="o",
        ax=axes[1]
    )
    axes[1].set_title("Recall por Pipeline Bioinformática")
    axes[1].set_xlabel("Pipeline Bioinformática")
    axes[1].set_ylabel("Recall")
    axes[1].tick_params(axis='x', rotation=0)
    axes[1].get_legend().remove() 
    
    plt.tight_layout()
    plt.show()


def plot_lines_by_algorithm_by_replica(data, replica_order):
    """
    Generate side-by-side line plots for Precision and Recall by replica.
    """
    sns.set_style("whitegrid")
    
    # Ensure 'Muestra' is categorical and ordered
    data['Muestra'] = pd.Categorical(data['Muestra'], categories=replica_order, ordered=True)
    
    # Create subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Precision lineplot
    sns.lineplot(
        data=data,
        x="Muestra",
        y="Precision",
        hue="Algoritmo",
        marker="o",
        ax=axes[0]
    )
    axes[0].set_title("Precision por Réplica")
    axes[0].set_xlabel("Réplica")
    axes[0].set_ylabel("Precision")
    axes[0].tick_params(axis='x', rotation=0)
    axes[0].legend(loc='upper right')

    # Recall lineplot
    sns.lineplot(
        data=data,
        x="Muestra",
        y="Recall",
        hue="Algoritmo",
        marker="o",
        ax=axes[1]
    )
    axes[1].set_title("Recall por Réplica")
    axes[1].set_xlabel("Réplica")
    axes[1].set_ylabel("Recall")
    axes[1].tick_params(axis='x', rotation=0)
    axes[1].get_legend().remove() 
    plt.tight_layout()
    plt.show()



def plot_metrics_boxplot_and_lines(data, metrics, algorithm_order, replica_order):
    """
    Plot a 3x2 figure: each row represents a metric (Precision, Recall, Fscore),
    with a boxplot and lineplot for each metric.
    Boxplot: x=Muestra, hue=Algoritmo.
    Lineplot: x=Muestra, line=Algoritmo.
    """
    sns.set_style("whitegrid")
    
    # Ensure 'Algoritmo' and 'Muestra' are categorical and ordered
    data['Algoritmo'] = pd.Categorical(data['Algoritmo'], categories=algorithm_order, ordered=True)
    data['Muestra'] = pd.Categorical(data['Muestra'], categories=replica_order, ordered=True)
    
    # Create subplots
    fig, axes = plt.subplots(len(metrics), 2, figsize=(15, 12))
    
    for i, metric in enumerate(metrics):
        # Boxplot
        sns.boxplot(
            data=data,
            x="Muestra",
            y=metric,
            hue="Algoritmo",
            ax=axes[i, 0],
            palette="bright"
        )
        axes[i, 0].set_title(f"{metric} - Boxplot (por réplica)")
        axes[i, 0].set_xlabel("Muestra")
        axes[i, 0].set_ylabel(metric)
        axes[i, 0].tick_params(axis='x', rotation=45)
        if i == 0:
            axes[i, 0].legend(title="Algoritmo", loc='upper left')
        else:
            axes[i, 0].legend_.remove()

        # Lineplot
        sns.lineplot(
            data=data,
            x="Muestra",
            y=metric,
            hue="Algoritmo",
            marker="o",
            ax=axes[i, 1]
        )
        axes[i, 1].set_title(f"{metric} - Lineplot")
        axes[i, 1].set_xlabel("Muestra")
        axes[i, 1].set_ylabel(metric)
        axes[i, 1].tick_params(axis='x', rotation=0)
        if i == 0:
            axes[i, 1].legend( loc='upper right')
        else:
            axes[i, 1].legend_.remove()

    plt.tight_layout()
    plt.show()




def main(filepath):
    """
    Main function to load data and generate all plots
    """
    # load data
    data = load_data(filepath)
    # define orders and metrics
    algorithm_order= ["cnmops", "cnvkit", "contra", "exomedepth", "manta", "laconv", "xhmm"]
    kit_order = ["IDT-V1", "IDT-V2", "ROCHE-V1"]
    replica_order = ["Réplica 1", "Réplica 2", "Réplica 3", "Réplica 4", "Réplica 5"]
    pipeline_order = ["bowtie2-picard-gatk3", "bowtie2-samtools-gatk4", "minimap2-samtools-gatk4"]
    metrics = ["Precision", "Recall", "Fscore"]
    # generate plots
    plot_boxplots_by_algorithm_by_goldstandard(data, algorithm_order, metrics)
    plot_lines_by_algorithm_by_goldstandard(data, algorithm_order)
    plot_lines_by_algorithm_by_pipeline(data, pipeline_order)
    plot_lines_by_algorithm_by_replica(data, replica_order)
    #plot_metrics_boxplot_and_lines(data, metrics, algorithm_order, replica_order)

# execute the script
if __name__ == "__main__":
    filepath = "/Users/claulara/Desktop/TFM/results_posprocess/metrics/overall_metrics.csv"
    main(filepath)
