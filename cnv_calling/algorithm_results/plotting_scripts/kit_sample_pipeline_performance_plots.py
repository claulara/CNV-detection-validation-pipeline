import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def load_data(filepath, fillna_value=0):
    data = pd.read_csv(filepath)
    data.fillna(0, inplace=True)
    return data


def plot_boxplots_by_kit(data, kit_order, metrics):
    """
    Generate boxplots for metrics, grouped by capture kit and gold standard.
    """
    sns.set_style("whitegrid")    
    # generate subplots
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 5))
    for ax, metric in zip(axes, metrics):
        sns.boxplot(
            x='kitdecaptura',
            y=metric,
            hue='gold_standard',
            data=data,
            ax=ax,
            palette="bright",
            order=kit_order)
        ax.set_title(f"{metric}")
        ax.set_xlabel("Kit de exoma")
        ax.set_ylabel(metric)
        ax.legend(title="Condition", loc='upper right')
    # remove duplicate legends
    for ax in axes[1:]:
        ax.legend().remove()
    plt.tight_layout()
    plt.show()


def plot_lines_by_kit_by_goldstandard(data, kit_order):
    """
    Generate line plots for metrics, grouped by capture kit and gold standard.
    """
    sns.set_style("whitegrid")
    custom_palette = {'general': 'blue', 'pass': 'orange'}
    # transform data to a long format for lineplot and sort by kit
    df = pd.melt(
        data,
        id_vars=['kitdecaptura', 'gold_standard'],
        value_vars=['Precision', 'Recall'],
        var_name='Métrica',
        value_name='Valor')
    df['kitdecaptura'] = pd.Categorical(df['kitdecaptura'], categories=kit_order, ordered=True)
    # generate plot
    plt.figure(figsize=(12, 6))
    sns.lineplot(
        data=df,
        x='kitdecaptura',
        y='Valor',
        hue='gold_standard',
        style='Métrica',
        markers=True,
        dashes={'Precision': '', 'Recall': (2, 2)},
        palette=custom_palette)
    plt.xlabel('Kit de exoma', fontsize=11)
    plt.ylabel('Valor', fontsize=11)
    plt.xticks(fontsize=12,rotation=0)
    plt.legend(fontsize=11, loc='upper center')
    plt.tight_layout()
    plt.show()


def plot_boxplots_by_sample(data, replica_order, metrics):
    """
    Generate boxplots for metrics, grrouped by sample.
    """
    sns.set_style("whitegrid")
    # generate subplots 
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 5))
    for ax, metric in zip(axes, metrics):
        sns.boxplot(
            x="Muestra",
            y=metric,
            data=data,
            palette="Set2",
            order=replica_order,
            ax=ax)
        ax.set_title(f"{metric}")
        ax.set_ylabel(metric)
        ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.show()


def plot_lines_by_sample(data, replica_order):
    """
    Generate line plots for metrics, grouped by sample.
    """
    sns.set_style("whitegrid")
    sns.set_style("whitegrid")
    # generate plot
    plt.figure(figsize=(12, 6))
    sns.lineplot(data=data, x='Muestra', y='Precision', label='Precision', marker='o', color='blue')
    sns.lineplot(data=data, x='Muestra', y='Recall', label='Recall', marker='o', color='orange')
    plt.ylabel('Valor', fontsize=10)
    # replica's order
    plt.xticks(ticks=range(len(replica_order)), fontsize=12, labels=replica_order, rotation=0)
    plt.legend(title='Métricas', fontsize=11)
    plt.tight_layout()
    plt.show()


def plot_boxplots_by_sample_by_goldstandard(data, replica_order, metrics):
    """
    Generate boxplots for metrics, grouped by sample and gold standard.
    """
    sns.set_style("whitegrid") 
    # generate subplots
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 5))
    for ax, metric in zip(axes, metrics):
        sns.boxplot(
            x="Muestra",
            y=metric,
            hue="gold_standard",
            data=data,
            palette="bright",
            order=replica_order,
            ax=ax)
        ax.set_title(f"{metric}")
        ax.set_ylabel(metric)
        ax.tick_params(axis='x', rotation=45)
        ax.legend(title="Condition", loc='upper right')
    # remove duplicate legends
    for ax in axes[1:]:
        ax.legend().remove()
    plt.tight_layout()
    plt.show()


def plot_lines_by_sample_by_goldstandard(data, replica_order):
    """
    Generate line plots for metrics, grouped by sample and gold standard.
    """
    sns.set_style("whitegrid")
    custom_palette = {'general': 'blue', 'pass': 'orange'}
    # transform data to a long format for lineplot and sort by sample
    df = pd.melt(
        data,
        id_vars=['Muestra', 'gold_standard'],
        value_vars=['Precision', 'Recall'],
        var_name='Métrica',
        value_name='Valor')
    df['Muestra'] = pd.Categorical(df['Muestra'], categories=replica_order, ordered=True)
    # generate plot
    plt.figure(figsize=(12, 6))
    sns.lineplot(
        data=df,
        x='Muestra',
        y='Valor',
        hue='gold_standard',
        style='Métrica',
        markers=True,
        dashes={'Precision': '', 'Recall': (2, 2)},
        palette=custom_palette)
    plt.xlabel('Muestra', fontsize=11)
    plt.ylabel('Valor', fontsize=11)
    plt.xticks(fontsize=12, rotation=0)
    plt.legend(fontsize=11, loc='upper center')
    plt.tight_layout()
    plt.show()


def plot_boxplots_by_pipeline(data, pipeline_order, metrics):
    """
    Generate boxplots for metrics, grouped by pipeline.
    """
    sns.set_style("whitegrid") 
    custom_palette = [
        (0.302, 0.686, 0.290),  # green
        (0.596, 0.306, 0.639),  # purple
        (1.000, 1.000, 0.200)   # yellow
    ]
    # generate subplots
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 5))
    for ax, metric in zip(axes, metrics):
        sns.boxplot(
            x="Preprocesado",
            y=metric,
            hue="gold_standard",
            data=data,
            palette=custom_palette,
            order=pipeline_order,
            ax=ax)
        ax.set_title(f"{metric}")
        ax.set_xlabel("Pipeline Bioinformática")
        ax.set_ylabel(metric)
        ax.tick_params(axis='x', rotation=45)
        ax.legend(title="Condition", loc='upper right')
    # remove duplicate legends
    for ax in axes[1:]:
        ax.legend().remove()
    plt.tight_layout()
    plt.show()


def plot_boxplots_by_sample_pipeline(data, replica_order, pipeline_order, metrics):
    """
    Generate boxplots for metrics, grouped by sample and pipeline.
    """
    sns.set_style("whitegrid") 
    custom_palette = [
        (0.302, 0.686, 0.290),  # green
        (0.596, 0.306, 0.639),  # purple
        (1.000, 1.000, 0.200)   # yellow
    ]
    # generate subplots
    fig, axes = plt.subplots(1, len(metrics), figsize=(6 * len(metrics), 6), sharey=False)
    for i, metric in enumerate(metrics):
        ax = axes[i]
        sns.boxplot(
            x='Muestra',
            y=metric,
            hue='Preprocesado',
            data=data,
            palette=custom_palette,  
            order=replica_order,
            hue_order=pipeline_order,
            ax=ax)
        ax.set_title(f"{metric.capitalize()}")
        ax.set_xlabel("Réplica")
        ax.set_ylabel(metric.capitalize() if i == 0 else "")
        ax.tick_params(axis='x', rotation=45)
        ax.legend(title="Pipeline Bioinfomática")
    # remove duplicate legends
    for ax in axes[1:]:
        ax.legend().remove()
    plt.tight_layout()
    plt.show()


def plot_lines_by_pipeline_by_goldstandard(data, pipeline_order):
    """
    Generate line plots for metrics, grouped by pipeline and gold standard.
    """
    sns.set_style("whitegrid")
    custom_palette = {'general': 'blue', 'pass': 'orange'}
    # transform data to a long format for lineplot and sort by sample
    df = pd.melt(
        data,
        id_vars=['Preprocesado', 'gold_standard'],
        value_vars=['Precision', 'Recall'],
        var_name='Métrica',
        value_name='Valor')
    df['Preprocesado'] = pd.Categorical(df['Preprocesado'], categories=pipeline_order, ordered=True)
    # generate plot
    plt.figure(figsize=(12, 6))
    sns.lineplot(
        data=df,
        x='Preprocesado',
        y='Valor',
        hue='gold_standard',
        style='Métrica',
        markers=True,
        dashes={'Precision': '', 'Recall': (2, 2)},
        palette=custom_palette)
    plt.xlabel('Pipeline Bioinformática', fontsize=11)
    plt.ylabel('Valor', fontsize=11)
    plt.xticks(fontsize=12, rotation=0)
    plt.legend(fontsize=11, loc='upper left')
    plt.tight_layout()
    plt.show()
    

def plot_boxplots_by_sample_pipeline_by_goldstandard(data, replica_order, pipeline_order, metrics):
    """
    Generate boxplots for metrics, grouped by sample and pipeline.
    """
    sns.set_style("whitegrid") 
    custom_palette = [
        (0.302, 0.686, 0.290),  # green
        (0.596, 0.306, 0.639),  # purple
        (1.000, 1.000, 0.200)   # yellow
    ]
    # generate subgroups by gold standard
    conditions = data['gold_standard'].unique()
    num_conditions = len(conditions)
    # generate subplots
    fig, axes = plt.subplots(num_conditions, len(metrics), figsize=(5 * len(metrics), 5 * num_conditions))
    for i, condition in enumerate(conditions):
        condition_data = data[data['gold_standard'] == condition]
        for j, metric in enumerate(metrics):
            ax = axes[i, j] if num_conditions > 1 else axes[j]
            sns.boxplot(
                x='Muestra',
                y=metric,
                hue='Preprocesado',
                data=condition_data,
                palette=custom_palette,
                order=replica_order,
                hue_order=pipeline_order,
                ax=ax
            )
            ax.set_title(f"{metric} - {condition.capitalize()}")
            ax.set_ylabel(metric if j == 0 else "")
            ax.tick_params(axis='x', rotation=45)
            if i == 0 and j == 0:
                ax.legend(title="Pipeline", loc='upper right')
            else:
                ax.legend().remove()
    plt.tight_layout()
    plt.show()


def heatmap_by_gold_standard(data, metrics):
    """
    Generate heatmap for metrics, grouped by sample and gold standard.
    """
    gold_standards = data["gold_standard"].unique()
    fig, axes = plt.subplots(1, len(gold_standards), figsize=(20, 5), sharey=True)
    for i, gold_standard in enumerate(gold_standards):
        subset = data[data["gold_standard"] == gold_standard]
        heatmap_data = subset.groupby("Muestra")[metrics].mean().T
        sns.heatmap(heatmap_data, annot=True, cmap="viridis", cbar_kws={'label': 'Media'}, ax=axes[i])
        axes[i].set_title(f"{gold_standard}")
        axes[i].set_ylabel("Métrica")
        axes[i].set_xlabel("Muestra")
    plt.tight_layout()
    plt.show()


def heatmaps_by_svtype_and_gold_standard(data, metrics, svtype_col, gold_standard_col, title_prefix="Heatmap"):
    svtypes = data["svtype"].unique()
    gold_standards = data["gold_standard"].unique()
    fig, axes = plt.subplots(len(svtypes), len(gold_standards), figsize=(16, 16), sharey=True)
    for i, svtype in enumerate(svtypes):
        for j, gold_standard in enumerate(gold_standards):
            subset = data[(data[svtype_col] == svtype) & (data[gold_standard_col] == gold_standard)]
            heatmap_data = subset.groupby("mMestra")[metrics].mean().T
            sns.heatmap(heatmap_data, annot=True, cmap="viridis", cbar_kws={'label': 'Media'}, ax=axes[i, j])
            axes[i, j].set_title(f"{title_prefix}: {svtype} | Gold standard {gold_standard}")
            axes[i, j].set_ylabel("Métrica")
            axes[i, j].set_xlabel("Muestra")
    plt.tight_layout()
    plt.show()


def main(filepath):
    """
    Main function to load data and Generate all plots
    """
    # load data
    data = load_data(filepath)
    sv_data = load_data(filepath_sv)
    # define orders and metrics
    kit_order = ["IDT-V1", "IDT-V2", "ROCHE-V1"]
    replica_order = ["Réplica 1", "Réplica 2", "Réplica 3", "Réplica 4", "Réplica 5"]
    pipeline_order = ["bowtie2-picard-gatk3", "bowtie2-samtools-gatk4", "minimap2-samtools-gatk4"]
    metrics = ["Precision", "Recall", "Fscore"]
    # generate plots
    #plot_boxplots_by_kit(data, kit_order, metrics)
    #plot_lines_by_kit_by_goldstandard(data, kit_order)
    plot_boxplots_by_sample(data, replica_order, metrics)
    plot_lines_by_sample(data, replica_order)
    plot_boxplots_by_sample_by_goldstandard(data, replica_order, metrics)
    plot_lines_by_sample_by_goldstandard(data, replica_order)
    plot_boxplots_by_pipeline(data, pipeline_order, metrics)
    plot_boxplots_by_sample_pipeline(data, replica_order, pipeline_order, metrics)
    plot_lines_by_pipeline_by_goldstandard(data, pipeline_order)
    plot_boxplots_by_sample_pipeline_by_goldstandard(data, replica_order, pipeline_order, metrics)
    heatmap_by_gold_standard(data, metrics)
    #heatmaps_by_svtype_and_gold_standard(sv_data, metrics)

# execute the script
if __name__ == "__main__":
    filepath = "/Users/claulara/Desktop/TFM/results_posprocess/metrics/overall_metrics.csv"
    filepath_sv = "/Users/claulara/Desktop/TFM/results_posprocess/metrics/svtype_metrics.csv"
    main(filepath)



