import os
import matplotlib.pyplot as plt
import pandas as pd 
import seaborn as sns
import numpy as np  


def df_bed(path):
    """
    Read bed and create a dataframe with the specific columns
    """
    columnas = ["Chromosome", "Start", "End", "Type", "Length", "Quality", "Support", "Genotype"]
    data = pd.read_csv(path, sep='\t', names=columnas)
    data.fillna(0, inplace=True)
    sample_name = os.path.basename(path)
    return data, sample_name


def sv_chrom_distrib(data, sample_name):
    """
    Plots the distribution of structural variants across chromosomes
    """
    plt.figure(figsize=(15,6))
    sns.countplot(x="Chromosome", data=data)
    plt.title(f"Distribución de variantes en los cromosomas - {sample_name}")
    plt.xlabel("Cromosoma")
    plt.ylabel("Número de variantes")
    plt.show()


def sv_type_distrib(data, sample_name, ax):
    """
    Plots a pie chart showing the distribution of variant types
    """
    type_counts = data['Type'].value_counts()
    wedges, texts, autotexts = ax.pie(
        type_counts,
        autopct='%1.1f%%',  
        startangle=140,
        colors=sns.color_palette("pastel"),
        wedgeprops={'edgecolor': 'black'},  
    )
    ax.legend(wedges, type_counts.index, title="Tipo de variante", loc="center left", bbox_to_anchor=(1, 0.5), fontsize=14, title_fontsize=16)  
    ax.axis('equal')


def sv_qual_distrib(data, sample_name, ax):
    """
    Plots a pie chart showing the distribution of variant quality
    """
    type_counts = data['Quality'].value_counts()
    wedges, texts, autotexts = ax.pie(
        type_counts,
        autopct='%1.1f%%',  
        startangle=140,
        colors=sns.color_palette("pastel"),
        wedgeprops={'edgecolor': 'black'},  
    )
    ax.legend(wedges, type_counts.index, title="Calidad", loc="center left", bbox_to_anchor=(1, 0.5), fontsize=14, title_fontsize=16)  
    ax.axis('equal')


def sv_genotype_distrib(data, sample_name, ax):
    """
    Plots a pie chart showing the distribution of variant genotypes
    """
    type_counts = data['Genotype'].value_counts()
    wedges, texts, autotexts = ax.pie(
        type_counts,
        autopct='%1.1f%%',  
        startangle=140,
        colors=sns.color_palette("pastel"),
        wedgeprops={'edgecolor': 'black'}, 
    )
    ax.legend(wedges, type_counts.index, title="Genotipo", loc="center left", bbox_to_anchor=(1, 0.5), fontsize=14, title_fontsize=16)  
    ax.axis('equal')


def plot_three_pie_charts(data, sample_name):
    """
    Creates three pie charts in a single row: distribution by type, quality, and genotype
    """
    fig, axes = plt.subplots(1, 3, figsize=(20, 7))  
    sv_type_distrib(data, sample_name, axes[0])
    sv_qual_distrib(data, sample_name, axes[1])
    sv_genotype_distrib(data, sample_name, axes[2])
    plt.tight_layout()
    plt.show()

def sv_length_distrib(data, sample_name):
    """
    Plots the frequency distribution of variant lengths by intervals.
    """
    bins = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 
    5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000,
    70000, 80000, 90000, 100000, 500000, 1000000, 10000000]
    bin_labels = ['0-100', '101-200', '201-300', '301-400', '401-500', '501-600', 
                '601-700', '701-800', '801-900', '901-1000', '1001-2000', '2001-3000', 
                '3001-4000', '4001-5000', '5001-6000', '6001-7000', '7001-8000', 
                '8001-9000', '9001-10000', '10001-20000', '20001-30000', '30001-40000', 
                '40001-50000', '50001-60000', '60001-70000', '70001-80000', '80001-90000', 
                '90001-100000', '100001-500000',  '500001-1000000', '1000001-10000000']
    
    freq, _ = np.histogram(data['Length'], bins=bins) # variant frequencies per bin
    x_positions = np.arange(len(freq))

    plt.figure(figsize=(12, 6))
    plt.bar(x_positions, freq, color='skyblue', edgecolor='black')
    plt.xticks(x_positions, bin_labels, rotation=45, ha='right')
    plt.xlabel('Tamaño de las variantes (bp)')
    plt.ylabel('Frecuencia')
    plt.title(f"Distribución de tamaño de variantes - {sample_name}")
    plt.tight_layout()
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.show()


def sv_length_distrib_svtype(data, sample_name):
    """
    Plots the frequency distribution of variant lengths by intervals and by variant type
    """
    bins = [0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 
            5000, 6000, 7000, 8000, 9000, 10000, 20000, 40000, 
            60000, 80000, 100000, 500000, 1000000, 10000000]

    # create a new df to count variants by type and bin
    data['Bin'] = pd.cut(data['Length'], bins=bins)
    grouped = data.groupby(['Bin', 'Type']).size().unstack(fill_value=0)

    plt.figure(figsize=(12, 6))
    width = 0.3  # Width of the bars
    x_positions = np.arange(len(grouped))

    # Palette colors for specific variant types
    colors = {"DEL": "blue", "INS": "orange", "INV": "green"}
    palette = sns.color_palette("Set2", n_colors=len(grouped.columns))

    # create a bar for each variant type
    for i, column in enumerate(grouped.columns):
        color = colors.get(column, palette[i % len(palette)])
        plt.bar(x_positions + i * width, grouped[column], width=width, label=column, color=color)

    plt.xlabel('Intervalos por tamaño de variante (bp)')
    plt.ylabel('Frecuencia')
    plt.xticks(x_positions + width * (len(grouped.columns) - 1) / 2, [str(interval) for interval in grouped.index])
    plt.xticks(rotation=45, ha='right')
    plt.title(f'Distribución de las variantes low quality filtradas del estudio {sample_name.split("_")[0]} en función de su tamaño')
    plt.legend(title='Tipo de Variante')
    plt.tight_layout()
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.show()



beds = {"complete_bed": ["/Users/claulara/Desktop/TFM/gold_standard/benchmarks/1000genomes/1000genomes.bed",
                         "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/giab/giab.bed", 
                         "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/metasv/metasv.bed",
                         "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/svclassify/svclassify.bed"], 
                         
        "pass_bed": ["/Users/claulara/Desktop/TFM/gold_standard/benchmarks/1000genomes/PASS/1000genomes_pass.bed", 
                     "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/giab/PASS/giab_pass.bed", 
                     "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/metasv/PASS/metasv_pass.bed"], 

        "pass_filtered_bed": ["/Users/claulara/Desktop/TFM/gold_standard/benchmarks/1000genomes/PASS/1000genomes_pass_filtered.bed", 
                      "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/giab/PASS/giab_pass_filtered.bed", 
                      "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/metasv/PASS/metasv_pass_filtered.bed",
                      "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/svclassify/svclassify_filtered.bed"], 

        'low_qual_bed': ["/Users/claulara/Desktop/TFM/gold_standard/benchmarks/1000genomes/LowQual/1000genomes_lowqual.bed", 
                         "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/giab/LowQual/giab_lowqual.bed",
                         "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/metasv/LowQual/metasv_lowqual.bed"],

        'low_qual_filtered_bed': ["/Users/claulara/Desktop/TFM/gold_standard/benchmarks/1000genomes/LowQual/1000genomes_lowqual_filtered.bed", 
                         "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/giab/LowQual/giab_lowqual_filtered.bed",
                         "/Users/claulara/Desktop/TFM/gold_standard/benchmarks/metasv/LowQual/metasv_lowqual_filtered.bed"],
                         }



for path in beds["complete_bed"]:
    data, sample_name = df_bed(path)
    print(sample_name)
    print(data['Length'].describe())
    print("")
    sv_chrom_distrib(data, sample_name)
    plot_three_pie_charts(data, sample_name)
    sv_length_distrib(data, sample_name)
    sv_length_distrib_svtype(data, sample_name)



