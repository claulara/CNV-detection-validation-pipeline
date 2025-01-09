import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def load_data(filepath, fillna_value=0):
    data = pd.read_csv(filepath)
    data.fillna(fillna_value, inplace=True)
    data['svtype'] = data['svtype'].astype('category')
    return data


def heatmap_by_goldstandard_bysvtype(data, output_dir):
    # Diccionario para traducir svtype
    svtype_translation = {
        "Deletion": "Deleciones",
        "Duplication": "Duplicaciones"}
    # combinatios of gold_standard and svtype
    combinations = [
        {"gold_standard": "pass", "svtype": "Deletion"},
        {"gold_standard": "pass", "svtype": "Duplication"},
        {"gold_standard": "general", "svtype": "Deletion"},
        {"gold_standard": "general", "svtype": "Duplication"}]
    for combo in combinations:
        gold_standard = combo["gold_standard"]
        svtype = combo["svtype"]
        # create subdataframes for the combinations
        filtered_data = data[(data["gold_standard"] == gold_standard) & (data["svtype"] == svtype)]
        # create subplots
        fig, axes = plt.subplots(2, 1, figsize=(10, 15), sharey=False)
        metrics = ["Precision", "Recall"]
        for ax, metric in zip(axes, metrics):
            pivot_data = filtered_data.pivot_table(
                values=metric, 
                index=["Algoritmo", "Muestra"], 
                columns="Preprocesado", 
                aggfunc='mean')
            sns.heatmap(pivot_data, annot=True, cmap="viridis", ax=ax)
            # change the sv type name for the title 
            translated_svtype = svtype_translation.get(svtype, svtype)  # Traduce si es posible
            ax.set_title(f"{metric} en {translated_svtype} para Gold Standard {gold_standard}")
            ax.set_ylabel("Algoritmo y Muestra")
            ax.set_xlabel("Preprocesado")
        plt.tight_layout()
        filename = f"{output_dir}results_heatmap_{gold_standard}_{svtype}.png"
        plt.savefig(filename)  
        plt.show()
        plt.close()  

data = load_data(filepath)
filepath = "/Users/claulara/Desktop/TFM/results_posprocess/metrics/svtype_metrics.csv"
output_dir = "/Users/claulara/Desktop/TFM/imagenes/"
heatmap_by_goldstandard_bysvtype(data, output_dir)



