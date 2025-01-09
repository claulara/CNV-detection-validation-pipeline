import os
import json
import pandas as pd


def extract_metrics(json_file):
    """
    Function to extract metrics from a JSON file
    """
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        # extract information of interest from the file path which is not into the JSON file
        path_parts = json_file.split(os.sep)
        algoritmo = path_parts[-6]
        kitdecaptura = path_parts[-5]
        preprocesado = path_parts[-4]
        muestra = path_parts[-3]
        gold_standard = path_parts[-2]

        # overall metrics from the JSON file
        overall_stats = data['PerSampleStats'][0]['OverallStats'][0]
        # deletion and duplication metrics from the JSON file
        deletion_stats = next(item for item in data['PerSampleStats'][0]['DetailedStats'] if item['VariantType'] == 'Deletion')['OverallStats'][0]
        duplication_stats = next(item for item in data['PerSampleStats'][0]['DetailedStats'] if item['VariantType'] == 'Duplication')['OverallStats'][0]

        # overall
        metrics_overall = {
            'algoritmo': algoritmo,
            'kitdecaptura': kitdecaptura,
            'preprocesado': preprocesado,
            'muestra': muestra,
            'gold_standard': gold_standard,
            'Precision': overall_stats['Precision'],
            'Recall': overall_stats['Recall'],
            'Fscore': overall_stats['Fscore'],
            'QueryTpCount': overall_stats['QueryTpCount'],
            'QueryFpCount': overall_stats['QueryFpCount'],
            'QueryTotalCount': overall_stats['QueryTotalCount'],
            'TruthTpCount': overall_stats['TruthTpCount'],
            'TruthFnCount': overall_stats['TruthFnCount'],
            'TruthTotalCount': overall_stats['TruthTotalCount']
        }
        
        # deletions
        metrics_deletion = metrics_overall.copy()
        metrics_deletion.update({
            'Precision': deletion_stats['Precision'],
            'Recall': deletion_stats['Recall'],
            'Fscore': deletion_stats['Fscore'],
            'QueryTpCount': deletion_stats['QueryTpCount'],
            'QueryFpCount': deletion_stats['QueryFpCount'],
            'QueryTotalCount': deletion_stats['QueryTotalCount'],
            'TruthTpCount': deletion_stats['TruthTpCount'],
            'TruthFnCount': deletion_stats['TruthFnCount'],
            'TruthTotalCount': deletion_stats['TruthTotalCount'],
            'svtype': 'Deletion'  # Specify this is for deletions
        })
        
        # duplications
        metrics_duplication = metrics_overall.copy()
        metrics_duplication.update({
            'Precision': duplication_stats['Precision'],
            'Recall': duplication_stats['Recall'],
            'Fscore': duplication_stats['Fscore'],
            'QueryTpCount': duplication_stats['QueryTpCount'],
            'QueryFpCount': duplication_stats['QueryFpCount'],
            'QueryTotalCount': duplication_stats['QueryTotalCount'],
            'TruthTpCount': duplication_stats['TruthTpCount'],
            'TruthFnCount': duplication_stats['TruthFnCount'],
            'TruthTotalCount': duplication_stats['TruthTotalCount'],
            'svtype': 'Duplication'  # Specify this is for duplications
        })
        
        return metrics_overall, metrics_deletion, metrics_duplication
    
    except Exception as e:
        print(f"Error processing file {json_file}: {e}")
        return None, None, None



def main(main_dir):
    """
    Function to process JSON files and generate CSVs
    """
    main_dir = main_dir
    overall_metrics = []
    detailed_metrics = []

    # find JSON files from witty.er
    for subdir, dirs, files in os.walk(main_dir):
        for file in files:
            if file == 'Wittyer.Stats.json':  # specific JSON file from witty.er that contains the information of interest
                json_file = os.path.join(subdir, file)
                print(f"Processing file {json_file}")
                metrics_overall, metrics_deletion, metrics_duplication = extract_metrics(json_file)
                if metrics_overall:
                    overall_metrics.append(metrics_overall)
                if metrics_deletion and metrics_duplication:
                    detailed_metrics.extend([metrics_deletion, metrics_duplication])

    # overall CSV
    if overall_metrics:
        df_overall = pd.DataFrame(overall_metrics)
        df_overall = df_overall.fillna(0)
        # replace sample names with replica names
        df_overall['muestra'] = df_overall['muestra'].str.replace('_align.realign', '', regex=False)
        replicas = {
            'NGS25906-HiSeq4000-exoma-Run200604-HG-0001_S287': 'Réplica 1',
            'NGS25603-HiSeq4000-exoma-Run200513-HG-0001_S358': 'Réplica 2',
            'NGS35137-Novaseq6000-Nextera-DNA-Exome-Panel-v2-IDT-Run220707-Control-NIST_S4': 'Réplica 3',
            'NGS39574-Novaseq6000-Nextera-DNA-Exome-Panel-v2-IDT-Run230711-HG-001-NIST_S16': 'Réplica 4',
            'NGS38792-Novaseq6000-Exoma-Run230512-CONTROL-HG001_S26': 'Réplica 5'
        }
        df_overall['muestra'] = df_overall['muestra'].replace(replicas)
        df_overall.to_csv('/ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/results_posprocess/metrics/overall_metrics.csv', index=False)

    # combined deletions and duplications CSV
    if detailed_metrics:
        df_detailed = pd.DataFrame(detailed_metrics)
        df_detailed = df_detailed.fillna(0)
        df_detailed['muestra'] = df_detailed['muestra'].str.replace('_align.realign', '', regex=False)
        # replace sample names with replica names
        df_detailed['muestra'] = df_detailed['muestra'].replace(replicas)
        df_detailed.to_csv('/ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/results_posprocess/metrics/svtype_metrics.csv', index=False)

    print("CSV files successfully created.")



if __name__ == "__main__":
    main(main_dir = '/ingemm/scratch/TFM/CNV/TFM_borrador/algorithms/results_posprocess')
