#!/usr/bin/python

import optparse
import pandas as pd
import gzip
import re
    

def parse_bed_file(file_path):

    with gzip.open(file_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', header=None, usecols=[0, 1, 2, 3],
                         names=['chr', 'start', 'end', 'mean_cov'], low_memory=False)
        df['chr'] = df['chr'].apply(lambda x: f'chr{x}' if not str(x).startswith('chr') else str(x))
    return df


def sort_key(region):
    """
    Función para ordenar las regiones cromosómicas en formato 'chr:start-end'.
    Ordena primero por cromosoma numérico, luego 'chrX', y luego 'chrY'.
    """
    chrom, positions = region.split(':')
    start, end = map(int, positions.split('-'))
    
    # Convertir el nombre del cromosoma en una clave para ordenar.
    # Cromosomas numéricos se mantienen como enteros, 'chrX' y 'chrY' se colocan al final.
    if chrom == 'chrX':
        chrom_key = 23
    elif chrom == 'chrY':
        chrom_key = 24
    else:
        chrom_key = int(chrom[3:])  # Asumimos que todos los otros cromosomas son numéricos
    
    return (chrom_key, start)


def generate_cov_mat(file_paths):

    '''
    Función para generar una matriz de cobertura a partir de una lista de archivos .bed.gz.
    
    Input:
    - file_paths: Lista de rutas de archivos .bed.gz, cada uno correspondiente a un sample.
    Output:
    - cov_mat: DataFrame con samples como filas, regiones como columnas, y cobertura media como valores.
    '''

    coverage_data = {}
    all_regions = set()  # Set con los índices de cada región 

    for file_path in file_paths:

        sample_name = file_path.split('/')[-1].replace('.regions.bed.gz', '')
        df = parse_bed_file(file_path)
        
        # Crear los índices de cada región cubierta en formato 'chr:start-end'
        df['region'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
        # Set con los índices de cada región 
        all_regions.update(df['region'])
        # Crear diccionario de diccionarios para cada muestra con las coberturas medias cubiertas
        coverage_data[sample_name] = df.set_index('region')['mean_cov'].to_dict()

    # Convertir el set de todas las regiones en una lista ordenada 
    all_regions = sorted(list(all_regions), key=sort_key)
    
    # Crear el DataFrame de cobertura
    cov_mat = pd.DataFrame.from_dict(coverage_data, orient='index', columns=all_regions)

    return cov_mat


def run(argv=None):

    parser = optparse.OptionParser()
    parser.add_option('--b', default=None, help='Comma-separated list of bed files with the region of interest', dest='bed_files')
    parser.add_option('--o', default=None, help='Path of output coverage matrix', dest='output_path')
    (options, args) = parser.parse_args(argv[1:])


    bed_files = options.bed_files.split(',') if options.bed_files else []
    print(bed_files)
    output_path = options.output_path


    cov_mat = generate_cov_mat(file_paths=bed_files)
    cov_mat.to_csv(output_path, sep='\t')

    print(f'Matriz de cobertura generada y guardada en {output_path}')


if __name__ == '__main__':
    import sys
    run(sys.argv)

