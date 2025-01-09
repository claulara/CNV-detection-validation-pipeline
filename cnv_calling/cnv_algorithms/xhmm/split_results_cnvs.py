import optparse
import pandas as pd
import os


def split_variants_by_sample(cnv_file, output_path):
    
    df = pd.read_csv(cnv_file, sep='\t')
    os.makedirs(output_path, exist_ok=True)


    for sample_name, sample_df in df.groupby(df.columns[0]):
        # carpeta para cada sample
        sample_dir = os.path.join(output_path, sample_name)
        os.makedirs(sample_dir, exist_ok=True)
        # .csv con las cnvs de cada sample
        output_file = os.path.join(sample_dir, f"{sample_name}_cnv_output.csv")
        sample_df.to_csv(output_file, sep='\t', index=False, header=True)



def run(argv=None):

    parser = optparse.OptionParser()
    parser.add_option('--f', default=None, help='CSV with all samples and their cnvs', dest='cnv_file')
    parser.add_option('--o', default=None, help='Path of output coverage matrix', dest='output_path')
    (options, args) = parser.parse_args(argv[1:])

    cnv_file = options.cnv_file
    output_path = options.output_path
    split_variants_by_sample(cnv_file, output_path)

    print(f'Variantes guardadas en {output_path}')


if __name__ == '__main__':
    import sys
    run(sys.argv)