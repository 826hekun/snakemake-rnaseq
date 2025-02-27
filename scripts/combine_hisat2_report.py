import os
import re
import argparse
import pandas as pd

def extract_stats(log_file):
    with open(log_file, 'r') as f:
        log_content = f.read()

    uniq_mapped_reads = int(re.search(r'Aligned concordantly 1 time: (\d+)', log_content).group(1))
    uniq_mapped_rate = float(re.search(r'Aligned concordantly 1 time: \d+ \((\d+\.\d+)%\)', log_content).group(1))

    multi_mapped_reads = int(re.search(r'Aligned concordantly >1 times: (\d+)', log_content).group(1))
    multi_mapped_rate = float(re.search(r'Aligned concordantly >1 times: \d+ \((\d+\.\d+)%\)', log_content).group(1))

    overall_mapped_rate = float(re.search(r'Overall alignment rate: (\d+\.\d+)%', log_content).group(1))

    return {
        'UniqMappedReads': uniq_mapped_reads,
        'UniqMappedRate': uniq_mapped_rate,
        'MultiMappedReads': multi_mapped_reads,
        'MultiMappedRate': multi_mapped_rate,
        'OverallMappedRate': overall_mapped_rate
    }

parser = argparse.ArgumentParser(description='Combine hisat2 log files into a summary table')
parser.add_argument('-i', '--input_folder', type=str, required=True, help='Path to the folder containing hisat2.log files')
parser.add_argument('-o', '--output_folder', type=str, default='report', help='Path to the output folder [default: report]')

args = parser.parse_args()

if not os.path.exists(args.input_folder):
    raise ValueError("Error: Input folder does not exist. Please provide a valid path.")

log_files = [os.path.join(args.input_folder, f) for f in os.listdir(args.input_folder) if f.endswith('.hisat2.log')]

sample_names = [re.search(r'^(.*?)\.hisat2\.log$', os.path.basename(log_file)).group(1) for log_file in log_files]

hisat2_summary = pd.DataFrame([extract_stats(log_file) for log_file in log_files])
hisat2_summary.insert(0, 'Sample', sample_names)

hisat2_summary = hisat2_summary.sort_values(by='Sample')

if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)

hisat2_summary.to_csv(os.path.join(args.output_folder, 'hisat2_summary.tsv'), sep='\t', index=False)

