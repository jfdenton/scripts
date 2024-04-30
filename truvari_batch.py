import os
import glob
import csv
import subprocess
import json

###
### samples.txt should be a list of samples to be compared, e.g.  vcf_path1    vcf_path2   replicate_group
###     samples.txt rows should contain tab-separated tabix zipped/index vcf files and any string name
### parameters.txt should be raw values for the 4 parameters, e.g. a single line text file, 2000 2000 .8 .6 1
###

def prepare_truvari_batch(samples_file, parameters_file, output_folder, bed_file):
    with open(samples_file, "r") as f:
        samples = [line.strip().split() for line in f.readlines()]

    with open(parameters_file, "r") as f:
        parameters = f.read().strip().split()

    truvari_commands_file = "truvari_commands.txt"
    with open(truvari_commands_file, "w") as truvari_file:
        for sample_row in samples:
            vcf1, vcf2, replicate_group = sample_row
            sample_name = "{}_vs_{}".format(os.path.basename(vcf1).split(".")[0], os.path.basename(vcf2).split(".")[0])
            output_subfolder = os.path.join(output_folder, "{}_{}_vs_{}".format(replicate_group, sample_name, os.path.basename(vcf2).split(".")[0]))
            output_folder_short = "_".join(["RD" + parameters[0], "LEN" + parameters[1], "OL" + parameters[2], "SEQ" + parameters[3]])
            output_name = "{}_{}_{}".format(replicate_group, sample_name, output_folder_short)

            truvari_command = ["truvari", "bench",
                            "-b", vcf1,
                            "-c", vcf2,
                            "-o", output_subfolder,
                            "--multimatch", "--passonly",
                            "--refdist", parameters[0], "--pctsize", parameters[1],
                            "--pctovl", parameters[2], "--pctsim", parameters[3],
                            "--includebed", bed_file, "-f", "/lgg_dev/lgg_data/PIPELINE/genomes/hg38.fa"]

            truvari_file.write(" ".join(truvari_command) + "\n")

            try:
                process = subprocess.Popen(truvari_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                if process.returncode != 0:
                    print("Error running Truvari for sample pair:", sample_row)
                    print("Error message:", stderr.decode("utf-8").strip())
                else:
                    # Rename the output folder to include parameter values
                    os.rename(output_subfolder, os.path.join(output_folder, output_name))
            except Exception as e:
                print("Error running Truvari for sample pair:", sample_row)
                print("Error message:", str(e))

def parse_truvari_output(folder):
    stats_file = os.path.join(folder, "summary.txt")
    if os.path.exists(stats_file):
        with open(stats_file, "r") as f:
            stats_data = json.load(f)
            stats = {
                'Precision': stats_data.get('precision'),
                'Recall': stats_data.get('recall'),
                'F1': stats_data.get('f1'),
                'Num_TP': stats_data.get('TP-call'),
                'Num_FP': stats_data.get('FP'),
                'Num_FN': stats_data.get('FN')
            }
            return stats
    else:
        return None
        
def aggregate_truvari_output(output_folder):
    output_file = "aggregate_stats.csv"
    with open(output_file, "w") as csvfile:
        fieldnames = ['Sample', 'Precision', 'Recall', 'F1', 'Num_TP', 'Num_FP', 'Num_FN']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Parse Truvari output for each sample
        for sample_folder in glob.glob(os.path.join(output_folder, "*")):
            if os.path.isdir(sample_folder):
                sample_name = os.path.basename(sample_folder)
                truvari_stats = parse_truvari_output(sample_folder)
                if truvari_stats:
                    row_data = {'Sample': sample_name}
                    for key in truvari_stats:
                        if key in ['Precision', 'Recall', 'F1']:
                            row_data[key] = "{:.2f}%".format(float(truvari_stats[key]) * 100)
                        else:
                            row_data[key] = truvari_stats[key]
                    writer.writerow(row_data)

        
def aggregate_truvari_output_old(output_folder):
    output_file = os.path.join(output_folder, "aggregate_stats.csv")
    output_file = "aggregate_stats.csv"
    with open(output_file, "w") as csvfile:
        fieldnames = ['Sample', 'Precision', 'Recall', 'F1', 'Num_TP', 'Num_FP', 'Num_FN']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # Parse Truvari output for each sample
        for sample_folder in glob.glob(os.path.join(output_folder, "*")):
            if os.path.isdir(sample_folder):
                sample_name = os.path.basename(sample_folder)
                truvari_stats = parse_truvari_output(sample_folder)
                if truvari_stats:
                    writer.writerow({'Sample': sample_name,
                                     'Precision': truvari_stats.get('Precision'),
                                     'Recall': truvari_stats.get('Recall'),
                                     'F1': truvari_stats.get('F1_Score'),
                                     'Num_TP': truvari_stats.get('num_TP'),
                                     'Num_FP': truvari_stats.get('num_FP'),
                                     'Num_FN': truvari_stats.get('num_FN')})

def main():
    samples_file = "raw_samples.txt"
    parameters_file = "parameters.txt"
    output_folder = "[output_folder_here]"
    bed_file = 'GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark.bed'  #2.79 Gb benchmark bed region
    #T2T Project bed file https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.015-20240215/

    prepare_truvari_batch(samples_file, parameters_file, output_folder, bed_file)
    aggregate_truvari_output(output_folder)

if __name__ == "__main__":
    main()
