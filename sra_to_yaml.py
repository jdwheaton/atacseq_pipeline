import argparse
import pandas as pd
import yaml

'''
This script takes a csv file, such as returned by SRA, and
converts the pertinent info to a YAML file to be used in
the config file of a snakemake workflow.
'''

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--metadata", type=str, metavar="<csv_file>",
                    nargs=1, required=True,
                    help="CSV file containing metadata")
parser.add_argument("-o", "--outfile", type=str, nargs=1, required=True)
parser.add_argument("-i", "--id_var", type=str, nargs=1, required=True,
                    help="Column name for ID variable identifying individual files (i.e. SRR Run)")
parser.add_argument("-g", "--group", type=str, nargs='+', required=True,
                    help="Column names for grouping variables")

if __name__ == "__main__":
    args = parser.parse_args()

    data = pd.read_csv(args.metadata[0]).fillna('None')
    data['condition'] = data[args.group].apply(lambda x: '-'.join(x), axis=1)
    # Remove spaces and the plus character
    data.condition = data.condition.str.replace(" ", "_").str.replace("+", "")
    # Switch to a multi-index with condition first, then Run ID (this enables aggregation)
    data_multiindex = data.sort_values("condition").set_index(
        ["condition", args.id_var[0]]).index

    d = {}
    for condition, sample in data_multiindex:
        if condition not in d.keys():
            d[condition] = []
        d[condition].append(sample)
    config = {'samples': d}

    with open(args.outfile[0], 'a') as f:
        f.write(yaml.dump(config, default_flow_style=False))
