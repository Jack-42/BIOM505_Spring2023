import subprocess
import argparse
import os
from pathlib import Path
import csv
import pandas as pd


def create_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-prof_dir",
        "--prof_dir",
        help="Path to input directory containing .prof files",
        type=Path,
        required=True, )

    parser.add_argument(
        "-ddg_csv",
        "--ddg_csv",
        help="Path to csv file containing ddG values and chain info",
        type=Path,
        required=True, )

    parser.add_argument(
        "-o",
        "--out_csv_path",
        help="Output file for csv file",
        type=Path,
        required=True, )

    return parser


def run(args):
    prof_dir = args.prof_dir
    ddg_csv = args.ddg_csv
    assert os.path.exists(
        prof_dir), "Given profile dir not found: %s" % prof_dir
    assert os.path.exists(ddg_csv), "Given ddg csv not found: %s" % ddg_csv

    prot_df = pd.read_csv(ddg_csv)

    out_csv_path = args.out_csv_path
    os.makedirs(os.path.dirname(out_csv_path), exist_ok=True)

    prof_files = [filename for filename in os.listdir(prof_dir) if
                  filename.endswith('.prof')]
    with open(out_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['id', 'mut_code', 'experimental_ddg',
                         'predicted_ddg'])  # Write the header row
        for prof_file in prof_files:
            prof_path = os.path.join(prof_dir, prof_file)
            uni_id = os.path.splitext(prof_file)[0]
            sub_df = prot_df.loc[prot_df['uniprot_id'] == uni_id]
            mut_codes = sub_df.apply(lambda x: "%s%d%s" % (
                x['wild_type'], x['position'], x['mutation']), axis=1)
            ddgs = sub_df['ddG']
            for mut_code, ddg in zip(mut_codes, ddgs):
                command = "acdc-nn seq %s %s" % (mut_code, prof_path)
                try:
                    result = subprocess.run(command, stdout=subprocess.PIPE,
                                            shell=True, check=True)
                except:
                    print(command)
                    continue
                output_ddg = result.stdout.decode('utf-8').strip()
                output_ddg = float(output_ddg)
                writer.writerow([uni_id, mut_code, ddg, output_ddg])


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
