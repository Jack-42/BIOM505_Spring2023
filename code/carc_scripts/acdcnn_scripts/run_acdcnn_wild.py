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
        "-pdb_dir",
        "--pdb_dir",
        help="Path to input directory containing wildtype pdbs",
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
    pdb_dir = args.pdb_dir
    ddg_csv = args.ddg_csv
    assert os.path.exists(
        prof_dir), "Given profile dir not found: %s" % prof_dir
    assert os.path.exists(pdb_dir), "Given pdb dir not found: %s" % pdb_dir
    assert os.path.exists(ddg_csv), "Given ddg csv not found: %s" % ddg_csv

    prot_df = pd.read_csv(ddg_csv)

    out_csv_path = args.out_csv_path
    os.makedirs(os.path.dirname(out_csv_path), exist_ok=True)

    prof_files = [filename for filename in os.listdir(prof_dir) if
                  filename.endswith('.prof')]
    pdb_files = [filename for filename in os.listdir(pdb_dir) if
                 filename.endswith('.pdb')]

    # Filter both lists to only include entries that are shared by the other list
    prof_base = {os.path.splitext(filename)[0] for filename in prof_files}
    pdb_base = {os.path.splitext(filename)[0] for filename in pdb_files}
    matching_bases = prof_base.intersection(pdb_base)
    prof_files = sorted([filename for filename in prof_files if
                         os.path.splitext(filename)[0] in matching_bases])
    pdb_files = sorted([filename for filename in pdb_files if
                        os.path.splitext(filename)[0] in matching_bases])
    with open(out_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['id', 'mut_code', 'chain', 'experimental_ddg',
                         'predicted_ddg'])  # Write the header row
        for prof_file, pdb_file in zip(prof_files, pdb_files):
            assert (os.path.splitext(prof_file)[0] ==
                    os.path.splitext(pdb_file)[0]), "%s, %s" % (
                prof_file, pdb_file)
            prof_path = os.path.join(prof_dir, prof_file)
            pdb_path = os.path.join(pdb_dir, pdb_file)
            uni_id = os.path.splitext(prof_file)[0]
            sub_df = prot_df.loc[prot_df['uniprot_id'] == uni_id]
            chains = sub_df['chain']
            mut_codes = sub_df.apply(lambda x: "%s%d%s" % (
                x['wild_type'], x['position'], x['mutation']), axis=1)
            ddgs = sub_df['ddG']
            for mut_code, chain, ddg in zip(mut_codes, chains, ddgs):
                chain = 'A'  # all AlphaFold/ESMFold pdbs just have 'A' chain
                command = "acdc-nn struct %s %s %s %s" % (
                    mut_code, prof_path, pdb_path, chain)
                try:
                    result = subprocess.run(command, stdout=subprocess.PIPE,
                                            shell=True, check=True)
                except:
                    print(command)
                    continue
                output_ddg = result.stdout.decode('utf-8').strip()
                output_ddg = float(output_ddg)
                writer.writerow([uni_id, mut_code, chain, ddg, output_ddg])


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
