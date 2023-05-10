import subprocess
import argparse
import os
from pathlib import Path
import csv
import pandas as pd


def create_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-wild_prof_dir",
        "--wild_prof_dir",
        help="Path to input directory containing .prof files for wildtypes",
        type=Path,
        required=True, )

    parser.add_argument(
        "-mut_prof_dir",
        "--mut_prof_dir",
        help="Path to input directory containing .prof files for mutants",
        type=Path,
        required=True, )

    parser.add_argument(
        "-wild_pdb_dir",
        "--wild_pdb_dir",
        help="Path to input directory containing wildtype pdbs",
        type=Path,
        required=True, )
    
    parser.add_argument(
        "-mut_pdb_dir",
        "--mut_pdb_dir",
        help="Path to input directory containing mutated pdbs",
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
        help="Output path to csv file",
        type=Path,
        required=True, )

    return parser


def get_inverse_substitution(substitution: str):
    """
    Get the inverse of a given substitution
    :param str substitution: mutation code, expected format {wild residue}{position}{mutated residue}, e.g. A43R
    :return: str, inverse substitution (e.g., R43A)
    """
    wild_residue = substitution[0]
    mut_residue = substitution[-1]
    position = substitution[1:-1]
    return mut_residue + position + wild_residue


def run(args):
    wild_prof_dir = args.wild_prof_dir
    wild_pdb_dir = args.wild_pdb_dir
    mut_prof_dir = args.mut_prof_dir
    mut_pdb_dir = args.mut_pdb_dir
    ddg_csv = args.ddg_csv
    assert os.path.exists(
        wild_prof_dir), "Given wild_prof_dir not found: %s" % wild_prof_dir
    assert os.path.exists(wild_pdb_dir), "Given wild_pdb_dir not found: %s" % wild_pdb_dir
    assert os.path.exists(
        mut_prof_dir), "Given mut_prof_dir not found: %s" % mut_prof_dir
    assert os.path.exists(mut_pdb_dir), "Given mut_pdb_dir not found: %s" % mut_pdb_dir
    assert os.path.exists(ddg_csv), "Given ddg_csv not found: %s" % ddg_csv

    prot_df = pd.read_csv(ddg_csv)

    out_csv_path = args.out_csv_path
    os.makedirs(os.path.dirname(out_csv_path), exist_ok=True)

    wild_prof_files = [filename for filename in os.listdir(wild_prof_dir) if
                  filename.endswith('.prof')]
    wild_pdb_files = [filename for filename in os.listdir(wild_pdb_dir) if
                 filename.endswith('.pdb')]

    # Filter both lists to only include entries that are shared by the other list
    prof_base = {os.path.splitext(filename)[0] for filename in wild_prof_files}
    pdb_base = {os.path.splitext(filename)[0] for filename in wild_pdb_files}
    matching_bases = prof_base.intersection(pdb_base)
    wild_prof_files = sorted([filename for filename in wild_prof_files if
                         os.path.splitext(filename)[0] in matching_bases])
    wild_pdb_files = sorted([filename for filename in wild_pdb_files if
                        os.path.splitext(filename)[0] in matching_bases])
    with open(out_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['id', 'mut_code', 'chain', 'experimental_ddg',
                         'predicted_ddg'])  # Write the header row
        for prof_file, pdb_file in zip(wild_prof_files, wild_pdb_files):
            assert (os.path.splitext(prof_file)[0] ==
                    os.path.splitext(pdb_file)[0]), "%s, %s" % (
                prof_file, pdb_file)
            wild_prof_path = os.path.join(wild_prof_dir, prof_file)
            wild_pdb_path = os.path.join(wild_pdb_dir, pdb_file)
            uni_id = os.path.splitext(prof_file)[0]
            sub_df = prot_df.loc[prot_df['uniprot_id'] == uni_id]
            chains = sub_df['chain']
            mut_codes = sub_df.apply(lambda x: "%s%d%s" % (
                x['wild_type'], x['position'], x['mutation']), axis=1)
            ddgs = sub_df['ddG']
            for mut_code, chain, ddg in zip(mut_codes, chains, ddgs):
                chain = 'A'  # all AlphaFold/ESMFold pdbs just have 'A' chain
                inv_mut = get_inverse_substitution(mut_code)
                mut_file_name = "%s_%s" % (uni_id, mut_code)
                mut_prof_path = os.path.join(mut_prof_dir, mut_file_name + ".prof")
                mut_pdb_path = os.path.join(mut_pdb_dir, mut_file_name + ".pdb")
                if not(os.path.exists(mut_prof_path)):
                    print("Mutant profile file not found:", mut_prof_path)
                    continue
                if not(os.path.exists(mut_pdb_path)):
                    print("Mutant pdb file not found:", mut_pdb_path)
                    continue
                command = "acdc-nn istruct %s %s %s %s %s %s %s %s" % (
                    mut_code, wild_prof_path, wild_pdb_path, chain, inv_mut, mut_prof_path, mut_pdb_path, chain)
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
