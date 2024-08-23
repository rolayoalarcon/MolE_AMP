import argparse
import pandas as pd
from rdkit import Chem


def read_arguments():
    """
    This function returns parsed command line arguments.
    """

    # Instantiate parser
    parser = argparse.ArgumentParser(prog="AA sequence to SMILE conversion.",
                                     description="This program recieves a fasta file as input and resturns the SMILE of each amino acid sequence.",
                                     usage="python fasta2smile.py [options]")
    
    parser.add_argument("-f", "--filepath", help="Complete path to input FASTA file.",
                        metavar="")
    
    parser.add_argument("-o", "--outpath", help="Complete path for output file")

    # Parse arguments
    args = parser.parse_args()

    return args

def read_fasta(infile):

    """
    This function reads and pre-processes fasta file sequences.
    Amino acid (AA) sequences are stored in a dictionary {header: sequence}, where all "X" amino acids are removed
    from sequence.
    """

    # Read all lines in file
    with open(infile, "r") as file:
        fasta_content = file.read()
    
    # Create list where each element is header\nsequence
    fasta_proteins = fasta_content.split(">")[1:] # First is empty

    # Separate header and sequence
    fasta_dict = {s.split("\n", 1)[0]: s.split("\n", 1)[1].replace("\n", "").replace("X", "") for s in fasta_proteins}

    return fasta_dict

def create_mol_objs_from_aa(fasta_dict):

    # Create Mol objects from fasta sequences
    molobj_dict = {k:Chem.MolFromSequence(fasta_dict[k]) for k in fasta_dict.keys()}

    return molobj_dict

def smiles_from_mol_objs(molobj_dict):

    # Create SMILES from Mol objects
    isomeric_smiles_dict = {k:Chem.MolToSmiles(molobj_dict[k], isomericSmiles=True) for k in molobj_dict.keys()}

    # Create canonical smiles object
    canonical_smiles_dict = {k:Chem.MolToSmiles(molobj_dict[k], isomericSmiles=True) for k in molobj_dict.keys()}

    return isomeric_smiles_dict, canonical_smiles_dict


def create_combined_df(fasta_dict, isomeric_dict, canonical_dict):

    # Create a dataframe with fasta sequences
    fasta_df = pd.DataFrame(fasta_dict, index=["aa_sequence"]).transpose()

    # Create a dataframe with isomeric smiles
    isomeric_df = pd.DataFrame(isomeric_dict, index=["isomeric_smiles"]).transpose()

    # Create a dataframe with canonical smiles
    canonical_df = pd.DataFrame(canonical_dict, index=["canonical_smiles"]).transpose()

    # Join all data frames
    combined_df = fasta_df.join(isomeric_df).join(canonical_df)

    return combined_df


def aa_to_smiles(fasta_file):

    # Parse FASTA file
    fasta_parsed = read_fasta(fasta_file)

    # Create mol objects from fasta sequences
    molobjs_aa = create_mol_objs_from_aa(fasta_parsed)

    # Gather SMILES from molobjects
    iso_smiles, can_smiles = smiles_from_mol_objs(molobjs_aa)

    # Create combined df
    total_df = create_combined_df(fasta_parsed, iso_smiles, can_smiles)


    return total_df


def main():

    # Read arguments
    args = read_arguments()

    # Read fasta file as dictionary
    smiles_df = aa_to_smiles(args.filepath)

    # Write dataframe to file
    smiles_df.to_csv(args.outpath, sep='\t')

if __name__ == "__main__":
    main()









