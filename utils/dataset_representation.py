import os
import yaml
import numpy as np
import pandas as pd

import torch
from torch_geometric.data import Data, Dataset, Batch

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem.rdchem import BondType as BT


RDLogger.DisableLog('rdApp.*')


# Global variables

# Element list
ATOM_LIST = list(range(1,119))
# Chirality
CHIRALITY_LIST = [
    Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW,
    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW,
    Chem.rdchem.ChiralType.CHI_OTHER
]
BOND_LIST = [
    BT.SINGLE, 
    BT.DOUBLE, 
    BT.TRIPLE, 
    BT.AROMATIC
]
BONDDIR_LIST = [
    Chem.rdchem.BondDir.NONE,
    Chem.rdchem.BondDir.ENDUPRIGHT,
    Chem.rdchem.BondDir.ENDDOWNRIGHT
]

def read_fasta(filepath):

    """
    This function reads and pre-processes fasta file sequences.
    
    Parameters:
    - filepath (str): Complete path to a fasta file of Amimo Acids. 

    Returns:
    - fasta_dict (dictionary): Amino acid (AA) sequences are stored in a dictionary {header: sequence}, where all "X" amino acids are removed from sequence.

    """

    # Read all lines in file
    with open(filepath, "r") as file:
        fasta_content = file.read()
    
    # Create list where each element is header\nsequence
    fasta_proteins = fasta_content.split(">")[1:] # First is empty

    # Separate header and sequence
    fasta_dict = {s.split("\n", 1)[0]: s.split("\n", 1)[1].replace("\n", "").replace("X", "") for s in fasta_proteins}

    return fasta_dict


class MoleculeDataset(Dataset):

    """
    Dataset class for creating molecular graphs.

    Attributes:
    - protein_dict (dict): Dictionary containing {ID: sequence} information.
    """

    def __init__(self, protein_dict):
        super(Dataset, self).__init__()

        # Gather the AA sequences and corresponding IDs
        self.aa_data = list(protein_dict.values())
        self.id_data = list(protein_dict.keys())

    def __getitem__(self, index):
        # Get the molecule object from the amino acid sequence.
        mol = Chem.MolFromSequence(self.aa_data[index])
        mol = Chem.AddHs(mol)

        #########################
        # Get the molecule info #
        #########################
        type_idx = []
        chirality_idx = []
        atomic_number = []

        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                print(self.id_data[index])

            type_idx.append(ATOM_LIST.index(atom.GetAtomicNum()))
            chirality_idx.append(CHIRALITY_LIST.index(atom.GetChiralTag()))
            atomic_number.append(atom.GetAtomicNum())

        x1 = torch.tensor(type_idx, dtype=torch.long).view(-1,1)
        x2 = torch.tensor(chirality_idx, dtype=torch.long).view(-1,1)
        x = torch.cat([x1, x2], dim=-1)

        row, col, edge_feat = [], [], []
        for bond in mol.GetBonds():
            start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            row += [start, end]
            col += [end, start]
            edge_feat.append([
                BOND_LIST.index(bond.GetBondType()),
                BONDDIR_LIST.index(bond.GetBondDir())
            ])
            edge_feat.append([
                BOND_LIST.index(bond.GetBondType()),
                BONDDIR_LIST.index(bond.GetBondDir())
            ])

        edge_index = torch.tensor([row, col], dtype=torch.long)
        edge_attr = torch.tensor(np.array(edge_feat), dtype=torch.long)

        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, 
                    chem_id=self.id_data[index])
        
        return data

    def __len__(self):
        return len(self.aa_data)
    
    def get(self, index):
        return self.__getitem__(index)

    def len(self):
        return self.__len__()

# Function to generate the molecular representation with MolE
def batch_representation(fasta_dictionary, dl_model, batch_size= 10_000, device="cuda:0"):

    """
    Generate molecular representations using a Deep Learning model.

    Parameters:
    - fasta_dictionary (dict): Dictionary with {ID: sequence} information.
    - dl_model: Deep Learning model for molecular representation.
    - batch_size (int, optional): Batch size for processing (default is 10,000).
    - device (str, optional): Device for computation (default is "cuda:0").

    Returns:
    - chem_representation (pandas.DataFrame): DataFrame containing molecular representations.
    """
    
    # First we create a list of graphs
    molecular_graph_dataset = MoleculeDataset(fasta_dictionary)
    graph_list = [g for g in molecular_graph_dataset]

    # Determine number of loops to do given the batch size
    n_batches = len(graph_list) // batch_size

    # Are all molecules accounted for?
    remaining_molecules = len(graph_list) % batch_size

    # Starting indices
    start, end = 0, batch_size

    # Determine number of iterations
    if remaining_molecules == 0:
        n_iter = n_batches
    
    elif remaining_molecules > 0:
        n_iter = n_batches + 1
    
    # A list to store the batch dataframes
    batch_dataframes = []

    # Iterate over the batches
    for i in range(n_iter):
        # Start batch object
        batch_obj = Batch()
        graph_batch = batch_obj.from_data_list(graph_list[start:end])
        graph_batch = graph_batch.to(device)

        # Gather the representation
        with torch.no_grad():
            dl_model.eval()
            h_representation, _ = dl_model(graph_batch)
            chem_ids = graph_batch.chem_id
        
        batch_df = pd.DataFrame(h_representation.cpu().numpy(), index=chem_ids)
        batch_dataframes.append(batch_df)

        # Get the next batch
        ## In the final iteration we want to get all the remaining molecules
        if i == n_iter - 2:
            start = end
            end = len(graph_list)
        else:
            start = end
            end = end + batch_size
    
    # Concatenate the dataframes
    chem_representation = pd.concat(batch_dataframes)

    return chem_representation

# Function to load a pre-trained model
def load_pretrained_model(pretrained_model, device="cuda:0"):

    """
    Load a pre-trained MolE model.

    Parameters:
    - pretrain_architecture (str): Architecture of the pre-trained model.
    - pretrained_model (str): The path to the directory of the pre-trained MolE model.
    - device (str, optional): Device for computation (default is "cuda:0").

    Returns:
    - model: Loaded pre-trained model.
    """

    # Read model configuration
    config = yaml.load(open(os.path.join(pretrained_model, "config.yaml"), "r"), Loader=yaml.FullLoader)
    model_config = config["model"]

    # Instantiate model
    from utils.ginet_concat import GINet
    model = GINet(**model_config).to(device)
    
    # Load pre-trained weights
    model_pth_path = os.path.join(pretrained_model, "model.pth")
    print(model_pth_path)

    state_dict = torch.load(model_pth_path, map_location=device)
    model.load_my_state_dict(state_dict)

    return model

def process_dataset(fasta_filepath, pretrained_model, device="cuda:0"):
    """
    Process the dataset to generate molecular representations.

    Parameters:
    - fasta_filepath (str): Path to the FASTA file.
    - pretrained_model (str): Name of the pre-trained model. Can also be "MolCLR" or "ECFP4".
    - device (str): Device to use for computation (default is "cuda:0"). Can also be "cpu".

    Returns:
    - udl_representation (pandas.DataFrame): DataFrame containing molecular representations.
    """

    # Read the FASTA file into a dictionary. 
    fasta_dict = read_fasta(fasta_filepath)

    # Now we load our pretrained model
    pmodel = load_pretrained_model(pretrained_model, device=device)
    
    # Obtain the requested representation
    udl_representation = batch_representation(fasta_dict, pmodel, device=device)

    return udl_representation



    
