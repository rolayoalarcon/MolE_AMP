# Using MolE for antimicrobial peptides.  
  
In this repository, we use the MolE pre-trained representation of molecular structures to make predictions of the antimicrobial acitivity of peptides. The MolE and XGBoost models are the same that are used in our pre-print:
  
- [Pre-trained molecular representations enable antimicrobial discovery](https://www.biorxiv.org/content/10.1101/2024.03.11.584456v2).

For predictions on small molecules look at our [mole_antimicrobial_potential](https://github.com/rolayoalarcon/mole_antimicrobial_potential) repository. Pre-training details can be found in our manuscript and the [MolE](https://github.com/rolayoalarcon/MolE) GitHub repository.

 
The rest of this tutorial goes over how to install and use this repository.

## Installation  
  
Firstly, you can create a conda environment with all the necessary dependencies.

```
# Create conda environment from yaml file
$ conda env create -f environment.yml

# Afterwards activate the environment
$ conda activate mole_amp
```

Next, make sure to download the MolE pre-trained model from [here](https://zenodo.org/records/10803099?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImI3NTg0OTU0LTI5YWItNDgxZS04OGYyLTU5MmM1MjcwYzJjZiIsImRhdGEiOnt9LCJyYW5kb20iOiIzNzgyNTE5ZGU5N2MzZWI3YjZiZjkwYTIzZjFiMmEwZSJ9.oL6G0WZKxIowSb-2qdP55cPhef1W4yG5iF4PFlsWPpuPROmzRhutJtySzs9q02ACltl0qy9YPJjzB7NvzRMyaw) and place the `model.pth` file in the mole_pretrained/model_ginconcat_btwin_100k_d8000_l0.0001/ subdirectory. The hyperparameters of the pre-trained model used in this work are specified in `pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001/config.yaml`.  
  
Once these two steps are completed, you are ready to get started!  
  
## Predictions of antimicrobial activity of peptides.  

### Input:  
  
The input consists of of a fasta file of amino acid sequences. An example can be seen in `example/input/test_fasta.faa`

```
>WP_002558051.1
MKVRASLKKRTPECKIVRRNGRLYVINKKNPKYKQRQG
>TFASTA.1
IPXQ
```
During pre-processing all **X** amino acids are removed.

Additional inputs required are:
 1. An XGBoost Classifier model that has been pickled (see `xgboost_models/MolE-XGBoost-08.03.2024_14.20.pkl`)
 2. Pre-trained MolE model consisting of `config.yaml` and `model.pth` (see `pretrained_model/model_ginconcat_btwin_100k_d8000_l0.0001/`)

Also, the files found in `maier_information` are also required to make strain-level and aggregated predictions. See further for details.
  
### Making predictions:  

In order to make predictions one can make use of the `inference_script.py` provided in this repository 

```{code}
$ python inference_script.py  -h
usage: python inference_script.py fasta_filepath outpath [options]

This program recieves a fasta file as input, featurizes the molecules using MolE, then makes predictions of antimicrobial activity

positional arguments:
  fasta_filepath        Complete path to input FASTA file.
  outpath               Complete path for output file

optional arguments:
  -h, --help            show this help message and exit
  -a, --aggregate_scores
                        Flag variable. If called, then prediction scores are aggregated by compound using as the antimicrobial potential of each compound.
  -x XGBOOST_MODEL, --xgboost_model XGBOOST_MODEL
                        Path to the pickled XGBoost model that makes predictions (.pkl). Default set to: xgboost_models/MolE-XGBoost-08.03.2024_14.20.pkl
  -m MOLE_MODEL, --mole_model MOLE_MODEL
                        Path to the directory containing the config.yaml and model.pth files of the pre-trained MolE chemical representation. Default set to:
                        mole_pretrained/model_ginconcat_btwin_100k_d8000_l0.0001
  -s STRAIN_CATEGORIES, --strain_categories STRAIN_CATEGORIES
                        Path to the Maier et.al. screening results. Default is set to ./maier_information/maier_screening_results.tsv.gz
  -g GRAM_INFORMATION, --gram_information GRAM_INFORMATION
                        Path to strain metadata. Default is set to ./maier_information/strain_info_SF2.xlsx
  -t APP_THRESHOLD, --app_threshold APP_THRESHOLD
                        threshold score to binarize compound-microbe predictions. Default from original publication.
  -k MIN_NKILL, --min_nkill MIN_NKILL
                        Minimum number of microbes predicted to be inhibited in order to consider broad spectrum antibiotic.
  -d DEVICE, --device DEVICE
                        Device where the pre-trained model is loaded. Can be one of ['cpu', 'cuda', 'auto']. If 'auto' (default) then cuda:0 device is selected if a GPU is detected.
```

The models provided are set as deaults. Default threshold values are those used in our manuscript. Therefore you only need to indicate the input fasta file and the output file names. 


## Example

In order to make strain-level predictions of antimicrobial activity, execute the following command.  


```
$ python inference_script.py example/input/test_fasta.faa example/output/prediction_strains.tsv

```

The output can be seen in `example/output/prediction_strains.tsv`.  
  

| pred_id | 0 | 1 | growth_inhibition |
| ------- | - | - | ----------------- |
| WP_002558051.1:Akkermansia muciniphila | 0.93703735 | 0.06296267 | 1 |
| WP_002558051.1:Escherichia coli ED1a (NT5078) | 0.98663974 | 0.013360276 | 1 |

- **pred_id**: Indicates the protein - microbe combination for which the prediction is made.
- **0**: Indicates the prediction score for the protein _not_ having activity against the microbe.
- **1**: Indicates the prediction score fot the protein _having_ antimicrobial activity against the microbe.
- **growth_inhibition**: A binary column indicating whether the compound is predicted to inhibit the microbe's growth (1) or not (0). This is acheived by threshold the values of column **1** using a pre-determined score threshold. By default, the same threshold used in our publication is applied. This might not be appropriate for AMPs.
  
  
Predictions can be aggregated to get **Antimicrobial Potential scores** for each compound by using the `-a` flag.  
The **Antimicrobial Potential** is described in our [pre-print](https://www.biorxiv.org/content/10.1101/2024.03.11.584456v2). Briefly, it is the geometric mean of the predictions scores for antimicrobial activity across all species. In this way, the Antimicrobial Potential indicates a probability of the compound having a broad-spectrum of activity.
Additionally, the predicted **number of inhibited strains** can also by 
  
```
$ python inference_script.py example/input/test_fasta.faa example/output/prediction_aggregated.tsv -a  
```  
  
The output is saved to `example/output/prediction_aggregated.tsv`.  
  
| chem_id | apscore_total | apscore_gnegative | apscore_gpositive | ginhib_total | ginhib_gnegative | ginhib_gpositive | broad_spectrum |
| ------- | ------------- | ----------------- | ----------------- | ------------ | ---------------- | ---------------- | -------------- |
| DBAASPR_11 | 0.17025694 | 0.10979764 | 0.24376883 | 36 | 14 |	22 |	1 |
| DBAASPR_14 | 0.27378127 | 0.1664661 | 0.41133404 | 36 | 15 | 21 |	1 |

- **chem_id**: Is identifier of each protein.
- **apscore_total**: Is the aggregated Antimicrobial Potential score of each protein.
- **apscore_gnegative**: Is the aggregated Antimicrobial Potential score of each protein on the subset of Gram Negative species.
- **apscore_gpositive**: Is the aggregated Antimicrobial Potential score of each protein on the subset of Gram Positive species.
- **ginhib_total**: Is the total amount strains that are predicted to be inhibited by each protein.
- **ginhib_gpositive**: Is the total amount Gram Positive strains that are predicted to be inhibited by each protein.
- **ginhib_gnegative**: Is the total amount Gram Negative strains that are predicted to be inhibited by each protein.
- **broad_spectrum**: Is a binary column indicating whether the protein is predicted to be a broad spectrum antimicrobial (1) or not (0). This is done by thresholding the total amount of strains predicted to be inhibited (**ginhib_total**). By default if a protein is predicted to inhibit $\geq$ 10 strains, then it is predicted to have broad spectrum activity.