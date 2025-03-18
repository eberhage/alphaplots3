__version_info__ = (1, 0, 0)
__version__ = '.'.join(map(str, __version_info__))
__author__ = 'Jan Eberhage, Institute for Biophysical Chemistry, Hannover Medical School (eberhage.jan@mh-hannover.de)'

import os
import re
import sys
try:
    import matplotlib
except ModuleNotFoundError:
    print('Module »matplotlib« is not installed.')
    sys.exit('Please try »python3 -m pip install matplotlib«.')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
import numpy as np
import json
import logging

def load_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def get_score(directory, score):
    """Fetches the ranking_score from summary_confidences.json in the given directory."""
    json_path = os.path.join(directory, "summary_confidences.json")
    try:
        with open(json_path, "r") as f:
            data = json.load(f)
        return data.get(score, float("-inf"))  # Default to -inf if key is missing
    except (FileNotFoundError, json.JSONDecodeError):
        log.error(f"Summary data file not found in »{directory}«. Please use a different sorting method.")
        sys.exit(1)

def plot_atom_plddts(data, output_file):
    atom_chain_ids = data["atom_chain_ids"]
    atom_plddts = data["atom_plddts"]
    
    indices = np.arange(len(atom_chain_ids))
    
    plt.figure(figsize=(10, 4))
    plt.plot(indices, atom_plddts, label='pLDDT Scores')
    
    # Add vertical lines when the chain ID changes
    prev_chain = atom_chain_ids[0]
    for i, chain in enumerate(atom_chain_ids):
        if chain != prev_chain:
            plt.axvline(i, color='black', linestyle='--', linewidth=0.8)
            prev_chain = chain
    
    plt.xlabel("Atom Index")
    plt.ylabel("pLDDT Score")
    plt.title("pLDDT Scores Along the Sequence Atoms")
    plt.legend()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_pae_heatmap(data, output_file):
    pae = np.array(data["pae"])
    token_chain_ids = data["token_chain_ids"]
    token_res_ids = data["token_res_ids"]
    
    # Identify unique chains and their corresponding indices
    unique_chains = sorted(set(token_chain_ids))
    chain_indices = {chain: [] for chain in unique_chains}
    
    for res, chain in zip(token_res_ids, token_chain_ids):
        chain_indices[chain].append(res)
    
    # Compute division positions
    division_positions = [len(chain_indices[chain]) for chain in unique_chains]
    split_indices = np.cumsum([0] + division_positions)
    num_subplots = len(unique_chains)
    
    # Create a GridSpec layout
    fig = plt.figure(figsize=(5, 5))
    gs = gridspec.GridSpec(num_subplots, num_subplots, figure=fig, 
                           width_ratios=division_positions,
                           height_ratios=division_positions,
                           wspace=0.05, hspace=0.05)
        
    # Plot each submatrix
    for i in range(num_subplots):
        for j in range(num_subplots):
            row_start, row_end = split_indices[i], split_indices[i+1]
            col_start, col_end = split_indices[j], split_indices[j+1]
            submatrix = pae[row_start:row_end, col_start:col_end]
            
            ax = fig.add_subplot(gs[i, j])
            im = ax.imshow(submatrix, cmap="bwr", origin='upper', vmin=0, vmax=30)
            
            # Add black border around each subplot
            for spine in ax.spines.values():
                spine.set_edgecolor('black')  
            
            # Show y-axis labels only for the first column
            if j == 0:
                ax.set_ylabel(f"Chain {unique_chains[i]}")
            else:
                ax.set_yticks([])
            
            # Show x-axis labels only for the last row
            if i == num_subplots - 1:
                ax.set_xlabel(f"Chain {unique_chains[j]}")
            else:
                ax.set_xticks([])
    
    cax = fig.add_axes([0.92, 0.11, 0.025, 0.77])  # [left, bottom, width, height]
    clb = fig.colorbar(im, cax=cax)
    clb.ax.set_title('Å')
    
    plt.savefig(output_file, dpi=300)
    plt.close()


def process_input_dir(input_dir, nameprefix, sortingmethod):
    name_plddt = nameprefix + ('_' if nameprefix else '') + "plddt.png"
    name_pae = nameprefix + ('_' if nameprefix else '') + "pae.png"
    name_master_pae = nameprefix + ('_' if nameprefix else '') + "master_pae.png"
    pattern = re.compile(r"seed-\d+_sample-\d+")
    parent_directories = {}
    
    for root, dirs, _ in os.walk(input_dir):
        for subdir in dirs:
            if pattern.match(subdir):
                parent_dir = os.path.abspath(root)
                if parent_dir not in parent_directories:
                    parent_directories[parent_dir] = []
                subdir_path = os.path.join(root, subdir)
                parent_directories[parent_dir].append(subdir_path)
                data_file = os.path.join(subdir_path, "confidences.json")
                
                if os.path.exists(data_file):
                    data = load_json(data_file)
                    plot_atom_plddts(data, os.path.join(subdir_path, name_plddt))
                    plot_pae_heatmap(data, os.path.join(subdir_path, name_pae))
                    log.info(f"Generating plots in »{subdir_path}«.")
                else:
                    log.error(f"Data file not found in »{subdir_path}«.")                
    
    for parent, subdirs in parent_directories.items():
        num_images = len(subdirs)
        horizontal = np.ceil(np.sqrt(num_images)).astype(int)
        vertical = np.ceil(num_images / horizontal).astype(int)
        _, axes = plt.subplots(vertical, horizontal, figsize=(5 * horizontal, 5 * vertical), gridspec_kw = {'wspace':0.01, 'hspace':0.03})
        axes = np.array(axes).reshape(vertical, horizontal)  # Ensure it's always a 2D array
        
        scores = {"rank": ("ranking", "ranking_score"), "iptm": ("iPTM", "iptm"), "ptm": ("PTM", "ptm")} # pairs of Label and Key
        if sortingmethod in scores:
            ordered_subdirs = sorted(subdirs, key=lambda x: get_score(x, scores[sortingmethod][1]), reverse=True)
        else:
            ordered_subdirs = sorted(subdirs)  # sort by seed and by sample (alphabetically)

        for ax, subdir in zip(axes.flatten(), ordered_subdirs):
            img_path = os.path.join(subdir, name_pae)
            img = plt.imread(img_path)
            ax.imshow(img)
            ax.axis("off")
            titletext = os.path.basename(subdir)
            if sortingmethod in scores:
                label, score_key = scores[sortingmethod]
                titletext += f" ({label}: {get_score(subdir, score_key)})"
            ax.text(0.5, 0.93, titletext, ha='center', va='top', fontweight="bold", transform=ax.transAxes)

        # Hide unused subplots
        for ax in axes.flatten()[len(ordered_subdirs):]:
            ax.axis("off")
        
        master_output_file = os.path.join(parent, name_master_pae)
        plt.savefig(master_output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Generated master plot »{master_output_file}« with sorting method »{sortingmethod}«.")


class CustomFormatter(logging.Formatter):
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    #format = '%(asctime)s [%(levelname)-7s] %(message)s'
    format = '%(asctime)s :: %(message)s'
    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)

def main():
    global log
    log = logging.getLogger(__name__)
    log.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(CustomFormatter())
    log.addHandler(ch)

    parser = argparse.ArgumentParser(add_help=False, description='This script will generate plots containing the pLDDT \
      distribution and Predicted Alignment Error (PAE) of a given AlphaFold output using the JSON-files.')
    parser.add_argument('input_dir', nargs="?", default=".", metavar='<input_dir>',
                          help='Relative or absolute path to the input directory (AlphaFold output)')
    parser.add_argument('-n', '--name', dest='name', default='', metavar='<prefix>',
                          help='Add custom name as prefix to your plots')
    parser.add_argument('-s', '--sort', dest='sort', default='seedsample', metavar='<sorting method>',
                          help='Sorting method for the master plots. Default: seedsample. Other options: iptm, ptm, rank.')
    parser.add_argument('-v', '--version', action='version',
                          version='%(prog)s (' + __version__ + ') ' + ' by ' + __author__)
    parser.add_argument('-h', '--help', action='help', help='show this help message and exit')
    args = parser.parse_args()

    if not os.path.exists(args.input_dir):
        log.error(f'»{os.path.abspath(args.input_dir)}« was not found. Aborting!')
        sys.exit(1)

    process_input_dir(args.input_dir, args.name, args.sort)

if __name__ == "__main__":
    main()
