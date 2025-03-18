# alphaplots3.py
This script will scan a folder and use the confidence.json files of an AlphaFold3 output to generate the per-Atom plDDT distribution and a predicted alignment error (PAE) plot.

## Example Usage
Make plots. That's it. If the positional argument for the directory is empty, your current directory is assumed.
```
python3 alphaplots3.py /path/to/alphafold/result/output/directory
```

## Parameters
### 1st Positional (optional)
```
<input_dir>
```
### Optional
To add a prefix to the output file names
```
-n <prefix>
--name <prefix>
(default = None)
```
To sort the PAE plots for the master_pae.png. Options: "seedsample" (default), "rank", "iptm" or "ptm".
```
-s <sorting method>
--sort <sorting method>
(default = seedsample)
```
Version
```
-v
--version
```
Help
```
-h
--help
```
