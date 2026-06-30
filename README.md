# Supplementary preprocessing code

This repository contains preprocessing scripts used to generate the main analysis objects for the associated manuscript. The scripts cover preprocessing steps for single-cell RNA-seq, lipidomics and spatial transcriptomics analyses.

The scripts are provided for peer review and reproducibility. Local file paths, access tokens, personal computer information, installation attempts and exploratory/debugging comments have been removed.

This repository provides analysis scripts rather than a standalone software package.

## Code availability

The code is available at:

<https://github.com/HuXuantao/SCI_supplementary_code>

Input data required to run the scripts are described in `docs/input_data_structure.md`. Deposited data associated with the study are available from the data repository cited in the manuscript Data availability statement.

## Repository structure

```text
SCI_supplementary_code/
├── README.md
├── LICENSE
├── .gitignore
├── renv_notes.md
├── docs/
│   └── input_data_structure.md
└── scripts/
    ├── 01_build_scRNA_seurat_object.R
    ├── 02_build_lipidomics_object.R
    └── 03_build_spatial_seurat_object.R
```

## System requirements

The scripts were tested in an R-based analysis environment on a Linux workstation. No non-standard hardware is required for the example preprocessing workflow, although full single-cell and spatial transcriptomic analyses may require a workstation with sufficient memory.

Recommended environment:

- R 4.5.1
- Seurat 5.4.0
- Cell Ranger 9.0.1 for 10x Genomics preprocessing
- Additional R packages are listed in `renv_notes.md`

The scripts should also run on macOS or Windows if the required R packages are installed and input file paths are adapted accordingly.

## Installation

1. Install R and the required R packages listed in `renv_notes.md`.
2. Clone this repository:

```bash
git clone https://github.com/HuXuantao/SCI_supplementary_code.git
cd SCI_supplementary_code
```

3. Prepare the input files following `docs/input_data_structure.md`.
4. Update the input and output paths at the beginning of each script.

Typical installation time on a standard desktop or workstation is approximately 10–30 minutes, depending on package availability and internet connection.

## Input data

The scripts expect preprocessed or formatted input files for:

1. single-cell RNA-seq count matrices and metadata;
2. lipidomics intensity tables and sample metadata;
3. spatial transcriptomics matrices, spot coordinates and image-associated metadata.

Expected input organization and file naming are described in:

```text
docs/input_data_structure.md
```

For peer review or demonstration, users may run the scripts on a small subset of deposited data or on local test data formatted according to the same structure.

## How to run

Run the scripts in numerical order when rebuilding the main preprocessing objects:

```r
source("scripts/01_build_scRNA_seurat_object.R")
source("scripts/02_build_lipidomics_object.R")
source("scripts/03_build_spatial_seurat_object.R")
```

Each script contains path variables near the beginning of the file. These should be edited to point to the local input and output directories before running.

## Expected outputs

The expected outputs include:

- processed Seurat objects for single-cell RNA-seq analysis;
- lipidomics analysis objects and normalized data tables;
- processed Seurat objects for spatial transcriptomics analysis;
- intermediate metadata tables used for downstream analyses and figure generation.

Output file names and directories can be adjusted by editing the path variables in each script.

## Expected run time

Run time depends on dataset size and hardware. For small demonstration or subsetted data, each script is expected to run within approximately 10–60 minutes on a standard desktop or workstation. Running the full single-cell and spatial transcriptomic preprocessing workflows may take longer and may require substantially more memory.

## Reproducing preprocessing analyses

This repository contains preprocessing scripts used to generate the main analysis objects. Downstream statistical analyses, visualization and figure assembly are described in the Methods section of the manuscript and in the code comments where relevant.

To reproduce the preprocessing stage:

1. download or prepare the input files described in `docs/input_data_structure.md`;
2. install the required software and R packages;
3. update local file paths in each script;
4. run the scripts in numerical order;
5. use the generated objects for the downstream analyses described in the manuscript.

## Using the scripts on other data

Users can apply the scripts to their own single-cell RNA-seq, lipidomics or spatial transcriptomics datasets by formatting their input files according to `docs/input_data_structure.md` and updating the input/output paths in the relevant scripts.

## License

The custom analysis scripts are released under the MIT License. See `LICENSE` for details.

## Contact

For questions about the code or deposited data, please contact the corresponding author listed in the manuscript.
