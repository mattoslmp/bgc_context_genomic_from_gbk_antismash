# BGC Toolkit with Interactive Visualizations

The **BGC Toolkit** is a command-line application designed for the analysis and visualization of Biosynthetic Gene Clusters (BGCs) from GenBank files. It filters for gene clusters containing the CAL_domain (and AMP-binding domains) and generates a variety of outputs, including interactive HTML visualizations, combined SVG figures, bar charts displaying product distribution by taxonomy, FASTA files of protein sequences, and metadata files.

## Overview

This tool performs the following tasks:
- **Data Collection & Filtering:** Searches specified input folders for GenBank (`.gbk`) files, extracts BGCs that contain a CAL_domain, and collects various related information.
- **Interactive Visualizations:** Generates interactive HTML visualizations using Plotly to display the gene clusters with color-coded domains.
- **SVG Generation:** Produces combined SVG images where each SVG file can contain up to 10 BGCs with graphical annotations and legends.
- **Taxonomic Analysis:** Retrieves taxonomy data via NCBI Entrez (e.g., Phylum, Order) and creates stacked bar charts in HTML to show the distribution of candidate products across taxonomic levels.
- **Additional Outputs:** Optionally saves GenBank files, protein sequences in FASTA format, and metadata text files summarizing the BGC information.

## Features

- **Multi-folder GenBank Processing:** Automatically scans multiple directories for `.gbk` files.
- **BGC Extraction:** Filters BGCs based on the presence of CAL_domain and AMP-binding domains.
- **Color-Coded Domain Annotation:** Uses a color generator (with a specific preset color for AMP-binding proteins) to assign unique colors for each asDomain description.
- **Interactive Plotting:** Leverages Plotly for dynamic, interactive HTML visualizations.
- **SVG Visualization:** Generates scalable vector graphics (SVG) for high-quality, publication-ready images.
- **Taxonomy Enrichment:** Uses NCBI Entrez to fetch taxonomy details and integrates this data into both the visualization and statistical summaries.
- **Flexible Output Options:** Supports exporting data in various formats (HTML, SVG, FASTA, and plain text metadata).

## Requirements

- **Python Version:** Python 3.6 or higher
- **Dependencies:**
  - [Biopython](https://biopython.org)
  - [Plotly](https://plotly.com/python/)
  - [lxml](https://lxml.de/)
  - [numpy](https://numpy.org)
  - [pandas](https://pandas.pydata.org)
  - Standard Python libraries: `os`, `argparse`, `pathlib`, `shutil`, `time`, `colorsys`

You can install the required packages using pip:

```bash
pip install biopython plotly lxml numpy pandas
```

# Installation
Clone or Download: Get the source code from your repository or download the script file.

Ensure Dependencies: Verify that all the required Python packages are installed.

Configure Email: Update the Entrez.email variable in the script with your actual email address to ensure proper usage of NCBI Entrez.

# Usage
Run the script from the command line by specifying the required arguments. For example:
```bash
python bgc_toolkit.py -i /path/to/folder1 /path/to/folder2 -o ./output --html --svg --fasta --metadata my_metadata
```

# Command-Line Arguments
-i, --inputfolders
Type: One or more folder paths

# Description: Folder(s) to search for GenBank (.gbk) files.

-o, --outputfolder
Type: Path

Default: ./output
# Description: Base folder where results (SVG, HTML, FASTA, metadata, etc.) will be stored.

--svg
Type: Flag
Description: Enables SVG output for each individual BGC.

--svgcfg
Type: Path
Default: SVG_arrow_options.cfg
Description: Configuration file containing SVG style options.

--genbank
Type: Flag
Description: Saves the GenBank files into the output folder.

--metadata
Type: String

# Description: Base name for metadata files that summarize BGC information.

--fasta
Type: Flag
Description: Saves protein sequences (with the CAL_domain) in FASTA format.

--html
Type: Flag
Description: Generates an interactive HTML visualization of the BGCs.

--svgall
Type: Flag
Description: Generates combined SVG images containing up to 10 BGCs each.

# How It Works
# Data Collection:

The script scans the specified input folders for GenBank files. It parses each file using Biopython and extracts relevant data such as genome ID, species name, and BGC features (filtered by presence of CAL_domain).

# Taxonomic Annotation:

Using NCBI Entrez, the tool retrieves taxonomy details (e.g., phylum, order) for each BGC based on the organism’s genus.

# Visualization & Output Generation:

Interactive HTML: Uses Plotly to construct interactive visualizations of BGCs with detailed annotations.

SVG Images: Generates scalable SVG graphics that display gene clusters and asDomain legends.

Bar Charts: Creates stacked bar charts to show product counts by Phylum, Order, and Genus.

Export Options: Optionally saves GenBank files, FASTA files for protein sequences, and metadata summaries.

Color Mapping:
Unique colors are automatically assigned to different domain descriptions. A fixed bright red color is reserved for AMP-binding proteins, while other colors are generated in the HSL color space to ensure visually distinct assignments.

# Code Structure
# Main Script:

The main block uses an argument parser (argparse) to handle command-line inputs. It then calls a series of functions to process the files and generate outputs.

Placeholder Classes:
The script includes placeholder classes for BGC and BGCCollection. If you have the BGClib module, replace these placeholders with the actual implementations.

Core Functions:

get_files(): Reads GenBank files and aggregates BGC information.

assign_colors_to_domains(): Generates and assigns colors for gene domain annotations.

generate_combined_svg_and_html(): Creates combined SVG images and an interactive HTML visualization.

generate_combined_product_bar_charts(): Produces bar charts for product distributions.

save_fasta(): Outputs protein sequences in FASTA format.

write_metadata(): Writes out summarized and detailed metadata files.

# Example Outputs
HTML Visualization:
An interactive web page (BGCs_interactive.html) will be generated in the output folder, showing each BGC with interactive annotations and gene track details.

SVG Files:
Combined SVG files (e.g., combined_BGCs_1.svg, combined_BGCs_2.svg, etc.) are created to visually represent clusters.

Bar Charts:
A bar chart HTML file (product_distribution_combined.html) summarizes the distribution of candidate products by taxonomic levels.

FASTA & Metadata:
If enabled, the script outputs a FASTA file (proteins_with_CAL_domain.fasta) and metadata files with detailed information regarding the BGCs.

# License

MIT License

# Contact

For questions or feedback, please contact mattoslmp@gmail.com

