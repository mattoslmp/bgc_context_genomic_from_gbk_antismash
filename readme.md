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
