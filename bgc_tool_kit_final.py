import os
import argparse
from pathlib import Path
from shutil import copyfile
from Bio import SeqIO
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from lxml import etree
import numpy as np
from Bio import Entrez
import colorsys
import pandas as pd
import time

# === Placeholder Classes for BGClib ===
# If you have the BGClib module, replace these classes with the actual implementations.
class BGC:
    def __init__(self):
        self.protein_list = []
        self.loci = []
        self.identifier = ""

class BGCCollection:
    def __init__(self):
        self.bgcs = {}
# === End of Placeholder Classes ===

# Specific color for AMP-binding Proteins
AMP_BINDING_COLOR = "rgb(255,0,0)"  # Bright red

# Color generator to ensure unique colors
def generate_unique_colors():
    """
    Generates an infinite sequence of unique colors using the HSL color space,
    excluding the AMP_BINDING_COLOR.

    Yields:
        str: Color in 'rgb(r,g,b)' format.
    """
    num_colors = 1000  # Set a large number to avoid exhausting colors
    for i in range(num_colors):
        hue = i / num_colors
        lightness = 0.5
        saturation = 0.7
        r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
        rgb = f'rgb({int(r*255)},{int(g*255)},{int(b*255)})'
        if rgb != AMP_BINDING_COLOR:
            yield rgb

# Configure your email to use Entrez
Entrez.email = "your_email@example.com"  # Replace with your real email

def CMD_parser():
    parser = argparse.ArgumentParser(
        description="BGC toolkit with interactive HTML visualization. Tools to facilitate biosynthetic gene cluster handling and visualization.")
    parser.add_argument("-i", "--inputfolders", nargs='+', type=Path, required=True,
                        help="Folder(s) to search for .gbk files.")
    parser.add_argument("-o", "--outputfolder", default="./output", type=Path,
                        help="Base folder where results will be put.")
    parser.add_argument("--svg", action="store_true",
                        help="Enable SVG output for each BGC.")
    parser.add_argument("--svgcfg", default="SVG_arrow_options.cfg", type=Path,
                        help="Configuration file with SVG style.")
    parser.add_argument("--genbank", action="store_true",
                        help="Save GenBank files.")
    parser.add_argument("--metadata", type=str,
                        help="Write metadata files.")
    parser.add_argument("--fasta", action="store_true",
                        help="Save protein sequences in FASTA format.")
    parser.add_argument("--html", action="store_true",
                        help="Generate interactive HTML visualization.")
    parser.add_argument("--svgall", action="store_true",
                        help="Generate combined SVGs with up to 10 BGCs each.")
    return parser

def contains_cal_domain(feature) -> bool:
    """
    Checks if a feature contains a CAL_domain.

    Args:
        feature (SeqFeature): A feature from the GenBank record.

    Returns:
        bool: True if CAL_domain is present, False otherwise.
    """
    cal_keywords = ['CAL_domain']
    for key, values in feature.qualifiers.items():
        for value in values:
            if any(kw in value for kw in cal_keywords):
                return True
    return False

class BGCLocus:
    def from_record(self, record):
        self.id = record.id
        self.description = record.description
        self.features = record.features

class Protein:
    def from_feature(self, feature, record):
        self.id = feature.qualifiers.get('protein_id', [''])[0]
        self.identifier = feature.qualifiers.get('locus_tag', [''])[0]
        self.nucleotide_sequence = str(feature.extract(record.seq))
        self.amino_acid_sequence = str(feature.extract(record.seq).translate())
        self.role = "biosynthetic" if contains_cal_domain(feature) else "other"

def extract_orig_start_end_from_comment(record):
    orig_start = None
    orig_end = None

    structured_comment = record.annotations.get("structured_comment", {})
    for key, value in structured_comment.items():
        if "antiSMASH-Data" in key:
            orig_start = value.get("Orig. start", None)
            orig_end = value.get("Orig. end", None)
            break

    return orig_start, orig_end

def extract_species_name_from_definition(record):
    """
    Extracts the genome ID and species name from the DEFINITION line in the GenBank file.

    Args:
        record (SeqRecord): Biopython SeqRecord object.

    Returns:
        tuple: Genome ID and species name.
    """
    definition = record.description
    parts = definition.split()

    # Assuming the genome ID is the first part and the species name is the rest
    genome_id = parts[0]
    species_name = " ".join(parts[1:])  # The rest of the words should be the species name

    return genome_id, species_name

def extract_pfam_annotations(record):
    """
    Extracts PFAM annotations from a GenBank record.

    Args:
        record (SeqRecord): Biopython SeqRecord object.

    Returns:
        dict: Maps locus_tag to a tuple (aSDomain, description).
    """
    pfam_annotations = {}
    for feature in record.features:
        if feature.type == "PFAM_domain":
            locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
            asdomain = feature.qualifiers.get("aSDomain", [""])[0]
            description = feature.qualifiers.get("description", [""])[0]
            pfam_annotations[locus_tag] = (asdomain, description)
    return pfam_annotations

def extract_candidate_cluster_product(record):
    """
    Extracts the candidate cluster product from the features in the record.

    Args:
        record (SeqRecord): Biopython SeqRecord object.

    Returns:
        str: Product name or "Unknown" if not found.
    """
    for feature in record.features:
        if feature.type == "region" and "candidate_cluster_numbers" in feature.qualifiers:
            return feature.qualifiers.get("product", ["Unknown"])[0]
    return "Unknown"

def get_taxonomy(genus, taxonomic_rank):
    """
    Retrieves the taxonomic rank (e.g., Phylum, Order) for a given genus using NCBI Entrez.

    Args:
        genus (str): Genus name.
        taxonomic_rank (str): Desired taxonomic rank (e.g., 'phylum', 'order').

    Returns:
        str: Taxonomic rank name or "Unknown" if not found.
    """
    try:
        handle = Entrez.esearch(db="taxonomy", term=genus + "[genus]")
        record = Entrez.read(handle)
        handle.close()
        if record["Count"] == "0":
            return "Unknown"
        tax_id = record["IdList"][0]
        handle = Entrez.efetch(db="taxonomy", id=tax_id)
        records = Entrez.read(handle)
        handle.close()
        for lineage in records[0]["LineageEx"]:
            if lineage["Rank"] == taxonomic_rank:
                return lineage["ScientificName"]
        return "Unknown"
    except Exception as e:
        print(f"Error fetching taxonomy for genus {genus}: {e}")
        return "Unknown"

def get_files(inputfolders):
    """
    Processes GenBank files in the input folders and collects information on BGCs with CAL_domain.

    Args:
        inputfolders (list of Path): List of folders to search for .gbk files.

    Returns:
        tuple: (BGCCollection, dict, dict, dict, dict, dict)
            - BGCCollection: Collection of found BGCs.
            - dict: Maps bgc_id to GenBank file path.
            - dict: Maps bgc_id to species name.
            - dict: Genus and product counts for bar chart.
            - dict: Phylum and product counts.
            - dict: Order and product counts.
    """
    bgc_col = BGCCollection()
    gbk_files = {}
    species_map = {}  # Maps bgc_id to species name
    unique_bgc_entries = set()
    genus_product_counts = {}    # Product counts per genus
    phylum_product_counts = {}   # Product counts per phylum
    order_product_counts = {}    # Product counts per order

    genus_to_phylum = {}  # Cache to store genus to phylum mapping
    genus_to_order = {}   # Cache to store genus to order mapping

    for folder in inputfolders:
        for gbk_file in Path(folder).rglob("*.gbk"):
            print(f"Processing file: {gbk_file}")
            try:
                records = list(SeqIO.parse(gbk_file, "genbank"))
                for record in records:
                    genome_id, species_name = extract_species_name_from_definition(record)  # Extract genome ID and species name

                    # Correctly extract genus name
                    genus = "Unknown"  # Default genus
                    if species_name:
                        if 'MAG:' in species_name:
                            # Extract genus after 'MAG:'
                            mag_split = species_name.split('MAG:')
                            if len(mag_split) > 1:
                                genus_candidate = mag_split[1].strip().split()
                                if len(genus_candidate) > 0:
                                    genus = genus_candidate[0]
                        else:
                            genus_candidate = species_name.strip().split()
                            if len(genus_candidate) > 0:
                                genus = genus_candidate[0]
                    else:
                        genus = "Unknown"

                    # Get Phylum and Order for the genus
                    if genus not in genus_to_phylum or genus not in genus_to_order:
                        phylum = get_taxonomy(genus, 'phylum')
                        order = get_taxonomy(genus, 'order')
                        genus_to_phylum[genus] = phylum
                        genus_to_order[genus] = order
                        # Be polite to NCBI servers
                        time.sleep(0.5)
                    else:
                        phylum = genus_to_phylum[genus]
                        order = genus_to_order[genus]

                    candidate_product = extract_candidate_cluster_product(record)
                    for feature in record.features:
                        if feature.type == "CDS" and contains_cal_domain(feature):
                            bgc_entry = (record.id, feature.location.start, feature.location.end)
                            if bgc_entry not in unique_bgc_entries:
                                unique_bgc_entries.add(bgc_entry)
                                print(f"Found CAL_domain in {record.id} from {gbk_file}")
                                bgc = BGC()
                                bgc.protein_list = []
                                locus = BGCLocus()
                                locus.from_record(record)
                                bgc.loci = [locus]
                                bgc.identifier = f"{record.id}_{feature.location.start}_{feature.location.end}"
                                protein = Protein()
                                protein.from_feature(feature, record)
                                bgc.protein_list.append(protein)
                                bgc_col.bgcs[bgc.identifier] = bgc
                                gbk_files[bgc.identifier] = gbk_file
                                species_map[bgc.identifier] = (genome_id, species_name)  # Map bgc_id to species and genome ID

                                # Update genus and product counts
                                if genus not in genus_product_counts:
                                    genus_product_counts[genus] = {}
                                if candidate_product not in genus_product_counts[genus]:
                                    genus_product_counts[genus][candidate_product] = 0
                                genus_product_counts[genus][candidate_product] += 1

                                # Update phylum and product counts
                                if phylum not in phylum_product_counts:
                                    phylum_product_counts[phylum] = {}
                                if candidate_product not in phylum_product_counts[phylum]:
                                    phylum_product_counts[phylum][candidate_product] = 0
                                phylum_product_counts[phylum][candidate_product] += 1

                                # Update order and product counts
                                if order not in order_product_counts:
                                    order_product_counts[order] = {}
                                if candidate_product not in order_product_counts[order]:
                                    order_product_counts[order][candidate_product] = 0
                                order_product_counts[order][candidate_product] += 1
            except Exception as e:
                print(f"Error processing file {gbk_file}: {e}")
    return bgc_col, gbk_files, species_map, genus_product_counts, phylum_product_counts, order_product_counts

def assign_colors_to_domains(bgc_col, gbk_files):
    """
    Assigns unique colors to each asDomain description found in the BGCs.
    Ensures that each unique asDomain description has a unique color, and identical descriptions share the same color.

    Args:
        bgc_col (BGCCollection): Collection of BGCs.
        gbk_files (dict): Maps bgc_id to GenBank file path.
    """
    global domain_color_mapping
    domain_color_mapping = {}  # Reset mapping

    color_generator = generate_unique_colors()
    used_colors = set()

    # Set to store unique descriptions
    unique_descriptions = set()

    # First, collect all unique descriptions, including AMP-binding
    for bgc_id, bgc in bgc_col.bgcs.items():
        gbk_file = gbk_files[bgc_id]
        for record in SeqIO.parse(gbk_file, "genbank"):
            pfam_annotations = extract_pfam_annotations(record)
            for locus_tag, (asdomain, description) in pfam_annotations.items():
                if description:
                    key = description.strip().lower()
                    unique_descriptions.add(key)

    # Assign unique colors to each unique description
    for description in unique_descriptions:
        if 'amp-binding' in description.lower():
            domain_color_mapping[description] = AMP_BINDING_COLOR
        else:
            try:
                color = next(color_generator)
                while color in used_colors or color == AMP_BINDING_COLOR:
                    color = next(color_generator)
            except StopIteration:
                # Generate a random color if the generator is exhausted
                color = f'rgb({np.random.randint(0,256)},{np.random.randint(0,256)},{np.random.randint(0,256)})'
            domain_color_mapping[description] = color
            used_colors.add(color)

def write_metadata(o, metadata_base, bgc_col):
    """
    Writes summarized and detailed metadata files.

    Args:
        o (Path): Output folder.
        metadata_base (str): Base name for metadata files.
        bgc_col (BGCCollection): Collection of BGCs.
    """
    with open(o / f"{metadata_base}.metadata.summary.txt", "w") as s:
        s.write(f"{metadata_base} summary file\n\n")
        s.write("This collection contains\n")
        s.write(f"* {len(bgc_col.bgcs)} BGCs\n")

    with open(o / f"{metadata_base}.metadata.BGCs.tsv", "w") as m:
        header = "BGC\tCore Biosynthetic Protein type\tCore Biosynthetic Protein IDs\tCore Biosynthetic Protein Identifiers\n"
        m.write(header)
        for bgc_id in sorted(bgc_col.bgcs):
            bgc = bgc_col.bgcs[bgc_id]
            list_core_types = []
            list_protein_ids = []
            list_protein_identifiers = []
            for protein in bgc.protein_list:
                if protein.role == "biosynthetic":
                    list_core_types.append("CAL_domain")
                    if protein.id:
                        list_protein_ids.append(protein.id)
                    list_protein_identifiers.append(protein.identifier)
            m.write("{}\t{}\t{}\t{}\n".format(bgc.identifier,
                ", ".join(list_core_types),
                ", ".join(list_protein_ids),
                ", ".join(list_protein_identifiers)))

    with open(o / f"{metadata_base}.metadata.CBPs.tsv", "w") as c:
        header = "BGC\tCore Biosynthetic Protein type\tProtein identifier\tProtein Id\tNucleotide Sequence\tAmino Acid Sequence\n"
        c.write(header)
        for bgc_id in sorted(bgc_col.bgcs):
            bgc = bgc_col.bgcs[bgc_id]
            for protein in bgc.protein_list:
                if protein.role == "biosynthetic":
                    c.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(bgc_id,
                        "CAL_domain", protein.identifier, protein.id, protein.nucleotide_sequence, protein.amino_acid_sequence))

def add_start_end_to_html(fig, row, start_value, end_value, start_pos, end_pos, y_position):
    """
    Adds start and end values to the HTML visualization below the genes.

    Args:
        fig (go.Figure): The Plotly figure where annotations will be added.
        row (int): The row/subplot number where the BGC is plotted.
        start_value (str): Start value.
        end_value (str): End value.
        start_pos (int): x-position of the start.
        end_pos (int): x-position of the end.
        y_position (float): Y-position for the annotations.
    """
    xref = 'x' if row == 1 else f'x{row}'
    yref = 'y' if row == 1 else f'y{row}'

    # Add Start value annotation below the genes
    if start_value:
        fig.add_annotation(
            dict(
                text=f"Start: {start_value}",
                x=start_pos,
                y=y_position,
                showarrow=False,
                font=dict(size=12, color="black"),
                align="center",
                xanchor='center',
                yanchor='top',
                xref=xref,
                yref=yref
            )
        )

    # Add End value annotation below the genes
    if end_value:
        fig.add_annotation(
            dict(
                text=f"End: {end_value}",
                x=end_pos,
                y=y_position,
                showarrow=False,
                font=dict(size=12, color="black"),
                align="center",
                xanchor='center',
                yanchor='top',
                xref=xref,
                yref=yref
            )
        )

def add_bgc_to_html(fig, row, bgc_id, gbk_file, genome_id, species_name, candidate_product, candidate_region, pfam_annotations, orig_start, orig_end, total_included_bgcs, gene_scaling_factor):
    """
    Adds a specific BGC to the HTML visualization using Plotly.

    Args:
        fig (go.Figure): The Plotly figure where the BGC will be added.
        row (int): The row/subplot number where the BGC is plotted.
        bgc_id (str): Unique identifier of the BGC.
        gbk_file (Path): Path to the GenBank file.
        genome_id (str): Genome ID.
        species_name (str): Species name.
        candidate_product (str): Candidate product.
        candidate_region (SeqFeature): Candidate region.
        pfam_annotations (dict): PFAM annotations.
        orig_start (str): Original start value.
        orig_end (str): Original end value.
        total_included_bgcs (int): Total number of included BGCs.
        gene_scaling_factor (float): Gene scaling factor.
    """
    gene_positions = []
    asdomain_descriptions = []
    seen_descriptions = set()

    # Adjust Y positions for annotations and genes
    # Defining a Y range from 0 to 1.2 for each subplot
    genome_species_y = 1.0      # Y position for the title (Genome and Species)
    product_name_y = 0.9        # Y position for product name
    gene_y_position = 0.5       # Y position for genes
    start_end_y = 0.2           # Y position for Start and End values

    amp_binding_present = False

    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and feature.location.start >= candidate_region.start and feature.location.end <= candidate_region.end:
                start = int((feature.location.start - candidate_region.start) * gene_scaling_factor)
                end = int((feature.location.end - candidate_region.start) * gene_scaling_factor)
                strand = feature.location.strand
                locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]

                asdomain, description = pfam_annotations.get(locus_tag, ("", ""))
                is_amp_binding = 'amp-binding' in description.lower()

                if is_amp_binding:
                    amp_binding_present = True

                if description:
                    asdomain_descriptions.append(description)
                else:
                    asdomain_descriptions.append("Unknown")

                # Define color
                if is_amp_binding:
                    color = AMP_BINDING_COLOR
                else:
                    key = description.strip().lower()
                    color = domain_color_mapping.get(key, "rgb(169,169,169)")  # Default gray

                hover_text = f"Contig: {locus_tag}<br>asDomain: {asdomain}<br>Description: {description if description else 'Unknown'}<br>Strand: {'+' if strand == 1 else '-'}"

                # Determine gene direction based on strand
                arrow_symbol = 'triangle-up' if strand == 1 else 'triangle-down'

                fig.add_trace(
                    go.Scatter(
                        x=[start, end],
                        y=[gene_y_position, gene_y_position],
                        mode='lines+markers',
                        line=dict(color=color, width=12),
                        marker=dict(symbol=arrow_symbol, size=10, color=color),
                        name=locus_tag,
                        text=hover_text,
                        hoverinfo="text"
                    ),
                    row=row, col=1
                )

                gene_positions.append((start, end))

    # If AMP-binding domain is not present, return False
    if not amp_binding_present:
        return False  # Do not include this BGC

    xref = 'x' if row == 1 else f'x{row}'
    yref = 'y' if row == 1 else f'y{row}'

    # Set Y-axis range to accommodate all elements
    y_range = [0, 1.2]  # Increase Y range to accommodate the title above
    fig.update_yaxes(range=y_range, row=row, col=1)

    # Add start and end values below the genes
    if gene_positions:
        # Position of the first gene
        start_pos = gene_positions[0][0]
        # Position of the last gene
        end_pos = gene_positions[-1][1]

        # Add start and end values below the genes in HTML
        add_start_end_to_html(
            fig,
            row=row,
            start_value=orig_start,
            end_value=orig_end,
            start_pos=start_pos,
            end_pos=end_pos,
            y_position=start_end_y
        )

    # Add title with genome ID and species
    fig.add_annotation(
        dict(
            text=f"Genome: {genome_id} | Species: {species_name}",
            x=0.5,
            y=genome_species_y,
            xref=f'{xref} domain',
            yref=f'{yref} domain',
            showarrow=False,
            font=dict(size=14, color="black"),
            xanchor='center',
            yanchor='top'
        ),
        row=row, col=1
    )

    # Add product text below the title
    fig.add_annotation(
        dict(
            text=f"Product: {candidate_product}",
            x=0.5,
            y=product_name_y,
            xref=f'{xref} domain',
            yref=f'{yref} domain',
            showarrow=False,
            font=dict(size=12, color="black"),
            xanchor='center',
            yanchor='top'
        ),
        row=row, col=1
    )

    # Update subplot layout to remove titles and axes
    fig.update_yaxes(title_text="", row=row, col=1, showgrid=False, zeroline=False, visible=False)
    fig.update_xaxes(title_text="", row=row, col=1, showgrid=False, zeroline=False, visible=False)
    fig.update_layout(showlegend=False)

    return asdomain_descriptions

def generate_combined_svg_and_html(o, bgc_col, gbk_files, species_map, bgcs_per_svg=10):
    """
    Generates combined SVGs and an interactive HTML visualization containing all BGCs.

    Args:
        o (Path): Output folder.
        bgc_col (BGCCollection): Collection of BGCs.
        gbk_files (dict): Maps bgc_id to GenBank file path.
        species_map (dict): Maps bgc_id to species name.
        bgcs_per_svg (int): Maximum number of BGCs per SVG.
    """
    assign_colors_to_domains(bgc_col, gbk_files)  # Assign colors based on asDomain descriptions
    margin = 50
    right_margin = 100
    svg_width = 1400

    # Collect all BGCs to include
    bgcs_to_include = []
    for bgc_id in bgc_col.bgcs.keys():
        gbk_file = gbk_files[bgc_id]
        candidate_region = None
        candidate_product = None
        amp_binding_present_in_bgc = False

        for record in SeqIO.parse(gbk_file, "genbank"):
            candidate_product = extract_candidate_cluster_product(record)

            # Skip BGC if the product is "Unknown"
            if candidate_product == "Unknown":
                print(f"Skipping BGC {bgc_id} with unknown product.")
                break  # Skip to the next record

            for feature in record.features:
                if feature.type == "region" and "candidate_cluster_numbers" in feature.qualifiers:
                    candidate_region = feature.location
                    break

            # Check for AMP-binding domain in the BGC
            pfam_annotations = extract_pfam_annotations(record)
            for feature in record.features:
                if feature.type == "CDS" and feature.location.start >= candidate_region.start and feature.location.end <= candidate_region.end:
                    locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    _, description = pfam_annotations.get(locus_tag, ("", ""))
                    if 'amp-binding' in description.lower():
                        amp_binding_present_in_bgc = True
                        break  # No need to check further

            if amp_binding_present_in_bgc:
                bgcs_to_include.append(bgc_id)
            else:
                print(f"Skipping BGC {bgc_id} without AMP-binding domain.")

    total_included_bgcs = len(bgcs_to_include)
    print(f"Total BGCs to include: {total_included_bgcs}")

    if total_included_bgcs == 0:
        print("No BGCs to include. Exiting.")
        return

    # Decide how many SVGs to generate
    if total_included_bgcs <= bgcs_per_svg:
        num_svgs = 1
        bgc_groups = [bgcs_to_include]
    else:
        num_svgs = (total_included_bgcs + bgcs_per_svg - 1) // bgcs_per_svg  # Ceiling division
        bgc_groups = [bgcs_to_include[i:i + bgcs_per_svg] for i in range(0, total_included_bgcs, bgcs_per_svg)]

    # Initialize the Plotly figure for all BGCs
    max_vertical_spacing = 1 / (total_included_bgcs - 1) if total_included_bgcs > 1 else 1
    desired_vertical_spacing = 0.005  # Reduced value to minimize overlap
    vertical_spacing = min(desired_vertical_spacing, max_vertical_spacing - 0.001)  # Subtract a small value to avoid exact limits

    fig = make_subplots(rows=total_included_bgcs, cols=1, subplot_titles=None,
                        vertical_spacing=vertical_spacing)  # Adjust to control space between subplots

    row_number = 1  # Start row number for the HTML
    svg_index = 1

    for group in bgc_groups:
        y_offset = margin  # Initial vertical position
        space_between_bgcs = 50  # Space between BGCs

        # First, calculate the total height needed for the SVG
        combined_svg_height = 0
        for bgc_id in group:
            # Estimate the height per BGC (adjust as needed)
            combined_svg_height += 400  # Estimated height per BGC

        combined_svg = etree.Element("svg", attrib={
            "version": "1.1",
            "baseProfile": "full",
            "width": str(svg_width + right_margin),
            "height": str(combined_svg_height + margin + 400),  # Increase height for extra spacing
            "xmlns": "http://www.w3.org/2000/svg",
            "style": "background-color:white;"  # White background
        })

        # Draw each BGC in the SVG and add to HTML
        for bgc_id in group:
            gbk_file = gbk_files[bgc_id]
            candidate_region = None
            orig_start = None
            orig_end = None
            pfam_annotations = {}
            candidate_product = None

            genome_id, species_name = species_map[bgc_id]
            genome_name = Path(gbk_file).stem  # Assuming genome name is the filename without extension

            for record in SeqIO.parse(gbk_file, "genbank"):
                if not orig_start or not orig_end:
                    orig_start, orig_end = extract_orig_start_end_from_comment(record)

                pfam_annotations = extract_pfam_annotations(record)

                if not candidate_product:
                    candidate_product = extract_candidate_cluster_product(record)

                for feature in record.features:
                    if feature.type == "region" and "candidate_cluster_numbers" in feature.qualifiers:
                        candidate_region = feature.location
                        break

            if candidate_region is None:
                continue

            total_gene_length = 0
            num_genes = 0
            for record in SeqIO.parse(gbk_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS" and feature.location.start >= candidate_region.start and feature.location.end <= candidate_region.end:
                        gene_length = len(feature.location)
                        total_gene_length += gene_length
                        num_genes += 1

            gene_scaling_factor = (svg_width - margin - right_margin) / total_gene_length if total_gene_length > 0 else 1
            arrow_height = 12
            arrow_spacing = 35
            x_offset = margin

            first_feature_start = None
            last_feature_end = None

            # Define y_gene_start BEFORE processing genes
            y_gene_start = y_offset + 116  # Start of genes (just below the product)

            # Add BGC to the HTML plot
            asdomain_descriptions = add_bgc_to_html(
                fig,
                row=row_number,
                bgc_id=bgc_id,
                gbk_file=gbk_file,
                genome_id=genome_id,
                species_name=species_name,
                candidate_product=candidate_product,
                candidate_region=candidate_region,
                pfam_annotations=pfam_annotations,
                orig_start=orig_start,
                orig_end=orig_end,
                total_included_bgcs=total_included_bgcs,
                gene_scaling_factor=gene_scaling_factor
            )

            # If the BGC does not have AMP-binding domain, skip it
            if not asdomain_descriptions:
                continue

            # Process again to add to the SVG
            descriptions_in_order = []  # To store descriptions in order
            gene_positions = []

            # Add title for each BGC in the SVG (2 cm above the BGC)
            title_text = etree.SubElement(combined_svg, "text", attrib={
                "x": str(svg_width // 2),
                "y": str(y_offset + 76),  # 2 cm in pixels (76 pixels)
                "text-anchor": "middle",
                "font-size": "16",
                "fill": "black",
                "font-weight": "bold"
            })
            title_text.text = f"Genome: {genome_id} | Species: {species_name}"

            # Add product text below the title
            product_text = etree.SubElement(combined_svg, "text", attrib={
                "x": str(svg_width // 2),
                "y": str(y_offset + 96),  # Slightly below the title
                "text-anchor": "middle",
                "font-size": "14",
                "fill": "black",
                "font-weight": "bold"
            })
            product_text.text = f"Product: {candidate_product}"

            for record in SeqIO.parse(gbk_file, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS" and feature.location.start >= candidate_region.start and feature.location.end <= candidate_region.end:
                        if first_feature_start is None:
                            first_feature_start = int(feature.location.start)
                        last_feature_end = int(feature.location.end)

                        start = int((feature.location.start - candidate_region.start) * gene_scaling_factor)
                        end = int((feature.location.end - candidate_region.start) * gene_scaling_factor)
                        strand = feature.location.strand
                        locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]

                        asdomain, description = pfam_annotations.get(locus_tag, ("", ""))
                        is_amp_binding = 'amp-binding' in description.lower()
                        if description:
                            descriptions_in_order.append(description)
                        else:
                            descriptions_in_order.append("Unknown")
                        description = description if description else "Unknown"

                        # Define color
                        if is_amp_binding:
                            color = AMP_BINDING_COLOR
                        else:
                            key = description.strip().lower()
                            color = domain_color_mapping.get(key, "rgb(169,169,169)")

                        group_elem = etree.SubElement(combined_svg, "g")
                        points = ""
                        if strand == 1:
                            points = f"{x_offset + start},{y_gene_start} {x_offset + end - 12},{y_gene_start} {x_offset + end},{y_gene_start + arrow_height / 2} {x_offset + end - 12},{y_gene_start + arrow_height} {x_offset + start},{y_gene_start + arrow_height}"
                        else:
                            points = f"{x_offset + end},{y_gene_start} {x_offset + start + 12},{y_gene_start} {x_offset + start},{y_gene_start + arrow_height / 2} {x_offset + start + 12},{y_gene_start + arrow_height} {x_offset + end},{y_gene_start + arrow_height}"

                        arrow = etree.SubElement(group_elem, "polygon", attrib={
                            "points": points,
                            "fill": color,
                            "stroke": "black",
                            "stroke-width": "1"
                        })

                        gene_positions.append((start, end))

            # Add start and end values below the genes in the SVG
            if gene_positions:
                # Position of the first gene
                start_pos = gene_positions[0][0]
                # Position of the last gene
                end_pos = gene_positions[-1][1]

                if orig_start:
                    # Add Start text in the SVG
                    start_text = etree.SubElement(combined_svg, "text", attrib={
                        "x": str(x_offset + start_pos),
                        "y": str(y_gene_start + arrow_height + 20),
                        "text-anchor": "middle",
                        "font-size": "12",
                        "fill": "black"
                    })
                    start_text.text = f"Start: {orig_start}"

                if orig_end:
                    # Add End text in the SVG
                    end_text = etree.SubElement(combined_svg, "text", attrib={
                        "x": str(x_offset + end_pos),
                        "y": str(y_gene_start + arrow_height + 20),
                        "text-anchor": "middle",
                        "font-size": "12",
                        "fill": "black"
                    })
                    end_text.text = f"End: {orig_end}"

            # Add asDomain legend below the BGC
            legend_x = margin
            legend_y = y_gene_start + arrow_height + 40  # Adjust vertical position for the legend

            num_legend_items = len(descriptions_in_order)
            # Adjust the legend to display all items in order
            for idx, description in enumerate(descriptions_in_order):
                key = description.strip().lower()
                if 'amp-binding' in key:
                    color = AMP_BINDING_COLOR
                else:
                    color = domain_color_mapping.get(key, "rgb(169,169,169)")

                # Positioning for each legend item
                legend_item_x = legend_x + (idx % 5) * 250  # Adjust horizontal spacing
                legend_item_y = legend_y + (idx // 5) * 50  # Adjust vertical spacing (increased to 50)

                # Add colored square
                rect = etree.SubElement(combined_svg, "rect", attrib={
                    "x": str(legend_item_x),
                    "y": str(legend_item_y),
                    "width": "15",
                    "height": "15",
                    "fill": color,
                    "stroke": "black",
                    "stroke-width": "1"
                })

                # Adjust description text if it has more than 3 words
                words = description.split()
                if len(words) > 3:
                    # Insert newline after the third word
                    first_line = ' '.join(words[:3])
                    second_line = ' '.join(words[3:])
                    description_lines = [first_line, second_line]
                else:
                    description_lines = [description]

                # Add description text with possible line breaks
                text_elem = etree.SubElement(combined_svg, "text", attrib={
                    "x": str(legend_item_x + 25),
                    "y": str(legend_item_y + 12),
                    "font-size": "12",
                    "fill": "black"
                })
                for i, line in enumerate(description_lines):
                    tspan = etree.SubElement(text_elem, "tspan", attrib={
                        "x": str(legend_item_x + 25),
                        "dy": str(12 if i > 0 else 0)  # Offset for subsequent lines
                    })
                    tspan.text = line

            # Update y_offset for the next figure
            # Adjust legend height based on the maximum number of lines in descriptions
            max_lines_in_legend = max(len(desc.split()) // 3 + 1 for desc in descriptions_in_order)
            legend_height = ((len(descriptions_in_order) + 4) // 5) * 50 + (max_lines_in_legend - 1) * 12  # Height of the legend area
            y_offset = legend_y + legend_height + space_between_bgcs  # Space before the next BGC

            # Update row number for the next BGC
            row_number += 1

        # Increment svg_index only if an SVG was actually generated
        if group:
            combined_svg_path = o / f"combined_BGCs_{svg_index}.svg"
            with open(combined_svg_path, "wb") as svg_file:
                svg_file.write(etree.tostring(combined_svg, pretty_print=True))
            print(f"Generated SVG file: {combined_svg_path}")
            svg_index += 1

    # Adjust overall layout of the Plotly figure
    fig.update_layout(
        height=300 * total_included_bgcs,  # Increase height for more vertical space
        width=1500,                        # Increase width
        title_text="Biosynthetic Gene Clusters (BGCs) Interactive Visualization",
        showlegend=False,                  # Remove legends from HTML
        hovermode='closest',
        paper_bgcolor="white",             # White background
        plot_bgcolor="white",              # White background
        margin=dict(l=50, r=50, t=50, b=50)  # Adjust margins
    )

    # Save the Plotly figure as HTML if there are included BGCs
    if total_included_bgcs > 0:
        html_output_path = o / "BGCs_interactive.html"
        fig.write_html(str(html_output_path))
        print(f"Generated interactive HTML file: {html_output_path}")

def save_fasta(o, bgc_col):
    """
    Saves protein sequences with CAL_domain in FASTA format.

    Args:
        o (Path): Output folder.
        bgc_col (BGCCollection): Collection of BGCs.
    """
    fasta_name = o / "proteins_with_CAL_domain.fasta"
    with open(fasta_name, "w") as fasta_file:
        for bgc_id, bgc in bgc_col.bgcs.items():
            for protein in bgc.protein_list:
                if protein.role == "biosynthetic":
                    fasta_file.write(f">{protein.identifier} [BGC={bgc_id}]\n{protein.amino_acid_sequence}\n")

def generate_combined_product_bar_charts(o, genus_counts, order_counts, phylum_counts):
    """
    Generates a single HTML file containing bar charts for Phylum, Order, and Genus.

    Args:
        o (Path): Output folder.
        genus_counts (dict): Product counts per genus.
        order_counts (dict): Product counts per order.
        phylum_counts (dict): Product counts per phylum.
    """
    # Create a figure with three subplots
    fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=("Phylum", "Order", "Genus"),
        vertical_spacing=0.15
    )

    # List of counts and taxonomic levels
    counts_list = [(phylum_counts, 'Phylum', 1), (order_counts, 'Order', 2), (genus_counts, 'Genus', 3)]

    for counts_dict, taxonomic_level, row in counts_list:
        # Collect all products
        all_products = set()
        for product_counts in counts_dict.values():
            all_products.update(product_counts.keys())
        all_products = sorted(all_products)

        # Create DataFrame
        df = pd.DataFrame(columns=all_products)
        for taxon, product_counts in counts_dict.items():
            df.loc[taxon] = [product_counts.get(product, 0) for product in all_products]
        df = df.fillna(0)

        # Sort taxa by total counts
        df['Total'] = df.sum(axis=1)
        df = df.sort_values('Total', ascending=False)
        df = df.drop('Total', axis=1)

        # Create the stacked bar chart
        for product in all_products:
            fig.add_trace(go.Bar(
                y=df.index,
                x=df[product],
                name=product,
                orientation='h',
                text=[f"{(count / df.loc[taxon].sum()) * 100:.1f}%" if df.loc[taxon].sum() > 0 else "0%" for taxon, count in zip(df.index, df[product])],
                textposition='inside',
                hovertemplate=f'<b>{taxonomic_level}:</b> %{{y}}<br><b>Product:</b> {product}<br><b>Count:</b> %{{x}}<br><b>Percentage:</b> %{{text}}<extra></extra>'
            ), row=row, col=1)

        # Update layout for each subplot
        fig.update_yaxes(title_text=taxonomic_level, row=row, col=1)
        fig.update_xaxes(title_text="Count", row=row, col=1)

    # Update the layout
    fig.update_layout(
        height=1800,
        width=900,
        title_text="Product Distribution by Phylum, Order, and Genus",
        barmode='stack',
        legend_title="Products",
        template='plotly_white',
        showlegend=True,
        hovermode='closest'
    )

    output_path = o / "product_distribution_combined.html"
    fig.write_html(str(output_path))
    print(f"Generated combined bar chart HTML file: {output_path}")

if __name__ == "__main__":
    args = CMD_parser().parse_args()

    output_folder = Path(args.outputfolder)
    output_folder.mkdir(parents=True, exist_ok=True)

    print("Collecting data and filtering for CAL_domain...")
    bgc_collection, gbk_files, species_map, genus_product_counts, phylum_product_counts, order_product_counts = get_files(args.inputfolders)
    print("...done")

    if args.svgall or args.html:
        print("SVG and HTML: Generating combined figure of all BGCs")
        generate_combined_svg_and_html(output_folder, bgc_collection, gbk_files, species_map)
        print("...done")

    # Generate combined bar chart of Products by Phylum, Order, and Genus
    print("Generating combined product distribution bar charts...")
    generate_combined_product_bar_charts(output_folder, genus_product_counts, order_product_counts, phylum_product_counts)
    print("...done")

    if args.genbank:
        print("Saving GenBank files...")
        for bgc_id, gbk_file in gbk_files.items():
            copyfile(gbk_file, output_folder / f"{bgc_id}.gbk")

    if args.fasta:
        print("Saving FASTA files...")
        save_fasta(output_folder, bgc_collection)

    if args.metadata:
        print("Writing metadata...")
        write_metadata(output_folder, args.metadata, bgc_collection)

    print("Processing complete.")
