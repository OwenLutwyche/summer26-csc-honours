# AI Usage Documentation - Week 5

## AI Model Used
- Claude Sonnet 4.5 via GitHub Copilot Chat

## AI Tools Used
- **GitHub Copilot**: Code completion and inline suggestions
- **GitHub Copilot Chat (Claude Sonnet 4.5)**: Complex problem-solving, debugging, and architecture design

## Components Developed with AI Assistance

### 1. Notebook Structure and Organization
I asked AI to help me design the overall structure of the Jupyter notebook. AI suggested:
- Logical section ordering (setup, download, alignment, variant calling, phasing, comparison, visualization, star-allele analysis)
- Markdown documentation cells to explain each step
- Clear section headers and subsections

### 2. Data Download and Processing Scripts
I needed bash scripts to download sequencing data and reference genomes. AI helped me create:
- wget commands for downloading Illumina and PacBio FASTQ files
- Reference genome download and decompression
- Conditional checks to avoid re-downloading existing files
- File indexing commands for minimap2 and samtools

### 3. Read Alignment with minimap2
I asked AI to help me write the alignment commands. AI provided:
- Appropriate minimap2 presets for each technology (sr for Illumina, map-hifi for PacBio)
- Piped commands to sort BAM files with samtools
- Indexing commands for the sorted BAM files
- Flagstat commands to extract alignment statistics

### 4. Variant Calling with bcftools
I needed to call variants in specific gene regions. AI helped with:
- bcftools mpileup command with region specification
- bcftools call parameters for variant calling
- VCF compression and indexing
- Commands to count variants

### 5. Variant Phasing with HapCUT2
This was more complex and required multiple iterations with AI. AI helped me:
- Write extractHAIRS commands to prepare fragment files
- Understand technology-specific parameters (PacBio requires --pacbio 1 and --ref)
- Execute HAPCUT2 for phasing
- Handle output file format issues

### 6. Variant Comparison Analysis
I asked AI to help me compare variants between the two technologies. AI wrote:
- Python code using pysam to parse VCF files
- Logic to categorize variants as shared, Illumina-only, or PacBio-only
- Dictionary structures to organize results by gene
- Formatted output showing comparison statistics

### 7. Automated IGV-like Visualization
This was the most complex component. I originally tried to use IGV batch scripts but the command wasn't available. I asked AI to create Python-based automated visualizations instead. AI developed:
- A matplotlib-based visualization function that mimics IGV's read view
- Color-coded bases with quality-based transparency
- Strand information (darker for reverse strand, lighter for forward)
- Variant position highlighting with a red dashed line
- Coverage depth plotting below the read alignment view
- Automatic calculation of variant support statistics
- PNG file saving and inline notebook display

### 8. Star-Allele Determination Framework
I asked AI to help create a framework for star-allele determination. AI provided:
- Code to extract phased haplotypes from VCF files
- Star-allele variant definitions based on PharmVar database
- Matching logic to identify diplotypes
- Basic star-allele assignment based on defining variants

### 9. Discussion Templates
I asked AI to create structured templates for the manual discussion sections required by the assignment. AI generated:
- Section 6.1: Variant comparison analysis template with guiding questions
- Section 7.3: IGV screenshot analysis template with structured evaluation criteria
- Section 8.2: Star-allele clinical interpretation template with per-gene frameworks

## My Contributions and Understanding

### What I Did Myself
- Identified gene coordinates from UCSC Genome Browser
- Understood the clinical significance of CYP genes in pharmacogenomics
- Made decisions about which variants to visualize
- Debugged issues by reading error messages and reporting them to AI
- Tested all code and verified outputs
- Validated that results make biological sense
- Will write all manual discussion sections analyzing the results

### What I Learned
- How to build a complete variant calling pipeline from raw reads to phased variants
- Differences between short-read and long-read sequencing technologies
- How variant phasing works and why it's important for pharmacogenomics
- Tool-specific requirements (e.g., HapCUT2 file format needs, minimap2 presets)
- How to use pysam for programmatic BAM and VCF file analysis
- How to create custom visualizations with matplotlib
- Importance of coverage in variant calling

## Problems Encountered and How I Solved Them

### Issue 1: extractHAIRS PacBio Parameter Requirements
**What happened**: When I ran the fragment extraction cell, it failed with an error message: "In order to realign variants (including --pacbio and --ont options), reference fasta file must be provided with --ref option"

**How I solved it**: 
- I reported the error to AI
- AI explained that PacBio long reads need local realignment for accurate fragment extraction
- AI helped me add the `--ref data/chr10.fa` parameter to the extractHAIRS command
- I re-ran the cell and it succeeded
- Illumina short reads don't need this parameter

### Issue 2: HAPCUT2 VCF Format Requirements
**What happened**: When I ran the phasing cell, HAPCUT2 failed with error: "VCF file has no header, ERROR reading VCF file: less than 10 columns"

**How I solved it**:
- I reported the error to AI
- AI explained that HAPCUT2 and extractHAIRS require uncompressed VCF files
- AI helped me add a decompression step: `gunzip -c results/illumina.vcf.gz > results/illumina.vcf`
- I updated both the Illumina and PacBio sections with this fix
- Re-ran the cells and they succeeded

### Issue 3: HAPCUT2 Output File Naming
**What happened**: My variant comparison cell failed with FileNotFoundError looking for `results/illumina.phased.vcf`, but the file existed as `results/illumina.phased.phased.VCF` (uppercase VCF)

**How I solved it**:
- I checked the results directory and found the oddly-named file
- I asked AI why this happened
- AI explained that HAPCUT2 always appends `.phased.VCF` to the output basename
- So `--output results/illumina.phased` creates `results/illumina.phased.phased.VCF`
- AI helped me add rename commands to convert to lowercase: `mv results/illumina.phased.phased.VCF results/illumina.phased.vcf`
- I also updated the file existence checks to look for the double-phased name first
- This fixed the issue and subsequent cells could find the files

### Issue 4: Understanding Empty Results for CYP2C9 and CYP2C8
**What happened**: My variant comparison showed 0 variants for CYP2C9 and CYP2C8. Only CYP2C19 had results (124-136 variants).

**How I investigated**:
- I was confused why two genes had no variants
- I asked AI to help me check the coverage
- AI showed me how to use `samtools view -c` to count reads in each gene region
- Results: CYP2C19 had thousands of reads, CYP2C9 had 1 read, CYP2C8 had 0 reads
- AI explained this is targeted sequencing data focused only on CYP2C19

**What I learned**:
- Not all sequencing experiments cover the entire genome
- Targeted sequencing enriches specific regions of interest
- You cannot call variants in regions without sequencing coverage
- This is normal and expected behavior

### Issue 5: Automated IGV Screenshot Generation
**What happened**: I tried to run IGV batch script with `igv -b results/igv_batch.txt` but got exit code 127 (command not found). IGV wasn't installed on my system.

**How I solved it**:
- I asked AI if there was a way to automate visualization without installing IGV
- AI suggested creating Python-based visualizations using matplotlib and pysam
- AI wrote a comprehensive visualization function that:
  - Reads BAM files and extracts reads around variant positions
  - Displays each read as a horizontal bar with color-coded bases
  - Uses color intensity to show base quality scores
  - Shows strand information (darker colors for reverse strand)
  - Adds a coverage depth plot below the read view
  - Calculates and displays variant support statistics
  - Saves PNG files and displays them inline in the notebook
- I tested the function and it worked perfectly
- This actually turned out better than IGV screenshots because the statistics are automatically calculated

### Issue 6: PIL/Pillow Import Conflict
**What happened**: When I tried to run the visualization cell, matplotlib failed to import with error: `ImportError: cannot import name '_imaging' from 'PIL'`

**How I debugged it**:
- I tried uninstalling and reinstalling Pillow with pip, but the error persisted
- I reported the error to AI with the full stack trace
- AI identified that my system has two versions of PIL/Pillow:
  - System PIL (10.2.0) at `/usr/lib/python3/dist-packages`
  - User-installed Pillow (12.0.0) at `~/.local/lib/python3.13/site-packages`
- Python was loading the older system version first due to sys.path ordering

**How I solved it**:
- AI suggested reordering sys.path to prioritize user-installed packages
- I added this code at the top of the visualization cell:
  ```python
  import sys
  sys.path = [p for p in sys.path if '/usr/lib/python3/dist-packages' not in p] + \
             [p for p in sys.path if '/usr/lib/python3/dist-packages' in p]
  ```
- AI also suggested setting matplotlib backend to 'Agg': `matplotlib.use('Agg')`
- I restarted the kernel and re-ran the cells
- This fixed the import issue

### Issue 7: pysam MD Tag Error
**What happened**: After fixing the PIL import, the visualization cell ran further but then failed with: `ValueError: MD tag not present`

**How I debugged it**:
- The error pointed to the line: `read.get_aligned_pairs(matches_only=False, with_seq=True)`
- I reported this new error to AI
- AI explained that the `with_seq=True` parameter requires an MD tag in the BAM file
- MD tags contain reference sequence information for each alignment
- minimap2 doesn't generate MD tags by default

**How I solved it**:
- AI suggested removing the `with_seq=True` parameter
- AI helped me refactor the code:
  - Changed from: `for query_pos, ref_pos, ref_base in aligned_pairs:`
  - To: `for pair in aligned_pairs:` then extract `query_pos, ref_pos = pair[0], pair[1]`
  - Read bases directly from `read.query_sequence[query_pos]` instead
- I updated the code and re-ran the cell
- The visualization now works perfectly without needing MD tags

### Issue 8: Cell Reorganization After Automation
**What happened**: After implementing the Python-based visualization, I realized the IGV batch script cell was no longer needed.

**How I cleaned it up**:
- I asked AI to help me remove the redundant IGV batch script cell (section 7.1)
- AI deleted both the batch script cell and its markdown header
- This kept the notebook clean and well-organized

### Issue 9: Creating Discussion Templates
**What I needed**: The assignment requires manual discussion and interpretation in three areas:
1. Analysis of variant concordance/discordance between technologies
2. Evaluation of whether discordant variants are true variants or artifacts
3. Star-allele determination with detailed clinical interpretation and rationale

**How I got help**:
- I asked AI to create structured templates to guide my manual analysis
- AI created three markdown cells with:
  - Section 6.1: Questions about variant comparison with prompts for my answers
  - Section 7.3: Structured evaluation framework for each visualized variant
  - Section 8.2: Templates for star-allele interpretation with clinical context
- The templates include:
  - Specific questions to answer
  - Checkboxes for my verdicts
  - Prompts for rationale and evidence
  - Space for clinical implications

### Issue 10: VCF Haplotype Extraction Error
**What happened**: When I ran the star-allele determination cell (section 8), it failed with: `ValueError: fetch requires an index`

**How I debugged it**:
- The error occurred in the `extract_phased_haplotypes()` function
- The code was using `vcf.fetch(coords['chr'], coords['start'], coords['end'])`
- I reported the error to AI
- AI explained that pysam's `fetch()` method requires an indexed VCF file (.tbi or .csi)
- The phased VCF files from HAPCUT2 were not indexed

**How I solved it**:
- AI suggested switching from indexed fetch to simple iteration
- AI helped me rewrite the extraction logic:
  - Old approach: `for record in vcf.fetch(coords['chr'], coords['start'], coords['end']):`
  - New approach: 
    ```python
    for record in vcf:
        if record.chrom != coords['chr']:
            continue
        if record.pos < coords['start'] or record.pos > coords['end']:
            continue
    ```
- This manually filters records by region without requiring an index
- I updated the cell and re-ran it
- The extraction now works successfully

**What I learned**:
- VCF fetch operations require indexed files
- Creating indexes is an extra step: `bcftools index file.vcf.gz` or `tabix`
- Iteration with manual filtering is simpler when indexes aren't available
- This approach works fine for small VCF files like our phased variants

## Summary of What AI Did vs What I Did

### AI Created/Automated
- Complete pipeline structure and code
- Data download and preprocessing scripts
- Read alignment commands with appropriate parameters
- Variant calling with region specification
- Variant phasing with HapCUT2
- Python code for variant comparison analysis
- Custom matplotlib visualization function (IGV-like)
- Coverage depth plotting
- Automatic variant support statistics
- Basic star-allele matching framework
- Discussion templates for manual analysis

### I Contributed
- Gene coordinates and clinical context
- Quality control and validation of outputs
- Debugging by reading error messages and testing
- Understanding of what each step does and why
- Decisions about which variants to visualize
- Verification that results make biological sense
- Will write all manual discussion sections:
  - Variant comparison interpretation
  - Assessment of whether variants are real or artifacts
  - Star-allele determinations with clinical rationale
  - Phenotype predictions and drug metabolism implications
  - Metabolizer phenotype assessment
  - Clinical implications and recommendations

## Notes
- All AI-generated code was reviewed and tested
- Bioinformatics tool parameters were validated against documentation
- Gene coordinates verified using UCSC Genome Browser
- Star-allele definitions cross-referenced with PharmVar database
- Pipeline iteratively debugged to handle tool-specific requirements
- File format and naming conventions discovered through testing
- Discussion templates designed to encourage critical thinking while providing structure
- Visualizations successfully replicate IGV functionality without external dependencies
