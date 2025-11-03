# AI Usage Documentation - Week 5

## Model
- Claude Sonnet 4.5

## AI Tools Used
- **GitHub Copilot**: Used for code completion and scaffolding
- **Assistant**: Used for:
  - Overall pipeline architecture and design
  - Jupyter notebook structure creation
  - Bash script generation for bioinformatics tools
  - Python code for VCF parsing and variant comparison
  - IGV batch script generation
  - GitHub Actions workflow configuration

## AI-Assisted Components

### 1. Notebook Structure
- AI helped design the logical flow of the notebook
- Suggested markdown cells for documentation
- Organized sections for clarity

### 2. Data Download and Processing
- Bash scripts for downloading sequencing data
- Reference genome download and indexing commands
- File management and conditional checks

### 3. Alignment with minimap2
- Technology-appropriate parameters (sr vs map-hifi)
- Samtools commands for sorting and indexing
- Alignment statistics extraction

### 4. Variant Calling with bcftools
- Region-based variant calling for CYP genes
- VCF compression and indexing
- Variant counting and statistics

### 5. Phasing with HapCUT2
- Fragment extraction commands
- Technology-specific parameters for PacBio
- HapCUT2 execution and VCF conversion

### 6. Variant Comparison Analysis
- Python code using pysam to parse VCF files
- Logic for identifying shared and unique variants
- Data structure design for per-gene comparison

### 7. IGV Visualization
- Automated batch script generation
- Screenshot naming conventions
- Manual analysis guidelines

### 8. Star-Allele Determination
- Phased haplotype extraction from VCF
- Star-allele variant definitions from PharmVar
- Matching logic for diplotype determination

### 9. CI/CD Configuration
- GitHub Actions workflow updates
- Dependency installation (minimap2, samtools, bcftools, HapCUT2)
- Jupyter notebook execution command

## Human Contributions
- Domain knowledge about pharmacogenomics
- Understanding of CYP gene clinical importance
- Star-allele significance interpretation
- Quality control decisions
- Final validation and testing

## Notes
- All AI-generated code was reviewed and tested
- Bioinformatics tool parameters were validated against documentation
- Gene coordinates verified using UCSC Genome Browser
- Star-allele definitions cross-referenced with PharmVar database
