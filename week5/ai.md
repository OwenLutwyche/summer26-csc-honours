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
- Fragment extraction commands with extractHAIRS
- Technology-specific parameters for PacBio (--pacbio 1 requires --ref)
- HapCUT2 execution and VCF output handling
- File format conversions and naming conventions

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

## Debugging and Problem Resolution

### Issue 1: extractHAIRS PacBio Parameter Requirements
**Problem**: Cell 5.1 failed with error "In order to realign variants (including --pacbio and --ont options), reference fasta file must be provided with --ref option"

**Solution**: 
- Added `--ref data/chr10.fa` parameter to extractHAIRS command for PacBio data
- PacBio long reads require local realignment for accurate fragment extraction
- Illumina short reads don't need this parameter

**AI Assistance**: Analyzed error message, researched HapCUT2 documentation, and identified the missing parameter requirement

### Issue 2: HAPCUT2 VCF Format Requirements
**Problem**: Cell 5.2 failed with "VCF file has no header, ERROR reading VCF file: less than 10 columns"

**Solution**:
- Changed input VCF from compressed (.vcf.gz) to uncompressed (.vcf)
- Added decompression step: `gunzip -c results/illumina.vcf.gz > results/illumina.vcf`
- HAPCUT2 and extractHAIRS require uncompressed VCF files

**AI Assistance**: Diagnosed the issue by recognizing HAPCUT2's VCF format requirements

### Issue 3: HAPCUT2 Output File Naming
**Problem**: Cell 6 failed with FileNotFoundError for `results/illumina.phased.vcf` - file existed as `results/illumina.phased.phased.VCF` (uppercase)

**Solution**:
- Discovered HAPCUT2 appends `.phased.VCF` (uppercase) to the output basename
- When using `--output results/illumina.phased`, it creates `results/illumina.phased.phased.VCF`
- Added rename commands: `mv results/illumina.phased.phased.VCF results/illumina.phased.vcf`
- Updated file existence checks to look for `.phased.phased.VCF` before renaming

**AI Assistance**: Investigated HAPCUT2 output behavior and implemented file renaming logic

### Issue 4: Understanding Empty Results for CYP2C9 and CYP2C8
**Problem**: Variant comparison showed 0 variants for CYP2C9 and CYP2C8, only CYP2C19 had results

**Investigation**:
- Used `samtools view -c` to check read coverage in each gene region
- Found: CYP2C19 has excellent coverage, CYP2C9 has 1 read, CYP2C8 has 0 reads
- Determined this is targeted sequencing data focused on CYP2C19 only

**Documentation**:
- Added "Data Coverage Observations" section to conclusions
- Explained that empty results are expected for regions without sequencing coverage
- Clarified this is typical of targeted sequencing panels

**AI Assistance**: Designed coverage analysis approach, ran diagnostic commands, and documented findings

## Notes
- All AI-generated code was reviewed and tested
- Bioinformatics tool parameters were validated against documentation
- Gene coordinates verified using UCSC Genome Browser
- Star-allele definitions cross-referenced with PharmVar database
- Pipeline iteratively debugged to handle tool-specific requirements
- File format and naming conventions discovered through testing
