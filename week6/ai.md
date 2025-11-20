# AI Usage Documentation - Week 6

## AI Model Used
- Claude Sonnet 4.5

## Components Developed with AI Assistance

### 1. System Requirements Check
I needed to verify that all dependencies were available before running the pipeline. Claude helped me create a comprehensive environment check cell that:
- Checks for system packages (cmake, libtbb12, wget, curl) and prints installation commands if missing
- Embeds a Python script within bash to verify and auto-install pip packages (scanpy, anndata<0.11, leidenalg, celltypist)
- Handles the anndata version constraint specifically for Python 3.13 compatibility
- Avoids running sudo commands directly in the notebook (WSL environment constraint)
- Provides clear feedback about what's installed vs what needs manual installation

### 2. Data Download and Setup
I asked Claude to help with the data acquisition workflow since the dataset was hosted on Box. Claude designed the data setup cell to:
- Check if data already exists at `data/toy_ref_read/` to avoid redundant downloads
- Display Box download instructions only when data is missing (with clear manual steps)
- Automatically download the 10x Chromium v3 whitelist from GitHub
- Normalize the data directory structure to a canonical `toy_ref_read` layout regardless of extraction naming
- Clean up the downloaded tarball after successful extraction
- Verify that all required files (genome.fa, genes.gtf, FASTQs) are present

### 3. Tool Installation
I needed to install salmon, Rust/cargo, alevin-fry, and simpleaf within the notebook. Claude assisted with:
- Downloading and installing salmon 1.10.0 from GitHub releases to `~/.local/bin/`
- Checking for cargo availability and auto-installing Rust toolchain via rustup if missing
- Keeping the Rust toolchain updated with `rustup update`
- Implementing version fallback logic for simpleaf (trying 0.18.4 → 0.18.3 → 0.18.2 → 0.18.1 → 0.17.3) to work around dependency issues
- Running cargo installations with visible output for transparency
- Checking for git availability (required by cargo for fetching crates)
- Setting ALEVIN_FRY_HOME environment variable for simpleaf configuration

### 4. Pipeline Structure
I asked Claude to help design the scRNA-seq analysis pipeline following Single-cell Best Practices 3.8.2. Claude assisted with:
- Structuring the workflow: USA-mode quantification with alevin-fry via simpleaf
- Organizing QC filtering, normalization, dimensionality reduction, and clustering steps
- Planning the integration of Leiden clustering and CellTypist annotation
- Suggesting adaptive parameters to handle small toy datasets

### 5. USA-Mode Quantification
I needed to implement quantification with spliced/unspliced read tracking. Claude helped me:
- Understand the USA mode (`-u` flag) and its purpose for RNA velocity analysis
- Configure simpleaf to use salmon backend via `--no-piscem` flag
- Set up t2g_3col.tsv mapping for proper transcript-to-gene resolution
- Aggregate multiple FASTQ read pairs in a single quantification command
- Add matrix orientation detection logic (rows=cells vs rows=genes)

### 6. Data Loading and USA Collapse
I asked Claude to help handle the USA variant format from alevin-fry output. With Claude's help, I:
- Built the `collapse_usa()` function to aggregate gene-S, gene-U, gene-A variants to gene-level
- Implemented auto-transpose logic when matrix orientation doesn't match expected format
- Hardened string handling for feature names (converting arrays to strings before manipulation)
- Added validation checks for cell and gene counts after loading

### 7. Adaptive QC and Normalization
I encountered errors with small datasets during QC. Claude helped me implement:
- QC metrics calculation without `percent_top` parameter (which caused IndexErrors)
- Conservative filtering thresholds: min_genes=10, min_cells=1 for toy data
- Adaptive HVG selection using `n_top_genes=min(2000, n_vars)` to prevent requesting more genes than available
- Fallback logic to use all genes if HVG selection returns zero genes
- Standard normalization workflow (10k counts per cell + log1p transformation)

### 8. Dimensionality Reduction and Clustering
I needed clustering to work with limited features. With Claude's assistance, I implemented:
- Capping PCA components at `n_comps=min(40, n_vars-1)` to avoid exceeding available dimensions
- Setting neighbors parameter to `k=min(10, n_obs-1)` to prevent errors with small cell counts
- Leiden clustering with igraph backend (flavor="igraph", n_iterations=2, directed=False)
- UMAP embeddings for visualization
- Plots colored by cluster assignments

### 9. Gene Symbol Mapping and CellTypist
I wanted to use CellTypist for annotation but needed gene symbols. Claude helped me:
- Parse the t2g_3col.tsv file to extract Ensembl ID to gene symbol mappings
- Implement fallback logic when gene_name is missing (use gene_id instead)
- Aggregate duplicate gene symbols by summing their expression values
- Normalize the symbol-mapped data appropriately for CellTypist (10k counts + log1p)
- Download and configure the Immune_All_High.pkl model
- Add graceful error handling with try-except for cases with insufficient feature overlap

## Problems Encountered and How I Solved Them

### Initial Implementation (Python 3.13)

### 1. simpleaf Backend Issues
**Problem**: simpleaf defaulted to piscem backend which wasn't available, causing execution failures.
**Solution**: I asked Claude how to force the salmon backend. Claude suggested adding the `--no-piscem` flag to the simpleaf index and quant commands.

### 2. Matrix Orientation Confusion
**Problem**: I wasn't sure if alevin-fry outputs rows=cells or rows=genes, leading to dimension mismatches.
**Solution**: Claude helped me add orientation detection logic that checks the matrix shape against expected barcode/feature counts and auto-transposes if needed.

### 3. USA Variant Handling
**Problem**: USA mode produces 3× features (gene-S, gene-U, gene-A), making the feature space large and harder to interpret.
**Solution**: I asked Claude how to collapse these variants. With Claude's help, I implemented the `collapse_usa()` function that sums S/U/A variants back to gene-level expression.

### 4. Small Dataset Dimension Errors
**Problem**: The toy dataset (19 genes, 109 cells) caused errors in PCA, neighbors, and HVG selection when using standard parameters.
**Solution**: Claude helped me identify the issue and implement adaptive parameters throughout: `min(2000, n_vars)` for HVGs, `min(40, n_vars-1)` for PCA, `min(10, n_obs-1)` for neighbors.

### 5. CellTypist Zero Feature Overlap
**Problem**: After gene symbol mapping, only 1 symbol was available versus the 2000+ markers required by the model, causing annotation to fail.
**Solution**: I asked Claude how to handle this gracefully. Claude suggested wrapping the annotation call in try-except to catch ValueError and assign "Unknown" labels with a clear error message for deliverable completeness.


### 6. Python 3.13 Compatibility Issues
**Problem**: I initially tried to implement the notebook using Python 3.13 to maintain consistency with the rest of the course. I attempted to work around pyroe (which doesn't support Python 3.13) by manually loading the alevin-fry output using anndata, but I couldn't get satisfactory results. The pyroe package provides essential functionality for properly handling USA-mode quantification data that's difficult to replicate manually.

**Solution**: I decided to switch to Python 3.12 to use pyroe properly. With Claude's assistance, I updated the notebook documentation to clearly specify Python 3.12 as a hard requirement, with pyroe's incompatibility being the primary reason. I also configured the GitHub Actions workflow to use dual Python environments, installing both versions and splitting package installations so weeks 1-5 use Python 3.13 while week 6 explicitly uses Python 3.12.

### Final Implementation (Python 3.12)

### 7. GitHub Actions CI Failure
**Problem**: The notebook executed successfully on my local machine but failed in GitHub Actions during the Reference Indexing cell. The issue was that `jupyter execute` runs the entire notebook non-interactively, so the Tool Installation cell's installations weren't available to PATH in subsequent cells.

**Solution**: With Claude's help, I updated `.github/workflows/actions.yml` to pre-install the bioinformatics tools system-wide before running the notebook. This involved installing Rust/Cargo toolchain, downloading salmon 1.10.0, building alevin-fry and simpleaf from source with version fallbacks, and properly setting environment variables like ALEVIN_FRY_HOME. The workflow also needed to handle dual Python environments (3.13 for weeks 1-5, 3.12 for week 6).

### 8. Analysis Summary Enhancement
**Problem**: The initial analysis summary was sparse and didn't provide enough insight into the results.
**Solution**: Claude helped me restructure the summary cell with detailed statistics covering dataset overview, clustering results, cell type distribution, and quality metrics.

