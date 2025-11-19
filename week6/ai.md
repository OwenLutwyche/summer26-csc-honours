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
I asked Claude to help handle the USA variant format from alevin-fry output. Claude assisted with:
- Building the `collapse_usa()` function to aggregate gene-S, gene-U, gene-A variants to gene-level
- Implementing auto-transpose logic when matrix orientation doesn't match expected format
- Hardening string handling for feature names (converting arrays to strings before manipulation)
- Adding validation checks for cell and gene counts after loading

### 7. Adaptive QC and Normalization
I encountered errors with small datasets during QC. Claude helped me implement:
- QC metrics calculation without `percent_top` parameter (which caused IndexErrors)
- Conservative filtering thresholds: min_genes=10, min_cells=1 for toy data
- Adaptive HVG selection using `n_top_genes=min(2000, n_vars)` to prevent requesting more genes than available
- Fallback logic to use all genes if HVG selection returns zero genes
- Standard normalization workflow (10k counts per cell + log1p transformation)

### 8. Dimensionality Reduction and Clustering
I needed clustering to work with limited features. Claude assisted with:
- Capping PCA components at `n_comps=min(40, n_vars-1)` to avoid exceeding available dimensions
- Setting neighbors parameter to `k=min(10, n_obs-1)` to prevent errors with small cell counts
- Implementing Leiden clustering with igraph backend (flavor="igraph", n_iterations=2, directed=False)
- Generating UMAP embeddings for visualization
- Creating plots colored by cluster assignments

### 9. Gene Symbol Mapping and CellTypist
I wanted to use CellTypist for annotation but needed gene symbols. Claude helped me:
- Parse the t2g_3col.tsv file to extract Ensembl ID to gene symbol mappings
- Implement fallback logic when gene_name is missing (use gene_id instead)
- Aggregate duplicate gene symbols by summing their expression values
- Normalize the symbol-mapped data appropriately for CellTypist (10k counts + log1p)
- Download and configure the Immune_All_High.pkl model
- Add graceful error handling with try-except for cases with insufficient feature overlap

## Problems Encountered and How I Solved Them

### 1. simpleaf Backend Issues
**Problem**: simpleaf defaulted to piscem backend which wasn't available, causing execution failures.
**Solution**: I asked Claude how to force the salmon backend. Claude suggested adding the `--no-piscem` flag to the simpleaf index and quant commands.

### 2. Matrix Orientation Confusion
**Problem**: I wasn't sure if alevin-fry outputs rows=cells or rows=genes, leading to dimension mismatches.
**Solution**: Claude helped me add orientation detection logic that checks the matrix shape against expected barcode/feature counts and auto-transposes if needed.

### 3. USA Variant Handling
**Problem**: USA mode produces 3× features (gene-S, gene-U, gene-A), making the feature space large and harder to interpret.
**Solution**: I asked Claude how to collapse these variants. Claude assisted in implementing the `collapse_usa()` function that sums S/U/A variants back to gene-level expression.

### 4. Small Dataset Dimension Errors
**Problem**: The toy dataset (19 genes, 109 cells) caused errors in PCA, neighbors, and HVG selection when using standard parameters.
**Solution**: Claude helped me identify the issue and implement adaptive parameters throughout: `min(2000, n_vars)` for HVGs, `min(40, n_vars-1)` for PCA, `min(10, n_obs-1)` for neighbors.

### 5. CellTypist Zero Feature Overlap
**Problem**: After gene symbol mapping, only 1 symbol was available versus the 2000+ markers required by the model, causing annotation to fail.
**Solution**: I asked Claude how to handle this gracefully. Claude suggested wrapping the annotation call in try-except to catch ValueError and assign "Unknown" labels with a clear error message for deliverable completeness.

### 6. Leiden Backend Warning
**Problem**: I saw a FutureWarning about leidenalg being deprecated in favor of the igraph backend.
**Solution**: I asked Claude how to update the clustering code. Claude helped me switch to `flavor="igraph"` with appropriate parameters (n_iterations=2, directed=False, random_state=0) to eliminate the warning.

### 7. GitHub Actions CI Failure
**Problem**: The notebook ran successfully on my local machine (WSL), but the GitHub Actions workflow failed during the Reference Indexing cell with errors:
- `could not find 'salmon' in your path: cannot find binary path`
- `Could not open JSON file /home/runner/.alevin_fry/simpleaf_info.json`
The issue was that `jupyter execute` runs the entire notebook without executing individual cells interactively, so the Tool Installation cell (which installs salmon, cargo, alevin-fry, and simpleaf) ran but its binaries weren't available to subsequent cells.

**Solution**: I asked Claude why the CI was failing when it worked locally. Claude identified that the GitHub Actions workflow needed to install the bioinformatics tools before running the notebook, not rely on the notebook's installation cell. Claude assisted me in updating `.github/workflows/actions.yml` to:
- Install Rust/Cargo toolchain via rustup
- Add build dependencies (build-essential, pkg-config, libssl-dev) required for Rust compilation
- Download and install salmon 1.10.0 to `/usr/local/bin/`
- Build alevin-fry and simpleaf from source using `cargo install` with version fallbacks
- Set the `ALEVIN_FRY_HOME` environment variable
- Source `$HOME/.cargo/env` to ensure cargo-installed binaries are in PATH

This ensures all tools are available system-wide before `jupyter execute` runs the notebook, making the CI environment match the local development environment.

