#!/bin/bash
set -e

echo "🚀 Building Scancodon Python Extension..."

# 1. CLEANUP (Crucial step to avoid false positives)
rm -f scancodon_native.so
rm -f src/scanpy/__init__.so

# 2. Configure Paths
PYTHON_LIB=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
PYTHON_INC=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('INCLUDEPY'))")
PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")

# Logic to find the .so file
if [ -f "$PYTHON_LIB/libpython${PYTHON_VERSION}.so" ]; then
    export CODON_PYTHON="$PYTHON_LIB/libpython${PYTHON_VERSION}.so"
else
    # Fallback for systems where LDLIBRARY isn't exact
    export CODON_PYTHON="$PYTHON_LIB/libpython3.so"
fi
export CPATH="$PYTHON_INC:$CPATH"

echo "Linker Path: $CODON_PYTHON"

# 3. Compile
# We DO NOT use -o here. We let Codon build to the default location to ensure linking happens.
echo "Compiling..."
codon build \
    -release \
    -pyext \
    -module scancodon_native \
    src/scanpy/__init__.codon

# 4. Move the output
# Codon defaults to naming the output after the input file
if [ -f "src/scanpy/__init__.so" ]; then
    echo "Found output at src/scanpy/__init__.so. Moving to project root..."
    mv src/scanpy/__init__.so scancodon_native.so
else
    echo "❌ Compile failed: Expected output file 'src/scanpy/__init__.so' not found."
    exit 1
fi

# 5. VERIFY (The actual check)
file_type=$(file scancodon_native.so)
echo "File Check: $file_type"

if [[ "$file_type" == *"shared object"* ]]; then
    echo "✅ SUCCESS: Valid Shared Object created."
    echo "To test: python3 -c 'import scancodon_native; print(dir(scancodon_native))'"
else
    echo "❌ FAILURE: File is not a shared object. Python cannot load it."
    exit 1
fi