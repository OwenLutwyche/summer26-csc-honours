#!/bin/bash
set -e

# Always run relative operations from the script directory so paths resolve
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "[INFO] Debug Build Started (Manual Linking + PIC + Runtime Support)..."

# 1. Setup Python Environment
PYTHON_LIB=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PY_LIB_NAME="python${PYTHON_VERSION}"

# 2. Setup Codon Environment (CRITICAL FIX)
# We need to find where libcodonrt.so lives to fix 'undefined symbol: seq_personality'
CODON_LIB_PATH="$HOME/.codon/lib/codon"

if [ ! -d "$CODON_LIB_PATH" ]; then
    # Fallback search if standard path fails
    CODON_BIN=$(which codon)
    if [ -n "$CODON_BIN" ]; then
        CODON_LIB_PATH=$(dirname $(dirname "$CODON_BIN"))/lib/codon
    fi
fi

echo "Python Lib: $PYTHON_LIB (-l$PY_LIB_NAME)"
echo "Codon Lib:  $CODON_LIB_PATH (-lcodonrt)"

if [ ! -f "$CODON_LIB_PATH/libcodonrt.so" ]; then
    echo "[ERROR] Could not find libcodonrt.so at $CODON_LIB_PATH"
    echo "        Please check your Codon installation."
    exit 1
fi

# Export for Codon compilation
if [ -f "$PYTHON_LIB/lib${PY_LIB_NAME}.so" ]; then
    export CODON_PYTHON="$PYTHON_LIB/lib${PY_LIB_NAME}.so"
else
    export CODON_PYTHON="$PYTHON_LIB/libpython3.so"
fi

# 3. Cleanup
rm -f scancodon_native.so
rm -f src/scanpy/__init__.o

# 4. Compile to Object File
echo "Compiling to Object Code (with PIC)..."
codon build \
    -release \
    -pyext \
    --relocation-model=pic \
    -module scancodon_native \
    src/scanpy/__init__.codon

if [ ! -f "src/scanpy/__init__.o" ]; then
    echo "[ERROR] Compilation failed."
    exit 1
fi

echo "[OK] Object file created."

# 5. MANUAL LINKING
# Added: -L... -lcodonrt (Links Codon Runtime)
# Added: -Wl,-rpath... (Tells the .so where to find libcodonrt at runtime)
echo "Linking Shared Object..."
gcc -shared -fPIC \
    -o scancodon_native.so \
    src/scanpy/__init__.o \
    -L"$PYTHON_LIB" -l"$PY_LIB_NAME" \
    -L"$CODON_LIB_PATH" -lcodonrt \
    -Wl,-rpath,"$CODON_LIB_PATH"

# 6. Verify
file_type=$(file scancodon_native.so)
echo "File Check: $file_type"

if [[ "$file_type" == *"shared object"* ]]; then
    echo "[SUCCESS] Library is ready."
    echo "Test command: python3 -c 'import scancodon_native; print(scancodon_native)'"
else
    echo "[ERROR] Manual linking produced invalid file."
    exit 1
fi