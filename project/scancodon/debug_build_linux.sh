#!/bin/bash
set -e

# Always run relative operations from the script directory so paths resolve
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "[INFO] Debug Build Started (Linux)..."

# 1. Setup Python Environment
PYTHON_LIB=$(python3 -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
PY_LIB_NAME="python${PYTHON_VERSION}"

# 2. Setup Codon Environment
CODON_BIN=$(which codon)
CODON_LIB_PATH="$(dirname "$(dirname "$CODON_BIN")")/lib/codon"




if [ ! -d "$CODON_LIB_PATH" ]; then
    CODON_BIN=$(which codon)
    if [ -n "$CODON_BIN" ]; then
        CODON_LIB_PATH="$(dirname "$(dirname "$CODON_BIN")")/lib/codon"
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
elif [ -f "$PYTHON_LIB/lib${PY_LIB_NAME}.so.1.0" ]; then
    export CODON_PYTHON="$PYTHON_LIB/lib${PY_LIB_NAME}.so.1.0"
else
    echo "[ERROR] Could not find Python shared library."
    exit 1
fi

# 3. Cleanup
rm -f scancodon_native.so
rm -f src/scanpy/__init__.o

# 4. Compile to Object File
echo "Compiling to Object Code (with PIC)..."
codon build \
    -release \
    -pyext \
    -module scancodon_native \
    --relocation-model=pic \
    src/scanpy/__init__.codon

if [ ! -f "src/scanpy/__init__.o" ]; then
    echo "[ERROR] Compilation failed."
    exit 1
fi

echo "[OK] Object file created."

# 5. Link Shared Object
echo "Linking Shared Object..."
gcc -shared -fPIC \
    -o scancodon_native.so \
    src/scanpy/__init__.o \
    -L"$PYTHON_LIB" -l"$PY_LIB_NAME" \
    -L"$CODON_LIB_PATH" -lcodonrt \
    -Wl,-rpath,"$CODON_LIB_PATH"

# 6. Verify
if [ ! -f "scancodon_native.so" ]; then
    echo "[ERROR] Failed to create shared library."
    exit 1
fi

