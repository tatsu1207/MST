#!/usr/bin/env bash
# ============================================================================
# MST-Pipeline — Docker entrypoint
# Initializes data directories on first run, then starts the uvicorn server.
# ============================================================================
set -e

# Activate conda
eval "$(conda shell.bash hook)"

# ── Ensure data directories exist (for fresh volume mounts) ──────────────
mkdir -p /app/data/uploads /app/data/datasets

PORT="${PORT:-8050}"

echo "============================================="
echo "  MST-Pipeline"
echo "  http://localhost:${PORT}"
echo "============================================="

exec conda run -n mst --no-capture-output \
    uvicorn app.main:app --host 0.0.0.0 --port "${PORT}"
