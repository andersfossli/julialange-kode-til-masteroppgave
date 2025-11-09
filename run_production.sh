#!/bin/bash
# Production simulation runner with logging
# This prevents VS Code timeouts and captures all output

# Set up logging
LOG_DIR="logs"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="$LOG_DIR/simulation_$TIMESTAMP.log"

echo "Starting production simulation at $(date)"
echo "Logging to: $LOG_FILE"
echo "==============================================="

# Run Julia with full output logging
julia smr-mcs.jl 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "==============================================="
    echo "Simulation completed successfully at $(date)"
    echo "Output saved to: $LOG_FILE"
else
    echo ""
    echo "==============================================="
    echo "ERROR: Simulation failed with exit code $EXIT_CODE"
    echo "Check log file: $LOG_FILE"
    exit $EXIT_CODE
fi
