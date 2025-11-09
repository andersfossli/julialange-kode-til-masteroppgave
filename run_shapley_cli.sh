#!/bin/bash
# Run Shapley analysis from command line (avoids VS Code crashes)
# Logs output to file to monitor progress

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="shapley_run_$TIMESTAMP.log"

echo "================================================"
echo "SHAPLEY ANALYSIS - COMMAND LINE MODE"
echo "================================================"
echo "Started: $(date)"
echo "Log file: $LOG_FILE"
echo ""
echo "This will take approximately 4-5 hours for 41 reactors"
echo "You can monitor progress in another terminal with:"
echo "  tail -f $LOG_FILE"
echo ""
echo "Press Ctrl+C to cancel, or wait 10 seconds to continue..."
sleep 10

echo "Starting Shapley analysis..."
julia run_3_shapley.jl 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Shapley analysis completed at $(date)"
    echo "Results saved to _output/shapley-*.csv"
else
    echo "ERROR: Shapley analysis failed with exit code $EXIT_CODE"
    echo "Check log file: $LOG_FILE"
    exit $EXIT_CODE
fi
echo "================================================"
