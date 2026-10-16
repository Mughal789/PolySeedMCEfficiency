#!/bin/bash
# Build analysis/bin/VtxEffAnalysis (needs only ROOT; no sPHENIX libraries).
# Usage: analysis/build.sh            (run from anywhere, once)
cd "$(dirname "$0")"
# Environment first: the sPHENIX setup script is not compatible with "set -u".
command -v root-config >/dev/null 2>&1 || source /opt/sphenix/core/bin/sphenix_setup.sh -n new.1 >/dev/null 2>&1
set -euo pipefail
mkdir -p bin
# -lEG provides TDatabasePDG (PDG code -> charge).
g++ -O2 -Wall -Wextra -std=c++17 -o bin/VtxEffAnalysis src/VtxEffAnalysis.cc \
    $(root-config --cflags) $(root-config --libs) -lEG
echo "built $(pwd)/bin/VtxEffAnalysis with ROOT $(root-config --version)"
# Step 1 of the tracking-efficiency extension: evaluator skim (EventTree/TruthTree/RecoTree).
g++ -O2 -Wall -Wextra -std=c++17 -o bin/EvalSkim src/EvalSkim.cc \
    $(root-config --cflags) $(root-config --libs) -lEG
echo "built $(pwd)/bin/EvalSkim with ROOT $(root-config --version)"
# Step 2: single-particle tracking efficiency from the skim.
g++ -O2 -Wall -Wextra -std=c++17 -o bin/TrkEffAnalysis src/TrkEffAnalysis.cc \
    $(root-config --cflags) $(root-config --libs) -lEG
echo "built $(pwd)/bin/TrkEffAnalysis with ROOT $(root-config --version)"
