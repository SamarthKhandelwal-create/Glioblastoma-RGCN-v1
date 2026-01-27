# Helper script to merge miRNA features (if present) and run graph builder
# Usage: python scripts/merge_mirna_and_rebuild.py [--inject]

"""Helper: ensure miRNA features exist and run graph builder.
Usage:
    python scripts/merge_mirna_and_rebuild.py [--inject]
"""
import subprocess
import os
import sys
from pathlib import Path

inject = '--inject' in sys.argv
if inject:
    os.environ['INJECT_MIRNA'] = 'true'

# Run build_graph
rc = subprocess.call([sys.executable, 'src/graph/build_graph.py'])
if rc != 0:
    print('build_graph.py failed with code', rc)
    sys.exit(rc)
print('Graph rebuild complete.')
