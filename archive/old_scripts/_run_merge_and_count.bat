@echo off
set INJECT_MIRNA=true
python scripts\merge_mirna_and_rebuild.py --inject
python scripts\count_graph.py
