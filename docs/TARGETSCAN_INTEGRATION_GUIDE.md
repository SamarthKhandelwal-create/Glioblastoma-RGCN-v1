"""
Guide: Downloading and Integrating Real TargetScan Predictions

TargetScan provides comprehensive miRNA target predictions. Follow these steps:

1. DOWNLOAD TARGETSCAN DATA
   - Go to: http://www.targetscan.org/
   - Select organism (e.g., "Human", "Homo sapiens")
   - Download: "Predicted Targets Info" CSV
   - Expected columns: miRNA, Gene Symbol, Cumulative Score (or similar)
   - Save to: results/targetscan_raw.csv

2. PREPARE FOR INTEGRATION
   The script expects columns: mirna, gene_symbol, score
   
   Standard TargetScan CSV format typically has:
   - Column 1: miRNA (e.g., "hsa-miR-1")
   - Column 2: Gene Symbol (e.g., "TP53")
   - Column 3: Cumulative Score (0-1, higher = more likely target)

3. INTEGRATE WITH GRAPH BUILDER
   python src/preprocessing/integrate_targetscan.py --input results/targetscan_raw.csv --score-cutoff 0.5

4. REBUILD GRAPH
   $env:INJECT_MIRNA='true'
   python src/graph/build_graph.py

EXAMPLE DATA STRUCTURE:
miRNA,Gene Symbol,Cumulative Score
hsa-miR-21,TP53,0.92
hsa-miR-21,PTEN,0.85
hsa-miR-155,FOXO3,0.78
hsa-miR-1,EGFR,0.65
...

Once you have real TargetScan data, place it in results/ and run the integration script above.
"""
