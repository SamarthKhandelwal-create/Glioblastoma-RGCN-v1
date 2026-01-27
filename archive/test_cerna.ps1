#!/usr/bin/env pwsh
<#
.SYNOPSIS
    Quick PowerShell script for running ceRNA inference tests
    
.DESCRIPTION
    Sets up environment and runs the complete ceRNA inference test workflow
    with miRNA injection enabled.
    
.PARAMETER Inject
    Enable INJECT_MIRNA mode (default: $true)
    
.PARAMETER Run
    Run full pipeline before testing (default: $true)
    
.PARAMETER Verbose
    Enable verbose output
    
.EXAMPLE
    # Basic test
    .\test_cerna.ps1
    
    # With verbose output
    .\test_cerna.ps1 -Verbose
    
    # Skip pipeline run (use existing graph)
    .\test_cerna.ps1 -Run $false

.NOTES
    Requires Python 3.8+ and pytorch_geometric
#>

param(
    [bool]$Inject = $true,
    [bool]$Run = $true,
    [switch]$Verbose
)

# Color helpers
function Write-Success { Write-Host "✓ $args" -ForegroundColor Green }
function Write-Error { Write-Host "✗ $args" -ForegroundColor Red }
function Write-Warn { Write-Host "⚠ $args" -ForegroundColor Yellow }
function Write-Info { Write-Host "→ $args" -ForegroundColor Cyan }
function Write-Header {
    Write-Host ""
    Write-Host "=" * 70 -ForegroundColor Cyan
    Write-Host "  $args" -ForegroundColor Cyan
    Write-Host "=" * 70 -ForegroundColor Cyan
}

# Main workflow
Write-Header "ceRNA INFERENCE TEST SUITE"

Write-Info "Configuration:"
Write-Info "  INJECT_MIRNA: $Inject"
Write-Info "  Run Pipeline: $Run"
Write-Info "  Verbose: $($PSBoundParameters.ContainsKey('Verbose'))"

# Step 1: Set environment
if ($Inject) {
    $env:INJECT_MIRNA = 'true'
    Write-Success "INJECT_MIRNA enabled"
} else {
    $env:INJECT_MIRNA = 'false'
    Write-Info "INJECT_MIRNA disabled"
}

# Step 2: Run pipeline if requested
if ($Run) {
    Write-Header "RUNNING PIPELINE"
    
    Write-Info "Step 2a: Building node features..."
    python src/preprocessing/build_node_features.py
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Feature building failed"
        exit 1
    }
    Write-Success "Features built"
    
    Write-Info "Step 2b: Building graph with ceRNA inference..."
    python src/graph/build_graph.py
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Graph building failed"
        exit 1
    }
    Write-Success "Graph built with ceRNA edges"
}

# Step 3: Run tests
Write-Header "RUNNING TEST SUITE"

$testCmd = "python test_cerna_inference.py --inject"
if ($Verbose) {
    $testCmd += " --verbose"
}

Write-Info "Executing: $testCmd"
Invoke-Expression $testCmd
$testExitCode = $LASTEXITCODE

# Step 4: Summary
Write-Header "WORKFLOW SUMMARY"

if ($testExitCode -eq 0) {
    Write-Success "All tests passed!"
    Write-Info ""
    Write-Info "Next steps:"
    Write-Info "  1. Review graph statistics above"
    Write-Info "  2. Check results/ directory for outputs"
    Write-Info "  3. Run model training:"
    Write-Info ""
    Write-Info "     python src/training/train_model.py \"
    Write-Info "         --graph-path results/hetero_graph_GBM.pt \"
    Write-Info "         --epochs 100"
    Write-Info ""
    Write-Info "  4. Or run cross-validation:"
    Write-Info ""
    Write-Info "     python src/training/train_cross_validation.py \"
    Write-Info "         --graph-path results/hetero_graph_GBM.pt \"
    Write-Info "         --folds 10 --epochs 100"
} else {
    Write-Error "Tests failed with exit code: $testExitCode"
    Write-Warn ""
    Write-Warn "Troubleshooting:"
    Write-Warn "  1. Check Python environment: python --version"
    Write-Warn "  2. Check dependencies: pip list | grep -E 'torch|pandas|numpy'"
    Write-Warn "  3. Run individual steps to isolate issue:"
    Write-Warn ""
    Write-Warn "     python src/preprocessing/build_node_features.py"
    Write-Warn "     python src/graph/build_graph.py"
    Write-Warn "     python test_cerna_inference.py --inject --verbose"
    exit 1
}

Write-Info "Complete!"
