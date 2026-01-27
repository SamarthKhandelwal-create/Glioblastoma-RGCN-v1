#!/usr/bin/env pwsh
<#
.SYNOPSIS
    Execute the exact workflow: $env:INJECT_MIRNA='true' + python test_cerna_inference.py
    
.DESCRIPTION
    This script reproduces the exact command sequence from the user request:
    $env:INJECT_MIRNA='true'
    python test_cerna_inference.py --inject
    
.EXAMPLE
    .\run_exact_test.ps1
    
.NOTES
    Shorthand for the full workflow mentioned in the issue
#>

param(
    [switch]$Run,
    [switch]$Verbose
)

Write-Host ""
Write-Host "╔════════════════════════════════════════════════════════════════════╗" -ForegroundColor Cyan
Write-Host "║           ceRNA INFERENCE TEST - EXACT WORKFLOW                    ║" -ForegroundColor Cyan
Write-Host "╚════════════════════════════════════════════════════════════════════╝" -ForegroundColor Cyan
Write-Host ""

# Set INJECT_MIRNA
Write-Host "Setting: " -NoNewline
Write-Host "`$env:INJECT_MIRNA='true'" -ForegroundColor Green
$env:INJECT_MIRNA = 'true'
Write-Host ""

# Optional: Run pipeline first
if ($Run) {
    Write-Host "Building pipeline..." -ForegroundColor Yellow
    Write-Host ""
    
    Write-Host "Step 1: Building node features..." -ForegroundColor Cyan
    python src/preprocessing/build_node_features.py
    Write-Host ""
    
    Write-Host "Step 2: Building graph with ceRNA inference..." -ForegroundColor Cyan
    python src/graph/build_graph.py
    Write-Host ""
}

# Run the exact command
Write-Host "Executing test command:" -ForegroundColor Yellow
$cmd = "python test_cerna_inference.py --inject"
if ($Verbose) { $cmd += " --verbose" }
Write-Host $cmd -ForegroundColor Green
Write-Host ""

Invoke-Expression $cmd
$exitCode = $LASTEXITCODE

Write-Host ""
Write-Host "╔════════════════════════════════════════════════════════════════════╗" -ForegroundColor Cyan
Write-Host "║                        TEST EXECUTION COMPLETE                     ║" -ForegroundColor Cyan
Write-Host "╚════════════════════════════════════════════════════════════════════╝" -ForegroundColor Cyan
Write-Host ""

if ($exitCode -eq 0) {
    Write-Host "✓ SUCCESS: All ceRNA inference tests passed" -ForegroundColor Green
    Write-Host "✓ Graph ready for model training" -ForegroundColor Green
} else {
    Write-Host "✗ FAILED: Exit code $exitCode" -ForegroundColor Red
}

Write-Host ""
exit $exitCode
