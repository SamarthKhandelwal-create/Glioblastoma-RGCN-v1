#!/usr/bin/env pwsh
<#
.SYNOPSIS
    Quick fix + re-run ceRNA inference test
    
.DESCRIPTION
    Runs the test suite with the PyTorch 2.6 fix applied
    
.EXAMPLE
    .\run_test_fixed.ps1 -Verbose
#>

param(
    [switch]$Verbose
)

Write-Host ""
Write-Host "╔════════════════════════════════════════════════════════════════════╗" -ForegroundColor Cyan
Write-Host "║       ceRNA INFERENCE TEST (PyTorch 2.6 Fix Applied)              ║" -ForegroundColor Cyan
Write-Host "╚════════════════════════════════════════════════════════════════════╝" -ForegroundColor Cyan
Write-Host ""

# Set environment
$env:INJECT_MIRNA = 'true'
Write-Host "✓ INJECT_MIRNA=true" -ForegroundColor Green

# Run test with fix
Write-Host ""
Write-Host "Running test suite..." -ForegroundColor Yellow
Write-Host ""

$cmd = "python test_cerna_inference.py --inject"
if ($Verbose) { $cmd += " --verbose" }

Invoke-Expression $cmd
$exitCode = $LASTEXITCODE

Write-Host ""
if ($exitCode -eq 0) {
    Write-Host "✓ SUCCESS: All tests passed!" -ForegroundColor Green
} else {
    Write-Host "⚠ Tests completed with issues (exit code: $exitCode)" -ForegroundColor Yellow
}

Write-Host ""
exit $exitCode
