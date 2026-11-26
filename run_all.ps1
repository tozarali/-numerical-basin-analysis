# One-shot runner for Windows + Miniconda/Anaconda.
#   run_all.cmd                 -> all exploratory experiments (resumable)
#   run_all.cmd --quick         -> smoke versions (a few minutes)
#   run_all.cmd --list          -> registry and rough durations
#   run_all.cmd --only regime_map,diagnostic_eval
#   run_all.cmd --confirmatory  -> also the LOCKED confirmatory test
# Any arguments are passed to `python -m basinlab.run_all`.
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
$runnerArgs = $args

$conda = $null
$command = Get-Command conda.exe -ErrorAction SilentlyContinue
if ($command) {
    $conda = $command.Source
} elseif ($env:CONDA_EXE -and (Test-Path -LiteralPath $env:CONDA_EXE -PathType Leaf)) {
    $conda = $env:CONDA_EXE
} else {
    $candidates = @(
        (Join-Path $env:USERPROFILE 'miniconda3\Scripts\conda.exe'),
        (Join-Path $env:LOCALAPPDATA 'miniconda3\Scripts\conda.exe'),
        'C:\ProgramData\miniconda3\Scripts\conda.exe',
        (Join-Path $env:USERPROFILE 'anaconda3\Scripts\conda.exe'),
        'C:\ProgramData\anaconda3\Scripts\conda.exe'
    )
    foreach ($candidate in $candidates) {
        if (Test-Path -LiteralPath $candidate -PathType Leaf) { $conda = $candidate; break }
    }
}
if (-not $conda) {
    throw 'Miniconda not found. Install Miniconda or open a Miniconda terminal and retry.'
}
Write-Host "Conda: $conda"

$environmentJson = & $conda env list --json
if ($LASTEXITCODE -ne 0) { throw 'Cannot list Conda environments.' }
$environments = ($environmentJson | ConvertFrom-Json).envs
$matches = @($environments | Where-Object { (Split-Path -Leaf $_) -eq 'basin-project' })
if ($matches.Count -gt 1) {
    throw 'Multiple basin-project environments exist; resolve the ambiguity before continuing.'
}
if ($matches.Count -eq 0) {
    Write-Host 'Creating basin-project environment from conda-forge (see environment.yml)...'
    & $conda create --yes --name basin-project --override-channels -c conda-forge 'python=3.12' numpy scipy matplotlib packaging psutil
    if ($LASTEXITCODE -ne 0) { throw 'Conda environment creation failed; see output above.' }
} else {
    Write-Host "Using existing environment: $($matches[0])"
}

& $conda run --no-capture-output -n basin-project python -m basinlab.experiments.check_environment
if ($LASTEXITCODE -ne 0) {
    throw 'Environment dependency check failed. Fix the existing environment before retrying.'
}

& $conda run --no-capture-output -n basin-project python -m unittest discover -s tests
if ($LASTEXITCODE -ne 0) { throw 'Tests failed; experiments not started.' }

& $conda run --no-capture-output -n basin-project python -m basinlab.run_all @runnerArgs
$code = $LASTEXITCODE
if ($code -ne 0) { throw "Some experiments failed or were blocked (exit $code); see data/_runs/LATEST.txt." }
Write-Host 'All selected experiments finished. Report: see data/_runs/LATEST.txt'
