<#
.SYNOPSIS
  CMM source installer for Windows PowerShell. Clone the repo, then run .\install.ps1

.DESCRIPTION
  Bootstraps uv without relying on a system Python, creates an isolated virtual environment,
  and installs CMM from this checkout (editable). New environments use a uv-managed Python
  3.12 by default. The macOS/Linux equivalent is install.sh.

.EXAMPLE
  .\install.ps1                       # desktop + strain design + gurobipy, into .\.venv
  .\install.ps1 -Dev                  # also install the test/lint tooling (pytest, ruff)
  .\install.ps1 -NoGurobi             # skip gurobipy (open GLPK runs LP/MILP, not QP/MIQP)
  .\install.ps1 -CoreOnly             # omit the desktop and strain-design extras
  .\install.ps1 -Python 3.11          # override the default Python (3.10-3.12 only)
  $env:CMM_PYTHON = 'C:\Python312\python.exe'; .\install.ps1
  .\install.ps1 -VenvDir C:\envs\cmm # install into a chosen venv directory (default: .\.venv)
#>
[CmdletBinding()]
param(
  [switch]$Dev,
  [switch]$NoGurobi,
  [switch]$CoreOnly,
  [string]$Python = "",
  [string]$VenvDir = ".venv"
)

$ErrorActionPreference = "Stop"
$defaultPython = "3.12"
$uvMinimumVersion = "0.8.0"
$uvBootstrapVersion = "0.12.5"

# Run from the repository root (the directory holding this script) regardless of CWD.
Set-Location -Path $PSScriptRoot

$extrasList = @()
if (-not $CoreOnly) { $extrasList += @("desktop", "design") }
if ($Dev)           { $extrasList += "dev" }
if (-not $NoGurobi) { $extrasList += "solver-gurobi" }
$extras = $extrasList -join ","

$pythonRequest = ""
$pythonSource = ""
if ($PSBoundParameters.ContainsKey("Python")) {
  if ([string]::IsNullOrWhiteSpace($Python)) {
    throw "-Python needs an interpreter path, command, or version request."
  }
  $pythonRequest = $Python.Trim()
  $pythonSource = "-Python"
} elseif (-not [string]::IsNullOrWhiteSpace($env:CMM_PYTHON)) {
  $pythonRequest = $env:CMM_PYTHON.Trim()
  $pythonSource = "CMM_PYTHON"
}

function Find-Uv {
  $command = Get-Command uv -ErrorAction SilentlyContinue
  if ($command) {
    return [string]$command.Source
  }

  $candidates = @()
  if (-not [string]::IsNullOrWhiteSpace($env:XDG_BIN_HOME)) {
    $candidates += (Join-Path $env:XDG_BIN_HOME "uv.exe")
  }
  if (-not [string]::IsNullOrWhiteSpace($env:USERPROFILE)) {
    $candidates += (Join-Path $env:USERPROFILE ".local\bin\uv.exe")
    $candidates += (Join-Path $env:USERPROFILE ".cargo\bin\uv.exe")
  }
  foreach ($candidate in $candidates) {
    if (Test-Path -LiteralPath $candidate -PathType Leaf) {
      return $candidate
    }
  }
  return $null
}

function Get-CompatibleUvVersion([string]$Path) {
  $versionOutput = ((& $Path --version 2>&1) -join "").Trim()
  $uvExitCode = $LASTEXITCODE
  if ($uvExitCode -ne 0) {
    throw "uv exists at '$Path' but '$Path --version' failed."
  }
  $versionMatch = [regex]::Match($versionOutput, '^uv\s+(\d+\.\d+\.\d+)(?:\s|$)')
  if (-not $versionMatch.Success) {
    throw "Could not parse uv version from '$versionOutput' ($Path)."
  }
  $installedVersion = $versionMatch.Groups[1].Value
  if ([version]$installedVersion -lt [version]$uvMinimumVersion) {
    throw "CMM source installation requires uv >= $uvMinimumVersion; found $installedVersion at '$Path'. Install the tested bootstrap release with 'irm https://astral.sh/uv/$uvBootstrapVersion/install.ps1 | iex', ensure '$env:USERPROFILE\.local\bin' is first on PATH, and rerun .\install.ps1."
  }
  return $versionOutput
}

function Stop-ManagedPython([string]$Message) {
  throw "$Message CMM supports Python 3.10-3.12. Retry after checking network/platform access, or set `$env:CMM_PYTHON to an existing supported interpreter, for example: `$env:CMM_PYTHON = 'C:\Python312\python.exe'."
}

# --- bootstrap uv before selecting Python -----------------------------------------------
$uvPath = Find-Uv
if ([string]::IsNullOrWhiteSpace($uvPath)) {
  Write-Host "==> uv was not found; installing pinned uv $uvBootstrapVersion"
  $oldNoModifyPath = [Environment]::GetEnvironmentVariable("UV_NO_MODIFY_PATH", "Process")
  try {
    $env:UV_NO_MODIFY_PATH = "1"
    Invoke-RestMethod "https://astral.sh/uv/$uvBootstrapVersion/install.ps1" | Invoke-Expression
  } finally {
    if ($null -eq $oldNoModifyPath) {
      Remove-Item Env:UV_NO_MODIFY_PATH -ErrorAction SilentlyContinue
    } else {
      $env:UV_NO_MODIFY_PATH = $oldNoModifyPath
    }
  }
  $uvPath = Find-Uv
  if ([string]::IsNullOrWhiteSpace($uvPath)) {
    throw "uv installation finished, but the uv executable could not be located. See https://docs.astral.sh/uv/getting-started/installation/."
  }
}
$uvVersionOutput = Get-CompatibleUvVersion $uvPath
Write-Host "==> Using $uvVersionOutput ($uvPath)"

function Find-VenvPython([string]$Directory) {
  $windowsPython = Join-Path $Directory "Scripts\python.exe"
  if (Test-Path -LiteralPath $windowsPython -PathType Leaf) {
    return $windowsPython
  }
  $unixPython = Join-Path $Directory "bin/python"
  if (Test-Path -LiteralPath $unixPython -PathType Leaf) {
    return $unixPython
  }
  return $null
}

function Get-SupportedPythonVersion([string]$Interpreter, [string]$Context) {
  $versionOutput = (& $Interpreter -c "import sys; print('.'.join(map(str, sys.version_info[:3]))); raise SystemExit(0 if (3, 10) <= sys.version_info < (3, 13) else 1)") -join ""
  $pythonExitCode = $LASTEXITCODE
  $versionOutput = $versionOutput.Trim()
  if ($pythonExitCode -ne 0) {
    if ([string]::IsNullOrWhiteSpace($versionOutput)) { $versionOutput = "unknown" }
    throw "$Context uses unsupported Python $versionOutput; CMM supports Python 3.10-3.12."
  }
  return $versionOutput
}

function Get-PythonMinor([string]$Version) {
  $parts = $Version.Split(".")
  return "$($parts[0]).$($parts[1])"
}

function Resolve-PythonRequest([string]$Request, [string]$Source) {
  if ($Request -match '^3\.(10|11|12)(\.\d+)?$') {
    Write-Host "==> Installing requested Python $Request with uv"
    & $uvPath python install $Request | Out-Host
    if ($LASTEXITCODE -ne 0) {
      Stop-ManagedPython "$Source requested Python '$Request', but uv could not install it."
    }
    $resolved = (& $uvPath python find --no-project --managed-python $Request) -join ""
    if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($resolved)) {
      Stop-ManagedPython "uv installed requested Python '$Request' but could not resolve its interpreter."
    }
  } elseif ($Request -match '^\d+\.\d+(\.\d+)?$') {
    throw "$Source requested unsupported Python '$Request'; CMM supports Python 3.10-3.12."
  } else {
    $executable = $null
    if (Test-Path -LiteralPath $Request -PathType Leaf) {
      $executable = (Resolve-Path -LiteralPath $Request).Path
    } else {
      $command = Get-Command $Request -ErrorAction SilentlyContinue
      if ($command) { $executable = [string]$command.Source }
    }
    if ([string]::IsNullOrWhiteSpace($executable)) {
      throw "$Source requested interpreter '$Request', but it was not found; use an executable path or a version such as 3.12."
    }
    $resolved = (& $uvPath python find --no-project --no-python-downloads $executable) -join ""
    if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($resolved)) {
      Stop-ManagedPython "$Source requested interpreter '$Request', but uv could not resolve it."
    }
  }
  return $resolved.Trim()
}

# --- select or install a supported Python -----------------------------------------------
$resolvedPython = ""
$requestedVersion = ""
if (-not [string]::IsNullOrWhiteSpace($pythonRequest)) {
  $resolvedPython = Resolve-PythonRequest $pythonRequest $pythonSource
  $requestedVersion = Get-SupportedPythonVersion $resolvedPython "$pythonSource '$pythonRequest'"
}

if (Test-Path -LiteralPath $VenvDir) {
  if (-not (Test-Path -LiteralPath $VenvDir -PathType Container)) {
    throw "virtual-environment path exists but is not a directory: $VenvDir"
  }
  $venvConfiguration = Join-Path $VenvDir "pyvenv.cfg"
  if (-not (Test-Path -LiteralPath $venvConfiguration -PathType Leaf)) {
    throw "$VenvDir has a Python executable or directory but no pyvenv.cfg, so it is not a reusable virtual environment. Rename it yourself or choose a new path such as -VenvDir .venv312; CMM will not delete or overwrite it."
  }
  $vpy = Find-VenvPython $VenvDir
  if ([string]::IsNullOrWhiteSpace($vpy)) {
    throw "$VenvDir exists but does not contain a Python virtual environment; choose another -VenvDir path or remove it yourself."
  }
  & $vpy -c "import sys; raise SystemExit(0 if sys.prefix != sys.base_prefix else 1)" | Out-Null
  if ($LASTEXITCODE -ne 0) {
    throw "$venvConfiguration exists, but '$vpy' reports sys.prefix == sys.base_prefix and is not running inside a virtual environment. Rename it yourself or choose a new path such as -VenvDir .venv312; CMM will not delete or overwrite it."
  }
  $venvVersion = Get-SupportedPythonVersion $vpy "existing virtual environment '$VenvDir'"
  if (-not [string]::IsNullOrWhiteSpace($resolvedPython)) {
    if ((Get-PythonMinor $venvVersion) -ne (Get-PythonMinor $requestedVersion)) {
      throw "$pythonSource requested Python $(Get-PythonMinor $requestedVersion), but $VenvDir already uses Python $(Get-PythonMinor $venvVersion); choose another -VenvDir path or remove the existing environment yourself."
    }
  }
  Write-Host "==> Reusing existing virtual environment at $VenvDir (Python $venvVersion)"
} else {
  if ([string]::IsNullOrWhiteSpace($resolvedPython)) {
    Write-Host "==> Installing preferred Python $defaultPython with uv"
    & $uvPath python install $defaultPython
    if ($LASTEXITCODE -ne 0) {
      Stop-ManagedPython "uv could not install preferred Python $defaultPython."
    }
    $resolvedPython = (& $uvPath python find --no-project --managed-python $defaultPython) -join ""
    if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($resolvedPython)) {
      Stop-ManagedPython "uv installed Python $defaultPython but could not resolve its interpreter."
    }
    $resolvedPython = $resolvedPython.Trim()
    $requestedVersion = Get-SupportedPythonVersion $resolvedPython "uv-managed Python '$defaultPython'"
  }
  Write-Host "==> Creating virtual environment at $VenvDir with Python $requestedVersion"
  & $uvPath venv --seed --python $resolvedPython $VenvDir
  if ($LASTEXITCODE -ne 0) {
    Stop-ManagedPython "uv could not create the virtual environment at $VenvDir with resolved Python $requestedVersion."
  }
  $vpy = Find-VenvPython $VenvDir
  if ([string]::IsNullOrWhiteSpace($vpy)) {
    throw "uv did not create a usable Python environment at $VenvDir."
  }
}

# --- install ---------------------------------------------------------------------------
if ($extras -ne "") {
  Write-Host "==> Installing CMM (editable) with extras: $extras"
  & $uvPath pip install --python $vpy -e ".[${extras}]"
} else {
  Write-Host "==> Installing CMM (editable, core only)"
  & $uvPath pip install --python $vpy -e .
}
if ($LASTEXITCODE -ne 0) {
  throw "uv could not install CMM into $VenvDir."
}

if ($NoGurobi) {
  Write-Host "==> Skipping gurobipy (GLPK supports LP/MILP; L2 MOMA/E-Flux2 need QP and published MTA/rMTA need MIQP)"
}

# --- done ------------------------------------------------------------------------------
$solver = (& $vpy -c "from cmm.core import solver_status; s = solver_status(); print(s.name, '(' + ', '.join(s.capabilities) + ')')" 2>$null)
if (-not $solver) { $solver = "unknown" }
Write-Host ""
Write-Host "==> CMM installed. Active solver: $solver"
Write-Host ""
Write-Host "   Launch the desktop application:"
Write-Host "     $vpy -m cmm.app"
Write-Host "   Or activate the environment first:"
Write-Host "     $VenvDir\Scripts\Activate.ps1"
Write-Host "     python -m cmm.app"
