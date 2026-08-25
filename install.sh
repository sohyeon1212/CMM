#!/usr/bin/env bash
#
# CMM source installer — clone the repo, then run ./install.sh
#
# Bootstraps uv without relying on a system Python, creates an isolated virtual environment,
# and installs CMM from this checkout (editable). New environments use a uv-managed Python
# 3.12 by default. Works on macOS, Linux, and Windows (Git Bash / WSL); for native PowerShell
# use install.ps1.
#
# Usage:
#   ./install.sh                      # desktop + strain design + gurobipy, into ./.venv
#   ./install.sh --dev                # also install the test/lint tooling (pytest, ruff)
#   ./install.sh --no-gurobi          # skip gurobipy (open GLPK runs LP/MILP, not QP/MIQP)
#   ./install.sh --core-only          # omit the desktop and strain-design extras
#   ./install.sh --python 3.11        # override the default Python (3.10-3.12 only)
#   CMM_PYTHON=/path/to/python ./install.sh
#   ./install.sh --venv /path/to/venv # install into a chosen venv directory (default: ./.venv)
#
set -euo pipefail

# Run from the repository root (the directory holding this script) regardless of CWD.
cd "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

WITH_DEV=0
WITH_GUROBI=1
CORE_ONLY=0
PYTHON_REQUEST="${CMM_PYTHON:-}"
PYTHON_SOURCE=""
VENV_DIR=".venv"
DEFAULT_PYTHON="3.12"
UV_MINIMUM_VERSION="0.8.0"
UV_BOOTSTRAP_VERSION="0.12.5"

if [ -n "$PYTHON_REQUEST" ]; then
  PYTHON_SOURCE="CMM_PYTHON"
fi

usage() {
  sed -n '2,18p' "$0"
}

while [ $# -gt 0 ]; do
  case "$1" in
    --dev)
      WITH_DEV=1
      ;;
    --no-gurobi)
      WITH_GUROBI=0
      ;;
    --core-only)
      CORE_ONLY=1
      ;;
    --python)
      PYTHON_REQUEST="${2:?--python needs an interpreter or version request}"
      PYTHON_SOURCE="--python"
      shift
      ;;
    --venv)
      VENV_DIR="${2:?--venv needs a path}"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "unknown option: $1 (try --help)" >&2
      exit 2
      ;;
  esac
  shift
done

append_extra() {
  if [ -n "${EXTRAS:-}" ]; then
    EXTRAS="$EXTRAS,$1"
  else
    EXTRAS="$1"
  fi
}

EXTRAS=""
if [ "$CORE_ONLY" -eq 0 ]; then
  append_extra "desktop"
  append_extra "design"
fi
if [ "$WITH_DEV" -eq 1 ]; then
  append_extra "dev"
fi
if [ "$WITH_GUROBI" -eq 1 ]; then
  append_extra "solver-gurobi"
fi

find_uv_after_bootstrap() {
  if command -v uv >/dev/null 2>&1; then
    command -v uv
    return 0
  fi
  if [ -n "${XDG_BIN_HOME:-}" ] && [ -x "$XDG_BIN_HOME/uv" ]; then
    printf '%s\n' "$XDG_BIN_HOME/uv"
    return 0
  fi
  if [ -n "${HOME:-}" ]; then
    for candidate in "$HOME/.local/bin/uv" "$HOME/.cargo/bin/uv"; do
      if [ -x "$candidate" ]; then
        printf '%s\n' "$candidate"
        return 0
      fi
    done
  fi
  return 1
}

version_at_least() {
  local actual_major actual_minor actual_patch
  local minimum_major minimum_minor minimum_patch
  IFS=. read -r actual_major actual_minor actual_patch <<< "$1"
  IFS=. read -r minimum_major minimum_minor minimum_patch <<< "$2"
  if [ "$actual_major" -ne "$minimum_major" ]; then
    [ "$actual_major" -gt "$minimum_major" ]
  elif [ "$actual_minor" -ne "$minimum_minor" ]; then
    [ "$actual_minor" -gt "$minimum_minor" ]
  else
    [ "$actual_patch" -ge "$minimum_patch" ]
  fi
}

require_compatible_uv() {
  local version_output
  local installed_version
  if ! version_output="$("$UV_BIN" --version 2>&1)"; then
    echo "error: uv exists at '$UV_BIN' but '$UV_BIN --version' failed." >&2
    exit 1
  fi
  if [[ "$version_output" =~ ^uv[[:space:]]+([0-9]+\.[0-9]+\.[0-9]+)([[:space:]]|$) ]]; then
    installed_version="${BASH_REMATCH[1]}"
  else
    echo "error: could not parse uv version from '$version_output' ($UV_BIN)." >&2
    exit 1
  fi
  if ! version_at_least "$installed_version" "$UV_MINIMUM_VERSION"; then
    echo "error: CMM source installation requires uv >= $UV_MINIMUM_VERSION; found $installed_version at '$UV_BIN'." >&2
    echo "Install the tested bootstrap release and make it first on PATH, then rerun ./install.sh:" >&2
    echo "  curl -LsSf https://astral.sh/uv/$UV_BOOTSTRAP_VERSION/install.sh | sh" >&2
    echo '  [ -f "$HOME/.local/bin/env" ] && . "$HOME/.local/bin/env"' >&2
    echo '  export PATH="$HOME/.local/bin:$PATH"' >&2
    exit 1
  fi
  printf '%s\n' "$version_output"
}

managed_python_error() {
  echo "error: $1" >&2
  echo "CMM supports Python 3.10-3.12. Retry after checking network/platform access, or point CMM_PYTHON at an existing supported interpreter:" >&2
  echo "  CMM_PYTHON=/path/to/python3.12 ./install.sh" >&2
  exit 1
}

# --- bootstrap uv before selecting Python -----------------------------------------------
if ! UV_BIN="$(find_uv_after_bootstrap)"; then
  echo "==> uv was not found; installing pinned uv $UV_BOOTSTRAP_VERSION"
  if command -v curl >/dev/null 2>&1; then
    curl -LsSf "https://astral.sh/uv/$UV_BOOTSTRAP_VERSION/install.sh" | env UV_NO_MODIFY_PATH=1 sh
  elif command -v wget >/dev/null 2>&1; then
    wget -qO- "https://astral.sh/uv/$UV_BOOTSTRAP_VERSION/install.sh" | env UV_NO_MODIFY_PATH=1 sh
  else
    echo "error: installing uv requires curl or wget; see https://docs.astral.sh/uv/getting-started/installation/." >&2
    exit 1
  fi
  if ! UV_BIN="$(find_uv_after_bootstrap)"; then
    echo "error: uv installation finished, but the uv executable could not be located." >&2
    exit 1
  fi
fi
UV_VERSION_OUTPUT="$(require_compatible_uv)"
echo "==> Using $UV_VERSION_OUTPUT ($UV_BIN)"

venv_python() {
  if [ -x "$VENV_DIR/bin/python" ]; then
    printf '%s\n' "$VENV_DIR/bin/python"
  elif [ -x "$VENV_DIR/Scripts/python.exe" ]; then
    printf '%s\n' "$VENV_DIR/Scripts/python.exe"
  else
    return 1
  fi
}

supported_python_version() {
  "$1" -c 'import sys; print(".".join(map(str, sys.version_info[:3]))); raise SystemExit(0 if (3, 10) <= sys.version_info < (3, 13) else 1)'
}

python_minor() {
  printf '%s\n' "$1" | awk -F. '{print $1 "." $2}'
}

resolve_python_request() {
  local request="$1"
  local resolved
  if [[ "$request" =~ ^3\.(10|11|12)(\.[0-9]+)?$ ]]; then
    echo "==> Installing requested Python $request with uv" >&2
    if ! "$UV_BIN" python install "$request" >&2; then
      managed_python_error "$PYTHON_SOURCE requested Python '$request', but uv could not install it."
    fi
    if ! resolved="$("$UV_BIN" python find --no-project --managed-python "$request")" || [ -z "$resolved" ]; then
      managed_python_error "uv installed requested Python '$request' but could not resolve its interpreter."
    fi
  elif [[ "$request" =~ ^[0-9]+\.[0-9]+(\.[0-9]+)?$ ]]; then
    echo "error: $PYTHON_SOURCE requested unsupported Python '$request'; CMM supports Python 3.10-3.12." >&2
    exit 1
  else
    local executable
    if [[ "$request" == */* || "$request" == *\\* ]]; then
      if [ ! -x "$request" ]; then
        echo "error: $PYTHON_SOURCE requested interpreter path '$request', but it is not an executable file." >&2
        exit 1
      fi
      executable="$request"
    else
      if ! executable="$(command -v "$request" 2>/dev/null)"; then
        echo "error: $PYTHON_SOURCE requested interpreter command '$request', but it was not found; use a path or a version such as 3.12." >&2
        exit 1
      fi
    fi
    if ! resolved="$("$UV_BIN" python find --no-project --no-python-downloads "$executable")" || [ -z "$resolved" ]; then
      managed_python_error "$PYTHON_SOURCE requested interpreter '$request', but uv could not resolve it."
    fi
  fi
  printf '%s\n' "$resolved"
}

validate_python() {
  local interpreter="$1"
  local context="$2"
  local version_output
  if ! version_output="$(supported_python_version "$interpreter" 2>&1)"; then
    echo "error: $context uses unsupported Python ${version_output:-unknown}; CMM supports Python 3.10-3.12." >&2
    exit 1
  fi
  printf '%s\n' "$version_output"
}

# --- select or install a supported Python -----------------------------------------------
RESOLVED_PYTHON=""
REQUESTED_VERSION=""
if [ -n "$PYTHON_REQUEST" ]; then
  RESOLVED_PYTHON="$(resolve_python_request "$PYTHON_REQUEST")"
  REQUESTED_VERSION="$(validate_python "$RESOLVED_PYTHON" "$PYTHON_SOURCE '$PYTHON_REQUEST'")"
fi

if [ -e "$VENV_DIR" ]; then
  if [ ! -d "$VENV_DIR" ]; then
    echo "error: virtual-environment path exists but is not a directory: $VENV_DIR" >&2
    exit 1
  fi
  if ! VPY="$(venv_python)"; then
    echo "error: $VENV_DIR exists but does not contain a Python virtual environment; choose another --venv path or remove it yourself." >&2
    exit 1
  fi
  if [ ! -f "$VENV_DIR/pyvenv.cfg" ]; then
    echo "error: $VENV_DIR has a Python executable but no pyvenv.cfg, so it is not a reusable virtual environment." >&2
    echo "Rename it yourself or choose a new path such as --venv .venv312; CMM will not delete or overwrite it." >&2
    exit 1
  fi
  if ! "$VPY" -c 'import sys; raise SystemExit(0 if sys.prefix != sys.base_prefix else 1)' >/dev/null 2>&1; then
    echo "error: $VENV_DIR/pyvenv.cfg exists, but '$VPY' reports sys.prefix == sys.base_prefix and is not running inside a virtual environment." >&2
    echo "Rename it yourself or choose a new path such as --venv .venv312; CMM will not delete or overwrite it." >&2
    exit 1
  fi
  VENV_VERSION="$(validate_python "$VPY" "existing virtual environment '$VENV_DIR'")"
  if [ -n "$RESOLVED_PYTHON" ] && [ "$(python_minor "$VENV_VERSION")" != "$(python_minor "$REQUESTED_VERSION")" ]; then
    echo "error: $PYTHON_SOURCE requested Python $(python_minor "$REQUESTED_VERSION"), but $VENV_DIR already uses Python $(python_minor "$VENV_VERSION"); choose another --venv path or remove the existing environment yourself." >&2
    exit 1
  fi
  echo "==> Reusing existing virtual environment at $VENV_DIR (Python $VENV_VERSION)"
else
  if [ -z "$RESOLVED_PYTHON" ]; then
    echo "==> Installing preferred Python $DEFAULT_PYTHON with uv"
    if ! "$UV_BIN" python install "$DEFAULT_PYTHON"; then
      managed_python_error "uv could not install preferred Python $DEFAULT_PYTHON."
    fi
    if ! RESOLVED_PYTHON="$("$UV_BIN" python find --no-project --managed-python "$DEFAULT_PYTHON")" || [ -z "$RESOLVED_PYTHON" ]; then
      managed_python_error "uv installed Python $DEFAULT_PYTHON but could not resolve its interpreter."
    fi
    REQUESTED_VERSION="$(validate_python "$RESOLVED_PYTHON" "uv-managed Python '$DEFAULT_PYTHON'")"
  fi
  echo "==> Creating virtual environment at $VENV_DIR with Python $REQUESTED_VERSION"
  if ! "$UV_BIN" venv --seed --python "$RESOLVED_PYTHON" "$VENV_DIR"; then
    managed_python_error "uv could not create '$VENV_DIR' with the resolved Python $REQUESTED_VERSION."
  fi
  if ! VPY="$(venv_python)"; then
    echo "error: uv did not create a usable Python environment at $VENV_DIR." >&2
    exit 1
  fi
fi

# --- install ---------------------------------------------------------------------------
if [ -n "$EXTRAS" ]; then
  echo "==> Installing CMM (editable) with extras: $EXTRAS"
  "$UV_BIN" pip install --python "$VPY" -e ".[${EXTRAS}]"
else
  echo "==> Installing CMM (editable, core only)"
  "$UV_BIN" pip install --python "$VPY" -e .
fi

if [ "$WITH_GUROBI" -eq 0 ]; then
  echo "==> Skipping gurobipy (GLPK supports LP/MILP; L2 MOMA/E-Flux2 need QP and published MTA/rMTA need MIQP)"
fi

# --- done ------------------------------------------------------------------------------
SOLVER="$("$VPY" -c 'from cmm.core import solver_status; s = solver_status(); print(s.name, "(" + ", ".join(s.capabilities) + ")")' 2>/dev/null || echo "unknown")"
cat <<EOF

==> CMM installed. Active solver: $SOLVER

   Launch the desktop application:
     $VPY -m cmm.app
   Or activate the environment first:
     source $VENV_DIR/bin/activate   # (Windows Git Bash: source $VENV_DIR/Scripts/activate)
     python -m cmm.app
EOF
