from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
INSTALL_SH = ROOT / "install.sh"
INSTALL_PS1 = ROOT / "install.ps1"
BASH = shutil.which("bash")
SHELL_INSTALLER_TEST = pytest.mark.skipif(
    os.name == "nt" or BASH is None,
    reason="the shell-installer behavior tests require a native Unix Bash environment",
)


def _write_executable(path: Path, contents: str) -> None:
    path.write_text(contents, encoding="utf-8")
    path.chmod(0o755)


def _fake_python(
    path: Path, version: str, *, supported: bool, in_venv: bool = False
) -> None:
    exit_code = 0 if supported else 1
    venv_exit_code = 0 if in_venv else 1
    _write_executable(
        path,
        f"""#!/usr/bin/env bash
if [ "${{1:-}}" = "-c" ] && [[ "${{2:-}}" == *sys.prefix* ]]; then
  exit {venv_exit_code}
fi
if [ "${{1:-}}" = "-c" ] && [[ "${{2:-}}" == *sys.version_info* ]]; then
  printf '%s\\n' '{version}'
  exit {exit_code}
fi
if [ "${{1:-}}" = "-c" ] && [[ "${{2:-}}" == *solver_status* ]]; then
  printf '%s\\n' 'stub (LP)'
  exit 0
fi
exit 0
""",
    )


def _installer_environment(tmp_path: Path) -> tuple[dict[str, str], Path, Path, Path]:
    stub_bin = tmp_path / "stub bin"
    stub_bin.mkdir()
    managed_312 = tmp_path / "managed python 3.12"
    managed_311 = tmp_path / "managed python 3.11"
    _fake_python(managed_312, "3.12.9", supported=True)
    _fake_python(managed_311, "3.11.8", supported=True)

    uv_log = tmp_path / "uv.log"
    system_python_log = tmp_path / "system-python.log"
    uv_stub = stub_bin / "uv"
    _write_executable(
        uv_stub,
        """#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$*" >> "$CMM_TEST_UV_LOG"
case "${1:-}" in
  --version)
    printf 'uv %s\\n' "${CMM_TEST_UV_VERSION:-0.12.5}"
    ;;
  python)
    case "${2:-}" in
      install)
        if [ "${CMM_TEST_FAIL_PYTHON_INSTALL:-0}" = "1" ]; then exit 71; fi
        ;;
      find)
        if [ "${CMM_TEST_FAIL_PYTHON_FIND:-0}" = "1" ]; then exit 72; fi
        request="${!#}"
        case "$request" in
          3.12|3.12.*) printf '%s\\n' "$CMM_TEST_PYTHON_312" ;;
          3.11|3.11.*) printf '%s\\n' "$CMM_TEST_PYTHON_311" ;;
          *) printf '%s\\n' "$request" ;;
        esac
        ;;
      *) exit 2 ;;
    esac
    ;;
  venv)
    interpreter=""
    destination="${!#}"
    for ((index = 1; index <= $#; index++)); do
      if [ "${!index}" = "--python" ]; then
        next=$((index + 1))
        interpreter="${!next}"
      fi
    done
    mkdir -p "$destination/bin"
    ln -s "$interpreter" "$destination/bin/python"
    ;;
  pip)
    ;;
  *) exit 2 ;;
esac
""",
    )
    _write_executable(
        stub_bin / "python3",
        """#!/usr/bin/env bash
printf '%s\\n' invoked >> "$CMM_TEST_SYSTEM_PYTHON_LOG"
printf '%s\\n' '3.9.6'
exit 1
""",
    )

    environment = os.environ.copy()
    environment.update(
        {
            "PATH": f"{stub_bin}{os.pathsep}{environment['PATH']}",
            "HOME": str(tmp_path / "home"),
            "CMM_TEST_UV_LOG": str(uv_log),
            "CMM_TEST_SYSTEM_PYTHON_LOG": str(system_python_log),
            "CMM_TEST_PYTHON_312": str(managed_312),
            "CMM_TEST_PYTHON_311": str(managed_311),
        }
    )
    environment.pop("CMM_PYTHON", None)
    return environment, uv_log, system_python_log, stub_bin


def _run_install_sh(
    tmp_path: Path, *arguments: str, environment: dict[str, str] | None = None
) -> tuple[subprocess.CompletedProcess[str], str, Path]:
    if environment is None:
        environment, uv_log, _, _ = _installer_environment(tmp_path)
    else:
        uv_log = Path(environment["CMM_TEST_UV_LOG"])
    venv = tmp_path / "CMM environment with spaces"
    assert BASH is not None
    completed = subprocess.run(
        [BASH, str(INSTALL_SH), "--venv", str(venv), *arguments],
        cwd=tmp_path,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    log = uv_log.read_text(encoding="utf-8") if uv_log.exists() else ""
    return completed, log, venv


@SHELL_INSTALLER_TEST
def test_shell_installer_bootstraps_managed_312_without_system_python(
    tmp_path: Path,
) -> None:
    environment, _, system_python_log, _ = _installer_environment(tmp_path)

    completed, log, venv = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode == 0, completed.stderr
    assert "python install 3.12" in log
    assert "python find --no-project --managed-python 3.12" in log
    assert f"venv --seed --python {environment['CMM_TEST_PYTHON_312']} {venv}" in log
    assert (
        f"pip install --python {venv}/bin/python -e .[desktop,design,solver-gurobi]"
    ) in log
    assert "Python 3.12.9" in completed.stdout
    assert "stub (LP)" in completed.stdout
    assert not system_python_log.exists(), "the macOS-style system Python was invoked"


@SHELL_INSTALLER_TEST
def test_shell_installer_bootstraps_truly_absent_uv_and_finds_home_binary(
    tmp_path: Path,
) -> None:
    environment, _, _, stub_bin = _installer_environment(tmp_path)
    uv_payload = tmp_path / "uv bootstrap payload"
    (stub_bin / "uv").replace(uv_payload)
    fetch_log = tmp_path / "fetch.log"
    environment.update(
        {
            "PATH": f"{stub_bin}{os.pathsep}/usr/bin{os.pathsep}/bin",
            "CMM_TEST_BOOTSTRAP_UV": str(uv_payload),
            "CMM_TEST_FETCH_LOG": str(fetch_log),
        }
    )
    _write_executable(
        stub_bin / "curl",
        """#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$*" > "$CMM_TEST_FETCH_LOG"
cat <<'CMM_BOOTSTRAP'
mkdir -p "$HOME/.local/bin"
cp "$CMM_TEST_BOOTSTRAP_UV" "$HOME/.local/bin/uv"
chmod 755 "$HOME/.local/bin/uv"
CMM_BOOTSTRAP
""",
    )
    assert shutil.which("uv", path=environment["PATH"]) is None

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    home_uv = Path(environment["HOME"]) / ".local" / "bin" / "uv"
    assert completed.returncode == 0, completed.stderr
    assert home_uv.is_file()
    assert "uv 0.12.5" in completed.stdout
    assert f"({home_uv})" in completed.stdout
    assert "https://astral.sh/uv/0.12.5/install.sh" in fetch_log.read_text(
        encoding="utf-8"
    )
    assert "python install 3.12" in log


@SHELL_INSTALLER_TEST
def test_shell_installer_accepts_documented_minimum_uv(tmp_path: Path) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    environment["CMM_TEST_UV_VERSION"] = "0.8.0"

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode == 0, completed.stderr
    assert "Using uv 0.8.0" in completed.stdout
    assert "python install 3.12" in log
    assert "pip install" in log


@SHELL_INSTALLER_TEST
def test_shell_installer_rejects_too_old_installed_uv_with_remediation(
    tmp_path: Path,
) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    environment["CMM_TEST_UV_VERSION"] = "0.7.9"

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode != 0
    assert "requires uv >= 0.8.0" in completed.stderr
    assert "found 0.7.9" in completed.stderr
    assert "https://astral.sh/uv/0.12.5/install.sh" in completed.stderr
    assert '. "$HOME/.local/bin/env"' in completed.stderr
    assert "python install" not in log


@SHELL_INSTALLER_TEST
def test_shell_installer_installs_supported_numeric_override(tmp_path: Path) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    environment["CMM_PYTHON"] = "3.11"

    completed, log, _ = _run_install_sh(
        tmp_path,
        "--core-only",
        "--no-gurobi",
        environment=environment,
    )

    assert completed.returncode == 0, completed.stderr
    assert log.index("python install 3.11") < log.index(
        "python find --no-project --managed-python 3.11"
    )
    assert "pip install" in log and " -e .\n" in log
    assert ".[desktop" not in log


@SHELL_INSTALLER_TEST
def test_shell_installer_rejects_unsupported_interpreter_path(tmp_path: Path) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    unsupported = tmp_path / "existing python 3.9"
    _fake_python(unsupported, "3.9.6", supported=False)
    environment["CMM_PYTHON"] = str(unsupported)

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode != 0
    assert "CMM_PYTHON" in completed.stderr
    assert "unsupported Python 3.9.6" in completed.stderr
    assert "CMM supports Python 3.10-3.12" in completed.stderr
    assert "--no-python-downloads" in log
    assert "python install 3.9" not in log
    assert "pip install" not in log


@SHELL_INSTALLER_TEST
def test_shell_installer_rejects_missing_interpreter_path_without_download(
    tmp_path: Path,
) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    missing = tmp_path / "missing Python" / "python"
    environment["CMM_PYTHON"] = str(missing)

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode != 0
    assert f"interpreter path '{missing}'" in completed.stderr
    assert "not an executable file" in completed.stderr
    assert "python install" not in log
    assert "python find" not in log
    assert "pip install" not in log


@SHELL_INSTALLER_TEST
def test_shell_installer_rejects_existing_venv_with_old_python(
    tmp_path: Path,
) -> None:
    environment, _, system_python_log, _ = _installer_environment(tmp_path)
    venv = tmp_path / "CMM environment with spaces"
    (venv / "bin").mkdir(parents=True)
    (venv / "pyvenv.cfg").write_text("version = 3.9.6\n", encoding="utf-8")
    _fake_python(venv / "bin" / "python", "3.9.6", supported=False, in_venv=True)

    completed, log, returned_venv = _run_install_sh(tmp_path, environment=environment)

    assert returned_venv == venv
    assert completed.returncode != 0
    assert f"existing virtual environment '{venv}'" in completed.stderr
    assert "unsupported Python 3.9.6" in completed.stderr
    assert "python install" not in log
    assert "python find" not in log
    assert "pip install" not in log
    assert not system_python_log.exists(), "the system Python was invoked"


@SHELL_INSTALLER_TEST
def test_shell_installer_reuses_only_a_verified_supported_venv(tmp_path: Path) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    venv = tmp_path / "CMM environment with spaces"
    (venv / "bin").mkdir(parents=True)
    (venv / "pyvenv.cfg").write_text("version = 3.12.9\n", encoding="utf-8")
    _fake_python(venv / "bin" / "python", "3.12.9", supported=True, in_venv=True)

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode == 0, completed.stderr
    assert "Reusing existing virtual environment" in completed.stdout
    assert "python install" not in log
    assert "python find" not in log
    assert f"pip install --python {venv}/bin/python" in log


@SHELL_INSTALLER_TEST
def test_shell_installer_requires_pyvenv_configuration_for_existing_directory(
    tmp_path: Path,
) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    venv = tmp_path / "CMM environment with spaces"
    (venv / "bin").mkdir(parents=True)
    _fake_python(venv / "bin" / "python", "3.12.9", supported=True, in_venv=True)

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode != 0
    assert "no pyvenv.cfg" in completed.stderr
    assert "--venv .venv312" in completed.stderr
    assert "pip install" not in log


@SHELL_INSTALLER_TEST
def test_shell_installer_rejects_plain_python_behind_pyvenv_marker(
    tmp_path: Path,
) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    venv = tmp_path / "CMM environment with spaces"
    (venv / "bin").mkdir(parents=True)
    (venv / "pyvenv.cfg").write_text("version = 3.12.9\n", encoding="utf-8")
    _fake_python(venv / "bin" / "python", "3.12.9", supported=True)

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode != 0
    assert "sys.prefix == sys.base_prefix" in completed.stderr
    assert "--venv .venv312" in completed.stderr
    assert "pip install" not in log


@SHELL_INSTALLER_TEST
@pytest.mark.parametrize(
    ("failure_variable", "expected_error"),
    [
        ("CMM_TEST_FAIL_PYTHON_INSTALL", "could not install preferred Python 3.12"),
        (
            "CMM_TEST_FAIL_PYTHON_FIND",
            "installed Python 3.12 but could not resolve its interpreter",
        ),
    ],
)
def test_shell_installer_explains_managed_python_failures(
    tmp_path: Path, failure_variable: str, expected_error: str
) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    environment[failure_variable] = "1"

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode != 0
    assert expected_error in completed.stderr
    assert "CMM supports Python 3.10-3.12" in completed.stderr
    assert "CMM_PYTHON=/path/to/python3.12" in completed.stderr
    assert "pip install" not in log


@SHELL_INSTALLER_TEST
def test_shell_installer_rejects_unsupported_numeric_override_before_download(
    tmp_path: Path,
) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    environment["CMM_PYTHON"] = "3.13"

    completed, log, _ = _run_install_sh(tmp_path, environment=environment)

    assert completed.returncode != 0
    assert "requested unsupported Python '3.13'" in completed.stderr
    assert "python install 3.13" not in log
    assert "python find" not in log


@SHELL_INSTALLER_TEST
def test_shell_python_flag_precedes_environment_override(tmp_path: Path) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)
    environment["CMM_PYTHON"] = "3.11"

    completed, log, _ = _run_install_sh(
        tmp_path,
        "--python",
        "3.12",
        "--no-gurobi",
        environment=environment,
    )

    assert completed.returncode == 0, completed.stderr
    assert "python install 3.12" in log
    assert "python install 3.11" not in log
    assert ".[desktop,design]" in log


@SHELL_INSTALLER_TEST
@pytest.mark.parametrize(
    ("arguments", "expected_extra"),
    [
        (("--dev", "--no-gurobi"), ".[desktop,design,dev]"),
        (("--core-only",), ".[solver-gurobi]"),
    ],
)
def test_shell_installer_preserves_extra_flags(
    tmp_path: Path, arguments: tuple[str, ...], expected_extra: str
) -> None:
    environment, _, _, _ = _installer_environment(tmp_path)

    completed, log, _ = _run_install_sh(tmp_path, *arguments, environment=environment)

    assert completed.returncode == 0, completed.stderr
    assert expected_extra in log


def test_installer_bootstrap_urls_are_version_pinned() -> None:
    shell = INSTALL_SH.read_text(encoding="utf-8")
    powershell = INSTALL_PS1.read_text(encoding="utf-8")

    assert 'UV_BOOTSTRAP_VERSION="0.12.5"' in shell
    assert '"https://astral.sh/uv/$UV_BOOTSTRAP_VERSION/install.sh"' in shell
    assert '$uvBootstrapVersion = "0.12.5"' in powershell
    assert '"https://astral.sh/uv/$uvBootstrapVersion/install.ps1"' in powershell
    assert "https://astral.sh/uv/install.sh" not in shell
    assert "https://astral.sh/uv/install.ps1" not in powershell
    assert 'UV_MINIMUM_VERSION="0.8.0"' in shell
    assert 'version_at_least "$installed_version" "$UV_MINIMUM_VERSION"' in shell
    assert '$uvMinimumVersion = "0.8.0"' in powershell
    assert "[version]$installedVersion -lt [version]$uvMinimumVersion" in powershell
    assert "pyvenv.cfg" in shell
    assert "pyvenv.cfg" in powershell
    assert "sys.prefix != sys.base_prefix" in shell
    assert "sys.prefix != sys.base_prefix" in powershell


@pytest.mark.skipif(
    shutil.which("pwsh") is None and shutil.which("powershell") is None,
    reason="PowerShell is not installed on this platform",
)
def test_powershell_installer_parses() -> None:
    executable = shutil.which("pwsh") or shutil.which("powershell")
    assert executable is not None
    command = (
        "$errors = $null; "
        f"[void][System.Management.Automation.Language.Parser]::ParseFile('{INSTALL_PS1}', "
        "[ref]$null, [ref]$errors); "
        "if ($errors.Count -gt 0) { $errors | ForEach-Object { Write-Error $_ }; exit 1 }"
    )
    subprocess.run(
        [executable, "-NoProfile", "-NonInteractive", "-Command", command],
        check=True,
        text=True,
    )
