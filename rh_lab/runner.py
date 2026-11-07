"""Recorded, bounded, deterministic subprocess runner for trusted repo scripts.

This is an audit/replay mechanism, NOT a security sandbox. Only locally
trusted code from the catalog should be invoked. It does not disable network.
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from .catalog import ROOT, CATALOG, load


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run(name: str, *, replay: bool = False, root: Path = ROOT,
        catalog: Path = CATALOG, output_dir: Path | None = None) -> dict:
    records = load(root, catalog)
    record = records[name]
    source = root / record["path"]
    argv = [sys.executable, str(source), *record.get("args", [])]
    env = dict(os.environ, PYTHONHASHSEED="0", TZ="UTC")
    producer = records.get(record.get("stdin_from"))

    def execute() -> dict:
        started = time.monotonic()
        try:
            input_data = None
            producer_hash = None
            if producer is not None:
                gen_source = root / producer["path"]
                gen = subprocess.run([sys.executable, str(gen_source), *producer.get("args", [])],
                                     cwd=root, capture_output=True, timeout=producer["timeout_s"],
                                     check=False, env=env)
                if gen.returncode != 0:
                    return {"returncode": gen.returncode, "timeout": False,
                            "error": "stdin producer failed", "stderr": gen.stderr.decode("utf-8", "replace")[:12000],
                            "elapsed_seconds": round(time.monotonic() - started, 6)}
                input_data = gen.stdout
                producer_hash = sha(input_data)
            p = subprocess.run(argv, cwd=root, capture_output=True, timeout=record["timeout_s"],
                               check=False, env=env, input=input_data)
            return {"returncode": p.returncode, "timeout": False,
                    "input_sha256": producer_hash,
                    "stdout_sha256": sha(p.stdout), "stderr_sha256": sha(p.stderr),
                    "stdout": p.stdout.decode("utf-8", "replace")[:12000],
                    "stderr": p.stderr.decode("utf-8", "replace")[:12000],
                    "elapsed_seconds": round(time.monotonic() - started, 6)}
        except subprocess.TimeoutExpired as exc:
            return {"returncode": None, "timeout": True,
                    "stdout": (exc.stdout or b"").decode("utf-8", "replace")[:12000],
                    "stderr": (exc.stderr or b"").decode("utf-8", "replace")[:12000],
                    "elapsed_seconds": round(time.monotonic() - started, 6)}

    first = execute()
    second = execute() if replay and first["returncode"] == 0 else None
    deterministic = (second is not None and second["returncode"] == 0 and
                     all(first.get(k) == second.get(k)
                         for k in ("returncode", "input_sha256", "stdout_sha256", "stderr_sha256"))) if replay else None
    result = {"schema_version": 1, "run_id": str(uuid.uuid4()),
              "timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "experiment": name, "claim": record["claim"],
              "evidence_declared": record["evidence"],
              "code_sha256": sha(source.read_bytes()),
              "stdin_source_sha256": sha((root / producer["path"]).read_bytes()) if producer else None,
              "catalog_sha256": sha(catalog.read_bytes()),
              "argv": [sys.executable, record["path"], *record.get("args", [])],
              "python": sys.version.split()[0], "first": first,
              "replay": second, "deterministic": deterministic,
              "passed": first["returncode"] == 0 and (not replay or deterministic is True),
              "warning": "Execution never upgrades the declared proof status; claims need independent source/error audit."}
    dest = output_dir if output_dir is not None else root / ".rh_runs"
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / f"{result['run_id']}.json"
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=dest, delete=False) as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")
        temp = Path(f.name)
    temp.replace(target)
    result["artifact_path"] = str(target)
    return result
