"""Strict, machine-readable catalog: only registered experiments may run."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CATALOG = Path(__file__).with_name("catalog.json")
EVIDENCE = {"hypothesis", "numerical", "interval_certified", "paper_proved", "lean_kernel_checked"}
METHODS = {"exact_rational", "interval", "ball_crosscheck", "floating", "conditional"}


def validate(record: dict, root: Path = ROOT) -> dict:
    mandatory = {"id", "path", "method", "evidence", "claim", "assumptions",
                 "falsifiers", "limitations", "timeout_s", "default_run"}
    absent = mandatory - set(record)
    if absent:
        raise ValueError(f"missing fields: {sorted(absent)}")
    if record["method"] not in METHODS or record["evidence"] not in EVIDENCE:
        raise ValueError("unrecognized method/evidence")
    if not isinstance(record["timeout_s"], int) or not 0 < record["timeout_s"] <= 300:
        raise ValueError("timeout must be an integer in 1..300 seconds")
    if not isinstance(record["default_run"], bool):
        raise ValueError("default_run must be boolean")
    for field in ("assumptions", "falsifiers", "limitations"):
        if not isinstance(record[field], list) or not record[field]:
            raise ValueError(f"{field} must be a nonempty list")
    candidate = root / record["path"]
    # Check containment before testing existence: never execute arbitrary paths.
    if not candidate.resolve().is_relative_to(root.resolve()):
        raise ValueError("experiment escapes repository root")
    if candidate.suffix != ".py" or not candidate.is_file():
        raise ValueError(f"missing Python experiment: {record['path']}")
    if "args" in record and (not isinstance(record["args"], list) or
                             not all(isinstance(v, str) for v in record["args"])):
        raise ValueError("args must be a list of strings")
    if "stdin_from" in record and not isinstance(record["stdin_from"], str):
        raise ValueError("stdin_from must identify a registered experiment")
    return record


def load(root: Path = ROOT, catalog: Path = CATALOG) -> dict[str, dict]:
    data = json.loads(catalog.read_text(encoding="utf-8"))
    if data.get("version") != 1 or not isinstance(data.get("experiments"), list):
        raise ValueError("unsupported catalog schema")
    results = {}
    for item in data["experiments"]:
        validate(item, root)
        if not isinstance(item["id"], str) or item["id"] in results or not item["id"].replace("-", "").replace("_", "").isalnum():
            raise ValueError("invalid or duplicate experiment id")
        results[item["id"]] = item
    for key, item in results.items():
        if "stdin_from" in item:
            dep = item["stdin_from"]
            if dep == key or dep not in results:
                raise ValueError("unknown or self-referencing stdin producer")
            if "stdin_from" in results[dep]:
                raise ValueError("nested stdin producers are not supported")
    return results
