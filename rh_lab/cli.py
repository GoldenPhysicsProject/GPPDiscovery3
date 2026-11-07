"""Command line: python -m rh_lab.cli list|check|run <id> [--replay]|smoke."""

import argparse
import json
import sys

from .catalog import load
from .runner import run


def main() -> None:
    p = argparse.ArgumentParser(description="Codex RH discovery laboratory")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list")
    sub.add_parser("check")
    r = sub.add_parser("run")
    r.add_argument("id")
    r.add_argument("--replay", action="store_true")
    sub.add_parser("smoke")
    args = p.parse_args()
    records = load()
    if args.command == "list":
        for key, value in records.items():
            print(f"{key:27} {value['evidence']:20} {value['method']}")
    elif args.command == "check":
        print(f"Validated {len(records)} catalog entries; no unregistered paths allowed")
    elif args.command == "run":
        if args.id not in records:
            p.error("unknown experiment identifier")
        result = run(args.id, replay=args.replay)
        print(json.dumps({k: v for k, v in result.items() if k != "replay"}, indent=2))
        if not result["passed"]:
            sys.exit(1)
    elif args.command == "smoke":
        for key, record in records.items():
            if record["default_run"]:
                result = run(key, replay=True)
                print(f"{key}: {'PASS' if result['passed'] else 'FAIL'}, record={result['artifact_path']}")
                if not result["passed"]:
                    sys.exit(1)


if __name__ == "__main__":
    main()
