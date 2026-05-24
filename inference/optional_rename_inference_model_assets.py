"""Rename timestamped assets under ``inference/models`` and rewrite JSON cards.

This utility is intentionally scoped to the local ``inference/models`` folder only.
It removes the ``_YYYYMMDD_HHMMSS`` timestamp segment from model asset filenames,
then updates the embedded filename fields in each JSON model card so the local
references continue to match the renamed files.

Examples
--------
    SVR_air_hsbd_20260516_145335_ad.npz -> SVR_air_hsbd_ad.npz
    SVR_sediment_hsbd_20260516_150639.joblib -> SVR_sediment_hsbd.joblib
    training_split_sediment_hsbd_20260516_150639.csv -> training_split_sediment_hsbd.csv

Run with ``--dry-run`` to preview changes. Without ``--dry-run``, the script
applies the renames and JSON rewrites.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

TIMESTAMP_PATTERN = re.compile(r"_(\d{8}_\d{6})(?=(_|\.))")
JSON_PATH_KEYS = ("model_file", "ad_artefact_file", "training_split_file", "t_half_meta_file")


def strip_timestamp(name: str) -> str:
    """Remove the timestamp segment from a filename.

    The function only removes the first timestamp-like segment of the form
    ``_YYYYMMDD_HHMMSS`` when it is followed by either another underscore or the
    file extension.
    """

    return TIMESTAMP_PATTERN.sub("", name, count=1)


def build_rename_plan(models_dir: Path) -> tuple[list[tuple[Path, Path]], dict[str, str]]:
    """Return a rename plan and basename mapping for files in ``models_dir``."""

    if not models_dir.is_dir():
        raise FileNotFoundError(f"Models directory does not exist: {models_dir}")

    rename_plan: list[tuple[Path, Path]] = []
    rename_map: dict[str, str] = {}
    sources = set()

    for path in sorted(models_dir.iterdir()):
        if not path.is_file():
            continue

        new_name = strip_timestamp(path.name)
        if new_name == path.name:
            continue

        if new_name in rename_map and rename_map[new_name] != path.name:
            raise ValueError(
                f"Multiple files would map to the same target name {new_name!r}: {rename_map[new_name]!r} and {path.name!r}"
            )

        rename_map[path.name] = new_name
        rename_plan.append((path, path.with_name(new_name)))
        sources.add(path)

    for source, target in rename_plan:
        if target.exists() and target not in sources:
            raise FileExistsError(f"Target already exists and is not being renamed: {target}")

    return rename_plan, rename_map


def rewrite_path_string(value: str, rename_map: dict[str, str]) -> tuple[str, bool]:
    """Rewrite a path-like string if its basename is in ``rename_map``."""

    basename = Path(value).name
    if basename not in rename_map:
        return value, False

    new_value = value[: len(value) - len(basename)] + rename_map[basename]
    return new_value, new_value != value


def rewrite_json_node(node: Any, rename_map: dict[str, str]) -> tuple[Any, bool]:
    """Recursively rewrite timestamped path strings inside a JSON structure."""

    if isinstance(node, dict):
        changed = False
        rewritten = {}
        for key, value in node.items():
            if key in JSON_PATH_KEYS and isinstance(value, str):
                new_value, value_changed = rewrite_path_string(value, rename_map)
                rewritten[key] = new_value
                changed |= value_changed
                continue

            new_value, value_changed = rewrite_json_node(value, rename_map)
            rewritten[key] = new_value
            changed |= value_changed

        return rewritten, changed

    if isinstance(node, list):
        changed = False
        rewritten_list = []
        for item in node:
            new_item, item_changed = rewrite_json_node(item, rename_map)
            rewritten_list.append(new_item)
            changed |= item_changed
        return rewritten_list, changed

    if isinstance(node, str):
        return rewrite_path_string(node, rename_map)

    return node, False


def rename_files(rename_plan: list[tuple[Path, Path]], dry_run: bool) -> None:
    """Rename files in two phases so target names do not collide mid-run."""

    temp_paths: list[tuple[Path, Path]] = []

    for index, (source, target) in enumerate(rename_plan):
        temp_path = source.with_name(f".{source.name}.rename_tmp_{index}")
        temp_paths.append((temp_path, target))
        if dry_run:
            print(f"DRY-RUN rename: {source.name} -> {target.name}")
            continue
        source.rename(temp_path)

    if dry_run:
        return

    for temp_path, target in temp_paths:
        temp_path.rename(target)


def rewrite_json_cards(models_dir: Path, rename_map: dict[str, str], dry_run: bool) -> None:
    """Rewrite the JSON cards in ``models_dir`` using the new filenames."""

    for json_path in sorted(models_dir.glob("*.json")):
        with json_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)

        rewritten_payload, changed = rewrite_json_node(payload, rename_map)
        if not changed:
            if dry_run:
                print(f"DRY-RUN JSON unchanged: {json_path.name}")
            continue

        if dry_run:
            print(f"DRY-RUN JSON rewrite: {json_path.name}")
            continue

        with json_path.open("w", encoding="utf-8") as handle:
            json.dump(rewritten_payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rename timestamped assets and rewrite JSON cards in inference/models.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Preview the rename plan without modifying files.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    models_dir = Path(__file__).resolve().parent / "models"

    rename_plan, rename_map = build_rename_plan(models_dir)
    if not rename_plan:
        print(f"No timestamped assets found in {models_dir}")
        return 0

    rename_files(rename_plan, dry_run=args.dry_run)
    rewrite_json_cards(models_dir, rename_map, dry_run=args.dry_run)

    if args.dry_run:
        print("Dry-run complete. Re-run without --dry-run to apply changes.")
    else:
        print(f"Renamed {len(rename_plan)} files and updated JSON cards in {models_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
