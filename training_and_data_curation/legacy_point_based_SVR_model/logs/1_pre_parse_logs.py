"""concatenates all log files into one, removing timestamps and lines with .png  and adding separators between logs for readability.
Further parsing possible with script 2_log_review.py which summarizes into .md file.

output: logs/combined_logs.txt
"""

#!/usr/bin/env python3
import re
from pathlib import Path

LOG_DIR = Path(__file__).parent
output_file = LOG_DIR / "combined_logs.txt"
SEPARATOR = "#" * 60


def find_log_files(root: Path) -> list[Path]:
    return list(root.rglob("SVR_point_model.log"))


def process_log_content(content: str) -> str:
    lines = content.splitlines()
    processed_lines = []

    timestamp_pattern = re.compile(r"^\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\]\s*")

    for line in lines:
        if ".png" in line:
            continue
        processed_line = timestamp_pattern.sub("", line)
        processed_lines.append(processed_line)

    return "\n".join(processed_lines)


def main():
    log_files = find_log_files(LOG_DIR)
    # ignore files in _log_summaries
    log_files = [f for f in log_files if "_log_by_occasion_manually_moved" not in str(f)]
    log_files.sort(key=lambda p: p.parent.name)

    all_logs = []

    for log_path in log_files:
        rel_path = log_path.parent.relative_to(LOG_DIR)
        content = log_path.read_text()
        processed = process_log_content(content)

        all_logs.append(SEPARATOR)
        all_logs.append(f"# {rel_path}")
        all_logs.append(SEPARATOR)
        all_logs.append("")
        all_logs.append(processed)
        all_logs.append("")

    output_file.write_text("\n".join(all_logs))
    print(f"Combined {len(log_files)} log files into {output_file}")


if __name__ == "__main__":
    main()
