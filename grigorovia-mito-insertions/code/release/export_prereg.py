"""Write the release copy of PREREG.md (protocol/frozen_rules.md in the release): a provenance header, PREREG.md verbatim, and commit tables from local git.

    uv run python code/release/export_prereg.py OUT_PATH
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HEADER = """# Pre-registered rules, predictions and changes

**Provenance.**
- This file is a copy of `PREREG.md` from the project's working repository. It is reproduced verbatim except for this
  header and the commit tables at the end.
- The working repository is a **local git repository**. The commit hashes and times below come from it and were **not
  externally timestamped**: there is no public remote and no third-party timestamping service. They show the order in
  which entries were committed, as recorded by local git; they are not independent proof of timing.
- **History rewrite.** On 2026-11-22 the local history after commit d6e06a9 was rewritten to remove a 611 MB
  intermediate file committed by mistake. Messages, authors and dates were preserved, and no hash cited in this file
  changed. The old-to-new hash table is in the log entry "Local history rewritten in this session".
- **Paths.** The log refers to files by their paths in the working repository. This release renamed several of
  them; `protocol/corrections.md` has the old-to-new table. The log itself is not edited.
- Times are US Eastern (EST), 2026-11-21 to 2026-11-22. Where an entry's typed time differed from its commit time, the
  log says so (entry "correction of timestamps").

---

"""


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def table(paths=()):
    out = git("log", "--reverse", "--date=format:%Y-%m-%d %H:%M", "--format=%h\t%ad\t%s", "--", *paths) if paths else \
        git("log", "--reverse", "--date=format:%Y-%m-%d %H:%M", "--format=%h\t%ad\t%s")
    rows = ["| commit | time | message |", "|---|---|---|"]
    for line in out.splitlines():
        h, t, m = line.split("\t", 2)
        rows.append(f"| {h} | {t} | {m.replace('|', '/')} |")
    return "\n".join(rows)


def main(out):
    body = (ROOT / "PREREG.md").read_text()
    bad = re.findall(r"/Users/[^/\s]+|/home/[^/\s]+|/private/tmp|/opt/homebrew", body)
    if bad:
        sys.exit(f"local paths in PREREG.md: {set(bad)}")
    text = (HEADER + body.rstrip() + "\n\n---\n\n## Commits touching PREREG.md or DESIGN.md (local git; hash, commit time EST, message)\n\n"
            + table(("PREREG.md", "DESIGN.md")) + "\n\n## All project commits (local git)\n\n" + table() + "\n")
    Path(out).write_text(text)
    print(f"wrote {out}: {len(text.splitlines())} lines")


if __name__ == "__main__":
    main(sys.argv[1])
