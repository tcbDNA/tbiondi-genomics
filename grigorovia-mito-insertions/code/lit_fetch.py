"""Find papers on Europe PMC and save their metadata, abstract and (if open access) full-text XML in lit/.

    uv run python code/lit_fetch.py search "<query>" [n]
    uv run python code/lit_fetch.py fetch <PMCID|PMID|DOI> <short_name>

Used for the prior-art step only; every file it saves is a source for PRIOR_ART.md.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

LIT = Path(__file__).resolve().parents[1] / "lit"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest"
UA = {"User-Agent": "grigorovia-mito-insertions/0.1 (literature search)"}


def get(url: str) -> bytes:
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read()
        except Exception:  # noqa: BLE001
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError


def search(query: str, n: int = 5) -> list[dict]:
    q = urllib.parse.quote(query)
    d = json.loads(get(f"{API}/search?query={q}&format=json&resultType=core&pageSize={n}"))
    return d["resultList"]["result"]


def show(rs: list[dict]) -> None:
    for r in rs:
        print(f"- {r.get('title')}\n  {r.get('authorString', '')[:90]} | {r.get('journalTitle', r.get('bookOrReportDetails', {}).get('publisher', ''))} "
              f"{r.get('pubYear')} | doi {r.get('doi')} | pmid {r.get('pmid')} | {r.get('pmcid') or 'no PMC'} | OA {r.get('isOpenAccess')}")


def fetch(ident: str, name: str) -> None:
    if ident.upper().startswith("PMC"):
        q = f"PMCID:{ident}"
    elif ident.isdigit():
        q = f"EXT_ID:{ident} AND SRC:MED"
    else:
        q = f'DOI:"{ident}"'
    rs = search(q, 1)
    if not rs:
        sys.exit(f"not found: {ident}")
    r = rs[0]
    LIT.mkdir(exist_ok=True)
    (LIT / f"{name}.json").write_text(json.dumps(r, indent=1))
    out = f"{name}: {r.get('title')} ({r.get('pubYear')}) doi {r.get('doi')}"
    if r.get("pmcid") and r.get("inEPMC") == "Y":
        try:
            xml = get(f"{API}/{r['pmcid']}/fullTextXML")
            (LIT / f"{name}.xml").write_bytes(xml)
            out += f" | full text {len(xml) // 1024} kB"
        except Exception as e:  # noqa: BLE001
            out += f" | full text failed ({e.__class__.__name__})"
    else:
        out += " | abstract only"
    print(out)


if __name__ == "__main__":
    if sys.argv[1] == "search":
        show(search(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 5))
    elif sys.argv[1] == "fetch":
        fetch(sys.argv[2], sys.argv[3])
