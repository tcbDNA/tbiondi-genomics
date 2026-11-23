"""Prior-art audit: follow seed papers forward to the papers citing them (Europe PMC + OpenAlex).

Writes lit/forward/<seed>_citing.tsv (every citing work) and lit/forward/screened.tsv (those whose
title or abstract mentions a frameshift-related term). Not part of the scan.

    uv run python code/audit/forward_citations.py
"""

from __future__ import annotations

import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "lit/forward"
UA = {"User-Agent": "grigorovia-mito-insertions/0.1 (literature audit)"}

SEEDS = {  # name: (PMID, DOI)
    "al_arab_2017": ("27693569", "10.1016/j.ympev.2016.09.024"),
    "gornik_2022": ("36108082", "10.1093/molbev/msac191"),
    "szaboova_2018": ("29528462", "10.1093/femsyr/foy028"),
    "rosengarten_2008": ("18215303", "10.1186/1471-2164-9-33"),
    "mitozoa_2012": ("22123747", "10.1093/nar/gkr1144"),
    # added in the audit (not in the user's list): key papers found during the audit
    "lang_2014": ("24711422", "10.1073/pnas.1322190111"),
    "haen_2014": ("24177232", "10.1016/j.gene.2013.10.041"),
    "masuda_2010": ("20507907", "10.1093/nar/gkq449"),
    "russell_beckenbach_2008": ("19030769", "10.1007/s00239-008-9179-0"),
}
TERMS = re.compile(r"frame-?\s?shift|ribosomal slippage|recoding|quadruplet|non-?triplet|"
                   r"\+1 (?:shift|insertion)|insertion editing|programmed (?:translational )?frameshift", re.I)


def get_json(url: str) -> dict:
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return json.loads(r.read())
        except Exception:  # noqa: BLE001
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError


def epmc_citing(pmid: str) -> list[dict]:
    out, page = [], 1
    while True:
        d = get_json(f"https://www.ebi.ac.uk/europepmc/webservices/rest/MED/{pmid}/citations"
                     f"?format=json&pageSize=1000&page={page}")
        cits = d.get("citationList", {}).get("citation", [])
        out += [{"src": "EPMC", "id": f"{c.get('source')}:{c.get('id')}", "doi": "", "year": c.get("pubYear"),
                 "type": c.get("source"), "title": c.get("title", ""), "abstract": ""} for c in cits]
        if len(cits) < 1000:
            return out
        page += 1


def openalex_citing(doi: str) -> list[dict]:
    w = get_json(f"https://api.openalex.org/works/doi:{urllib.parse.quote(doi)}")
    wid = w["id"].rsplit("/", 1)[-1]
    out, cursor = [], "*"
    while cursor:
        d = get_json(f"https://api.openalex.org/works?filter=cites:{wid}&per-page=200&cursor={cursor}"
                     f"&select=id,doi,title,publication_year,type,abstract_inverted_index")
        for x in d["results"]:
            inv = x.get("abstract_inverted_index") or {}
            words = sorted(((p, wd) for wd, ps in inv.items() for p in ps))
            out.append({"src": "OpenAlex", "id": x["id"].rsplit("/", 1)[-1], "doi": (x.get("doi") or "").replace("https://doi.org/", ""),
                        "year": x.get("publication_year"), "type": x.get("type"), "title": x.get("title") or "",
                        "abstract": " ".join(wd for _, wd in words)})
        cursor = d["meta"].get("next_cursor")
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    screened = []
    print(f"{'seed':18s} {'EPMC':>5s} {'OpenAlex':>8s} {'union':>6s} {'preprints':>9s} {'2025-26':>7s} {'term hits':>9s}")
    for name, (pmid, doi) in SEEDS.items():
        rows = epmc_citing(pmid) + openalex_citing(doi)
        seen, uniq = set(), []
        for r in rows:
            key = r["doi"].lower() or re.sub(r"\W+", "", r["title"].lower())[:80]
            if key in seen:
                continue
            seen.add(key)
            uniq.append(r)
        with open(OUT / f"{name}_citing.tsv", "w") as f:
            f.write("src\tid\tdoi\tyear\ttype\ttitle\n")
            for r in uniq:
                f.write(f"{r['src']}\t{r['id']}\t{r['doi']}\t{r['year']}\t{r['type']}\t{r['title']}\n")
        hits = [r for r in uniq if TERMS.search(r["title"] + " " + r["abstract"])]
        for r in hits:
            m = TERMS.search(r["title"] + " " + r["abstract"])
            ctx = (r["title"] + " | " + r["abstract"])
            i = ctx.lower().find(m.group(0).lower())
            screened.append({"seed": name, **{k: r[k] for k in ("year", "type", "doi", "title")},
                             "term": m.group(0), "context": ctx[max(0, i - 120): i + 160].replace("\t", " ")})
        pre = sum(1 for r in uniq if r["type"] in ("preprint", "posted-content", "PPR"))
        recent = sum(1 for r in uniq if str(r["year"]) in ("2025", "2026"))
        n_e = sum(r["src"] == "EPMC" for r in rows)
        n_o = sum(r["src"] == "OpenAlex" for r in rows)
        print(f"{name:18s} {n_e:5d} {n_o:8d} {len(uniq):6d} {pre:9d} {recent:7d} {len(hits):9d}")
    with open(OUT / "screened.tsv", "w") as f:
        f.write("seed\tyear\ttype\tdoi\ttitle\tterm\tcontext\n")
        for s in screened:
            f.write("\t".join(str(s[k]) for k in ("seed", "year", "type", "doi", "title", "term", "context")) + "\n")
    print(f"\nscreened (term in title/abstract): {len(screened)} rows -> lit/forward/screened.tsv")


if __name__ == "__main__":
    main()
