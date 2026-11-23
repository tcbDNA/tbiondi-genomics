"""Prior-art audit: find GenBank accessions for known frameshift cases that are not annotated in RefSeq.

    uv run python code/audit/ncbi_lookup.py
Prints each query (verbatim), its hit count and up to 14 accessions with titles. NCBI E-utilities,
no API key, so it pauses between calls and retries on HTTP 429.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
QUERIES = {
    "Crassostrea virginica mt genome": 'Crassostrea virginica[Organism] AND mitochondrion[filter] AND complete genome[Title]',
    "Physarum polycephalum mt genome": 'Physarum polycephalum[Organism] AND mitochondrion[filter] AND complete[Title]',
    "Trypanosoma brucei maxicircle": 'Trypanosoma brucei[Organism] AND maxicircle[Title]',
    "Crithidia fasciculata maxicircle/coxII": 'Crithidia fasciculata[Organism] AND (maxicircle[Title] OR cytochrome oxidase subunit II[Title])',
    "Lophophysema eversa mt": 'Lophophysema eversa[Organism] AND mitochondrion[filter]',
    "Docosaccus maculatus mt": 'Docosaccus maculatus[Organism] AND mitochondrion[filter]',
    "Hexactinellida mt genomes (all GenBank)": 'Hexactinellida[Organism] AND mitochondrion[filter] AND genome[Title]',
    "Polyrhachis cytb 2004-2005": 'Polyrhachis[Organism] AND (cytochrome b[Title] OR cytb[Title]) AND 2004:2005[PDAT]',
    "Kazachstania/Nakaseomyces mt complete": '(Kazachstania[Organism] OR Nakaseomyces[Organism] OR Huiozyma[Organism] OR Monosporozyma[Organism]) AND mitochondrion[filter] AND complete[Title]',
    "Perkinsus mt complete": 'Perkinsus[Organism] AND mitochondrion[filter] AND complete[Title]',
    "Phaeodactylum tricornutum mt": 'Phaeodactylum tricornutum[Organism] AND mitochondrion[filter] AND complete genome[Title]',
}


def get(url: str) -> dict:
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "grigorovia-mito-insertions/0.1"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code != 429 or attempt == 5:
                raise
            time.sleep(3 * (attempt + 1))
    raise AssertionError


def main() -> None:
    for name, q in QUERIES.items():
        d = get(E + "esearch.fcgi?db=nuccore&retmode=json&retmax=40&term=" + urllib.parse.quote(q))
        ids, n = d["esearchresult"]["idlist"], d["esearchresult"]["count"]
        print(f"## {name}: {n} hits\n   query: {q}")
        if ids:
            time.sleep(1.1)
            s = get(E + "esummary.fcgi?db=nuccore&retmode=json&id=" + ",".join(ids))
            for i in s["result"]["uids"][:14]:
                x = s["result"][i]
                print(f"   {x['accessionversion']} ({x['slen']} bp) {x['title'][:80]}")
        time.sleep(1.1)


if __name__ == "__main__":
    main()
