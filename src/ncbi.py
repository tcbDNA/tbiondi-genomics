"""Minimal rate-limited NCBI E-utilities helper (stdlib only; no contact details are sent)."""
import json, time, urllib.parse, urllib.request, urllib.error
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
_last = [0.0]
def call(tool, _tries=6, **params):
    data = urllib.parse.urlencode(params).encode()
    for i in range(_tries):
        wait = 0.4 - (time.time() - _last[0])
        if wait > 0: time.sleep(wait)
        _last[0] = time.time()
        try:
            with urllib.request.urlopen(urllib.request.Request(BASE + tool + ".fcgi", data=data), timeout=180) as r:
                return r.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError, ConnectionError) as e:
            err = e; time.sleep(2 * (i + 1))
    raise RuntimeError(f"NCBI {tool} failed after {_tries} tries: {err}")
def search(db, term, retmax=100):
    d = json.loads(call("esearch", db=db, term=term, retmode="json", retmax=retmax))["esearchresult"]
    return int(d.get("count", 0)), d.get("idlist", []), d
