"""Print plain-text sentences of a JATS XML matching a regex (verbatim apart from whitespace)."""
import re, sys, xml.etree.ElementTree as ET
t = " ".join(" ".join(ET.parse(sys.argv[1]).getroot().itertext()).split())
for s in re.split(r"(?<=[.?!])\s+(?=[A-Z(])", t):
    if re.search(sys.argv[2], s, re.I): print("-", s[:900])
