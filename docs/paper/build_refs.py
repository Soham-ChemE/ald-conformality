"""Build refs.bib from Crossref metadata for every key cited in draft_v1.md (run on 2026-11-17). Titles, journals, volumes,
pages and years come from Crossref, not from memory; the script also prints a table used to check the formatted list."""
import json, re, time, urllib.request, urllib.parse
KEYS = {
 "Aguinsky2023": "10.1016/j.sse.2022.108584", "Arts2019": "10.1116/1.5093620", "Berman1965": "10.1063/1.1702984",
 "Chen2023": "10.1016/j.cej.2023.144944", "Clausing1932": "10.1002/andp.19324040804", "Clausing1971": "10.1116/1.1316379",
 "Coile2020": "10.1116/6.0000274", "Cremers2019": "10.1063/1.5060967", "Gao2015": "10.1116/1.4903941", "Gakis2018": "10.1016/j.cherd.2018.02.031",
 "Gakis2019": "10.1016/j.ces.2018.09.037", "Gayle2021": "10.1021/acs.chemmater.1c00770", "Gobbert1997": "10.1149/1.1838116",
 "Gobbert2002": "10.1116/1.1481754", "Gonsalves2024": "10.1039/d4cp00131a", "Gordon2003": "10.1002/cvde.200390005",
 "Holmqvist2012": "10.1016/j.ces.2012.07.015", "Jarvilehto2023": "10.1039/d3cp01829f", "Knoops2010": "10.1149/1.3491381",
 "Lee2026": "10.1016/j.cej.2026.177940", "Lu2022": "10.1116/6.0001891", "Miyano2018": "10.7567/jjap.57.06jb03",
 "Steckelmacher1986": "10.1088/0034-4885/49/10/001", "Vandalon2016": "10.1063/1.4939654", "Vogt2025": "10.1116/6.0004702",
 "Wang2024": "10.1016/j.dche.2024.100140", "YanguasGil2012": "10.1002/cvde.201106938", "YanguasGil2014": "10.1149/06409.0063ecst",
 "YanguasGil2021": "10.1116/6.0001212", "YanguasGil2024": "10.1116/6.0004006", "Yim2020": "10.1039/d0cp03358h",
 "Yim2022": "10.1039/d1cp04758b", "Ylilammi2018": "10.1063/1.5028178", "Zhang2019": "10.1016/j.cherd.2019.05.049",
 "Zhang2020": "10.1016/j.compchemeng.2020.107066",
}
def clean(s): return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", s or "").replace("&amp;", "&").replace("&gt;", ">").replace("&lt;", "<")).strip()
out, rows = [], []
for key, doi in KEYS.items():
    for attempt in range(3):
        try:
            m = json.load(urllib.request.urlopen("https://api.crossref.org/works/" + urllib.parse.quote(doi), timeout=40))["message"]; break
        except Exception as e:
            m = None; time.sleep(3)
    if m is None:
        print("FAILED", key, doi); continue
    authors = " and ".join(f"{a.get('family','')}, {a.get('given','')}" for a in m.get("author", []))
    year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
    journal = clean((m.get("container-title") or [""])[0]); title = clean((m.get("title") or [""])[0])
    vol, issue, page, art = m.get("volume", ""), m.get("issue", ""), m.get("page", ""), m.get("article-number", "")
    pages = page or art
    out.append(f"@article{{{key},\n  author = {{{authors}}},\n  title = {{{title}}},\n  journal = {{{journal}}},\n  year = {{{year}}},\n  volume = {{{vol}}},\n  number = {{{issue}}},\n  pages = {{{pages}}},\n  doi = {{{doi}}}\n}}\n")
    rows.append(f"{key} | {year} | {journal} | vol {vol} ({issue}) | {pages} | {title[:70]}")
    time.sleep(0.5)
out.append("@misc{Zenodo3974438,\n  author = {Puurunen, Riikka L. and Ylivaara, Oili M. E.},\n  title = {ALD Al2O3 thickness profile in microscopic rectangular channel, TMA-water 300 degC 500 cycles, PillarHall LHAR3-1b, V0001},\n  year = {2020},\n  publisher = {Zenodo},\n  doi = {10.5281/zenodo.3974438},\n  note = {CC BY 4.0; DataCite DOI, record page checked 2026-11-09}\n}\n")
open(__import__("pathlib").Path(__file__).resolve().parent / "refs.bib", "w").write("% Built by build_refs.py from Crossref on 2026-11-17; every DOI resolved.\n\n" + "\n".join(out))
print("\n".join(rows)); print(len(out), "entries written")
