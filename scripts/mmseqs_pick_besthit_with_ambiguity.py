#!/usr/bin/env python3
"""This script has been written for metaJAM pipeline, it aims to:
- Reads the nt_query.resultDB_nucl_taxid.m8 (semicolon-separated, with header)
- Per query:
    Picks best hit = lowest evalue (tie-breaker: highest bits)
    Finds best hit from a different genus
    If that different-genus hit has an evalue 5% greater → genus is ambiguous, so output genus taxid as 0 (and keep family/order/kingdom from the best hit) writes a TSV with 1 line per query:query_id, lca_taxid, genus;family;order;kingdom, assigned_level

How to use it:
python3 mmseqs_pick_besthit_with_ambiguity.py \
  --m8 nt_query.resultDB_nucl_taxid.m8 \
  --ambig-frac 0.05 \
  --out nt_query.besthit.assigned.tsv
"""

import argparse
import csv

def main():
    ap = argparse.ArgumentParser(
        description="Pick one best MMseqs hit per query (min evalue, tie-break by max bits) and mark genus as ambiguous if best alt-genus is within +X% evalue."
    )
    ap.add_argument("--m8", required=True, help="MMseqs convertalis output with taxonomy columns (semicolon-separated)")
    ap.add_argument("--out", required=True, help="Output TSV (tab-separated)")
    ap.add_argument("--ambig-frac", type=float, default=0.05, help="Ambiguity fraction (0.05 = 5%%)")
    args = ap.parse_args()

    m8 = args.m8
    out = args.out
    ambig = args.ambig_frac

    best = {}      # q -> (e, -bits, taxid, genus, family, order, kingdom)
    best_alt = {}  # q -> best hit with different genus

    with open(m8, encoding="utf-8", errors="replace") as fh:
        r = csv.reader(fh, delimiter=";")
        h = [x.strip() for x in next(r)]
        iq, ie, ib = h.index("query"), h.index("evalue"), h.index("bits")
        it = h.index("Taxid")
        ig, ifa, io, ik = h.index("genus"), h.index("family"), h.index("order"), h.index("kingdom")

        for row in r:
            if len(row) <= max(iq, ie, ib, it, ig, ifa, io, ik):
                continue

            q = row[iq].split("|", 1)[0]
            try:
                e = float(row[ie])
            except Exception:
                continue
            try:
                bits = float(row[ib])
            except Exception:
                bits = -1.0

            taxid = row[it]
            genus, fam, order, king = row[ig], row[ifa], row[io], row[ik]
            rec = (e, -bits, taxid, genus, fam, order, king)

            if q not in best or rec < best[q]:
                old = best.get(q)
                best[q] = rec
                if old and old[3] != "NA" and genus != "NA" and old[3] != genus:
                    if q not in best_alt or old < best_alt[q]:
                        best_alt[q] = old
            else:
                bg = best[q][3]
                if bg != "NA" and genus != "NA" and genus != bg:
                    if q not in best_alt or rec < best_alt[q]:
                        best_alt[q] = rec

    with open(out, "w", newline="", encoding="utf-8") as fo:
        w = csv.writer(fo, delimiter="\t")
        w.writerow([
            "query_id", "hit_taxid",
            "assigned_genus", "assigned_family", "assigned_order", "assigned_kingdom",
            "assigned_level"
        ])

        for q, rec in best.items():
            e1, _, taxid, genus, fam, order, king = rec
            alt = best_alt.get(q)

            ambiguous = False
            if alt:
                e2 = alt[0]
                ambiguous = (e1 == 0.0 and e2 == 0.0) or (e1 > 0.0 and e2 <= e1 * (1.0 + ambig))

            if ambiguous:
                w.writerow([q, taxid, "NA", fam, order, king, "family"])
            else:
                w.writerow([q, taxid, genus, fam, order, king, "genus"])

    print(f"[OK] wrote {out} ({len(best)} queries)")

if __name__ == "__main__":
    main()