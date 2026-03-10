#!/usr/bin/env python3
"""This script has been written for metaJAM pipeline, it aims to compare the expected taxonomy from bamdam with the mmseqs2 taxonomy file
It requires first mmseqs_pick_besthit_with_ambiguity.py to run to produce the nt_query.besthit.assigned.tsv file
How to use it:
python3 postprocess_mmseqs_taxonomy_compare.py \
  --expected-map  $SAMPLE_$LAST_MAP_DB_TAG/$SAMPLE.mmseqs_expected.tsv \
  --assigned      $SAMPLE_$LAST_MAP_DB_TAG/nt_query.besthit.assigned.tsv \
  --out-summary   $SAMPLE_$LAST_MAP_DB_TAG/$SAMPLE/nt_query.taxcmp.summary.tsv
"""

import argparse
import csv
from collections import defaultdict


def norm(x: str) -> str:
    if x is None:
        return ""
    x = str(x).strip().strip('"').strip("'")
    if x.upper() == "NA":
        return ""
    return x.lower()


def pick_expected_group_taxon(e):
    """
    Keep output shape the same, but choose the grouping taxon as:
      genus -> family -> order -> kingdom
    Returns:
      (taxid_string, name_string)
    """
    gt = str(e.get("gt", "0") or "0").strip()
    ft = str(e.get("ft", "0") or "0").strip()
    ot = str(e.get("ot", "0") or "0").strip()
    kt = str(e.get("kt", "0") or "0").strip()

    gn = e.get("gn_raw", "") or "NA"
    fn = e.get("fn_raw", "") or "NA"
    on = e.get("on_raw", "") or "NA"
    kn = e.get("kn_raw", "") or "NA"

    if e.get("g", "") and gt not in ("", "0", "NA"):
        return gt, gn
    if e.get("f", "") and ft not in ("", "0", "NA"):
        return ft, fn
    if e.get("o", "") and ot not in ("", "0", "NA"):
        return ot, on
    return kt, kn


def main():
    ap = argparse.ArgumentParser(
        description="Compare expected taxonomy vs MMseqs assigned taxonomy using NAMES (genus/family/order/kingdom)."
    )
    ap.add_argument("--expected-map", required=True, help="*.mmseqs_expected.tsv from filter_lca_top10_subset_reads.py")
    ap.add_argument("--assigned", required=True, help="*.besthit.assigned_names.tsv from mmseqs_pick_besthit_assigned_names.py")
    ap.add_argument("--out-summary", required=True, help="Output TSV summary")
    args = ap.parse_args()

    # Load expected map
    exp = {}
    with open(args.expected_map, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        need = {
            "query_id",
            "exp_kingdom_name", "exp_genus_name", "exp_family_name", "exp_order_name",
            "exp_kingdom_taxid", "exp_genus_taxid", "exp_family_taxid", "exp_order_taxid"
        }
        miss = need - set(r.fieldnames or [])
        if miss:
            raise SystemExit(f"Expected map missing columns: {sorted(miss)}")

        for row in r:
            q = row["query_id"]
            exp[q] = {
                "k": norm(row.get("exp_kingdom_name", "")),
                "g": norm(row.get("exp_genus_name", "")),
                "f": norm(row.get("exp_family_name", "")),
                "o": norm(row.get("exp_order_name", "")),
                "kt": row.get("exp_kingdom_taxid", "0"),
                "gt": row.get("exp_genus_taxid", "0"),
                "ft": row.get("exp_family_taxid", "0"),
                "ot": row.get("exp_order_taxid", "0"),
                "kn_raw": row.get("exp_kingdom_name", "NA"),
                "gn_raw": row.get("exp_genus_name", "NA"),
                "fn_raw": row.get("exp_family_name", "NA"),
                "on_raw": row.get("exp_order_name", "NA"),
            }

    # Load assigned
    asn = {}
    with open(args.assigned, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        need = {"query_id", "hit_taxid", "assigned_genus", "assigned_family", "assigned_order", "assigned_kingdom"}
        miss = need - set(r.fieldnames or [])
        if miss:
            raise SystemExit(f"Assigned file missing columns: {sorted(miss)}")

        for row in r:
            q = row["query_id"]
            asn[q] = {
                "taxid": row.get("hit_taxid", "0"),
                "g": norm(row.get("assigned_genus", "")),
                "f": norm(row.get("assigned_family", "")),
                "o": norm(row.get("assigned_order", "")),
                "k": norm(row.get("assigned_kingdom", "")),
            }

    def bucket(e, a):
        if a is None:
            return "no_hit"
        if (a.get("taxid", "0") in ("0", "", "NA")) and not (a["g"] or a["f"] or a["o"] or a["k"]):
            return "no_hit"
        if e["g"] and a["g"] and a["g"] == e["g"]:
            return "same_genus"
        if e["f"] and a["f"] and a["f"] == e["f"]:
            return "same_family"
        if e["o"] and a["o"] and a["o"] == e["o"]:
            return "same_order"
        if e["k"] and a["k"] and a["k"] == e["k"]:
            return "same_kingdom"
        return "other"

    stats = defaultdict(lambda: defaultdict(int))

    for q, e in exp.items():
        b = bucket(e, asn.get(q))
        grp_taxid, grp_name = pick_expected_group_taxon(e)
        key = (e["kt"], e["kn_raw"], grp_taxid, grp_name)
        stats[key]["n"] += 1
        stats[key][b] += 1

    with open(args.out_summary, "w", newline="", encoding="utf-8") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow([
            "exp_kingdom_taxid",
            "exp_kingdom_name",
            "exp_taxa",
            "exp_genus_name",
            "n_queries",
            "same_genus",
            "same_family",
            "same_order",
            "same_kingdom",
            "other",
            "no_hit",
            "pct_same_genus",
            "pct_same_family",
            "pct_same_order",
            "pct_same_kingdom",
            "pct_other",
            "pct_no_hit",
        ])

        for key in sorted(stats.keys(), key=lambda k: (k[1], k[3])):
            n = stats[key]["n"]
            sg = stats[key].get("same_genus", 0)
            sf = stats[key].get("same_family", 0)
            so = stats[key].get("same_order", 0)
            sk = stats[key].get("same_kingdom", 0)
            ot = stats[key].get("other", 0)
            nh = stats[key].get("no_hit", 0)
            pct = lambda x: 0.0 if n == 0 else 100.0 * x / n

            w.writerow([
                key[0], key[1], key[2], key[3], n,
                sg, sf, so, sk, ot, nh,
                f"{pct(sg):.2f}",
                f"{pct(sf):.2f}",
                f"{pct(so):.2f}",
                f"{pct(sk):.2f}",
                f"{pct(ot):.2f}",
                f"{pct(nh):.2f}",
            ])

    print(f"[OK] wrote {args.out_summary} (expected={len(exp)} assigned={len(asn)})")


if __name__ == "__main__":
    main()
