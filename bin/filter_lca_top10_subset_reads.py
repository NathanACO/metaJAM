#!/usr/bin/env python3
import argparse
import csv
import random
import sys
import re
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


@dataclass
class TaxNode:
    taxid: int
    name: str
    rank: str


@dataclass
class LcaRecord:
    read_id: str
    seq: str
    length: int
    count: int  # collapsed multiplicity
    ranks: Dict[str, TaxNode]  # e.g. {"genus": TaxNode(...), "family": ...}


def parse_first_field(field: str) -> Tuple[str, str, int, int]:
    """
    Parse:
      <read_id>:<SEQ>:<LEN>:<COUNT>
    read_id itself may contain ':' so we split from the right.
    """
    parts = field.rstrip("\n").rsplit(":", 3)
    if len(parts) != 4:
        raise ValueError(f"Cannot parse first field into read_id/seq/len/count: {field}")
    read_id, seq, length_s, count_s = parts
    return read_id, seq, int(length_s), int(count_s)


def parse_tax_chain(fields: List[str]) -> Dict[str, TaxNode]:
    """
    Parse tab-separated taxonomy tokens like:
      9606:Homo sapiens:species
    Return dict keyed by rank.
    """
    ranks: Dict[str, TaxNode] = {}
    for tok in fields:
        tok = tok.strip()
        if not tok:
            continue
        # taxid:name:rank (name may contain spaces, but ':' separators are fixed)
        # Split only twice.
        p = tok.split(":", 2)
        if len(p) != 3:
            continue
        taxid_s, name, rank = p
        try:
            taxid = int(taxid_s)
        except ValueError:
            continue
        ranks[rank] = TaxNode(taxid=taxid, name=name, rank=rank)
    return ranks


def parse_lca_line(line: str) -> Optional[LcaRecord]:
    line = line.rstrip("\n")
    if not line:
        return None
    # LCA file appears tab-separated: first field + many tax tokens
    fields = line.split("\t")
    if len(fields) < 2:
        return None
    read_id, seq, length, count = parse_first_field(fields[0])
    ranks = parse_tax_chain(fields[1:])
    return LcaRecord(read_id=read_id, seq=seq, length=length, count=count, ranks=ranks)


def pick_kingdom_rank(ranks: Dict[str, TaxNode]) -> Tuple[int, str, str]:
    """
    Your LCA has e.g. Metazoa:kingdom, Viridiplantae:kingdom, etc.
    But many lineages might not have 'kingdom' consistently. We fall back to 'superkingdom'.
    Returns: (taxid, name, rank_used)
    """
    if "kingdom" in ranks:
        k = ranks["kingdom"]
        return k.taxid, k.name, "kingdom"
    if "superkingdom" in ranks:
        k = ranks["superkingdom"]
        return k.taxid, k.name, "superkingdom"
    return 0, "unknown", "unknown"


def get_rank(ranks: Dict[str, TaxNode], rank: str) -> Tuple[int, str]:
    if rank in ranks:
        t = ranks[rank]
        return t.taxid, t.name
    return 0, "NA"



def load_bamdam_mean_damage_by_taxid(path: str) -> dict[int, float]:
    """Load bamdam TSV and compute mean damage = (Damage+1 + Damage-1)/2 per TaxNodeID."""
    dmg: dict[int, float] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        # bamdam TSV has a header; values may be quoted
        reader = csv.DictReader((ln for ln in fh if ln.strip() and not ln.startswith("#")), delimiter="\t")
        if reader.fieldnames is None:
            return dmg
        for row in reader:
            taxid_s = (row.get("TaxNodeID") or "").strip()
            if not taxid_s:
                continue
            try:
                taxid = int(taxid_s)
            except ValueError:
                continue
            try:
                d1 = float((row.get("Damage+1") or "0").strip().strip('"'))
                d2 = float((row.get("Damage-1") or "0").strip().strip('"'))
            except ValueError:
                d1, d2 = 0.0, 0.0
            dmg[taxid] = (d1 + d2) / 2.0
    return dmg

def sample_reads_without_expanding(entries: List[LcaRecord], target_n: int, rng: random.Random) -> List[LcaRecord]:
    """
    Sample up to target_n reads accounting for record.count without expanding huge lists.
    We sample *without replacement* from the total multiplicity mass:
      - each draw selects an entry proportional to remaining count
      - decrement that entry's remaining count
    Returns a list (length <= target_n) of LcaRecord references (may repeat same record multiple times).
    """
    remaining = [max(0, e.count) for e in entries]
    total = sum(remaining)
    if total <= 0 or target_n <= 0:
        return []
    out: List[LcaRecord] = []

    # For up to 100 draws this is totally fine O(N*draws)
    for _ in range(min(target_n, total)):
        r = rng.randint(1, total)
        cum = 0
        chosen_idx = None
        for i, c in enumerate(remaining):
            if c <= 0:
                continue
            cum += c
            if cum >= r:
                chosen_idx = i
                break
        if chosen_idx is None:
            break
        out.append(entries[chosen_idx])
        remaining[chosen_idx] -= 1
        total -= 1

    return out


def main():
    ap = argparse.ArgumentParser(
        description="From a .lca file: pick top N genera per kingdom, sample reads per genus, output combined FASTA + mapping TSV."
    )
    ap.add_argument("--lca", required=True, help="Input .lca file")
    ap.add_argument("--top-genera", type=int, default=10, help="Top genera per kingdom (default: 10)")
    ap.add_argument("--genera-list", default="", help="Comma/space-separated list of genera to investigate (taxid or name). If set, overrides --top-genera.")
    ap.add_argument("--genera-file", default="", help="File with one genus per line (name or taxid). Lines starting with # are ignored. If set, overrides --top-genera (can be combined with --genera-list).")
    ap.add_argument("--max-reads", type=int, default=100, help="Max reads to sample per genus (default: 100)")
    ap.add_argument("--min-reads", type=int, default=30, help="Min reads required per genus (default: 30). Below: skipped.")
    ap.add_argument("--bamdam-tsv", default=None, help="Bamdam per-sample TSV (used to filter genera by damage).")
    ap.add_argument("--min-dmg", type=float, default=0.0, help="Minimum mean damage required to keep a genus. Uses mean(Damage+1, Damage-1). If >1, treated as percent and divided by 100.")
    ap.add_argument("--seed", type=int, default=1, help="Random seed (default: 1), if we set up 42, the same reads will be picked up")
    ap.add_argument("--out-fasta", required=True, help="Output combined FASTA")
    ap.add_argument("--out-map", required=True, help="Output mapping TSV (expected taxonomy per query), links each MMseqs2 query sequence you create back to: 1) which original read it came from; 2) which expected taxonomy (from the LCA) it should belong to (genus / family / order / kingdom taxids + names)")
    args = ap.parse_args()

    # Optional genus damage filter (bamdam TSV)
    min_dmg = float(getattr(args, "min_dmg", 0.0) or 0.0)
    # Allow users to specify percent (e.g. 3.5) or proportion (e.g. 0.035)
    if min_dmg > 1.0:
        min_dmg = min_dmg / 100.0

    dmg_by_taxid: dict[int, float] = {}
    if args.bamdam_tsv and min_dmg > 0.0:
        try:
            dmg_by_taxid = load_bamdam_mean_damage_by_taxid(args.bamdam_tsv)
        except Exception as e:
            print(f"[WARN] Could not load bamdam TSV '{args.bamdam_tsv}': {e}", file=sys.stderr)
            dmg_by_taxid = {}

    rng = random.Random(args.seed)

    # 1) read all records
    records: List[LcaRecord] = []
    with open(args.lca, "r", encoding="utf-8") as fh:
        for ln, line in enumerate(fh, 1):
            try:
                rec = parse_lca_line(line)
            except Exception as e:
                print(f"[WARN] Skipping unparsable line {ln}: {e}", file=sys.stderr)
                continue
            if rec is None:
                continue
            # We need genus to do anything meaningful
            genus_taxid, genus_name = get_rank(rec.ranks, "genus")
            if genus_taxid == 0:
                continue
            records.append(rec)

    if not records:
        print("[ERROR] No usable records with genus rank found.", file=sys.stderr)
        sys.exit(2)

    # 2) abundance per (kingdom, genus) using multiplicity counts
    genus_abund: Dict[Tuple[int, int], int] = defaultdict(int)  # (kingdom_taxid, genus_taxid) -> sum(count)
    genus_meta: Dict[Tuple[int, int], Tuple[str, str, str]] = {}  # -> (kingdom_name, genus_name, kingdom_rank_used)

    for rec in records:
        k_taxid, k_name, k_rank_used = pick_kingdom_rank(rec.ranks)
        g_taxid, g_name = get_rank(rec.ranks, "genus")
        key = (k_taxid, g_taxid)
        genus_abund[key] += max(0, rec.count)
        genus_meta[key] = (k_name, g_name, k_rank_used)
    # 3) select genera (either user-provided list or top genera per kingdom)
    # If --genera-list is provided, it overrides --top-genera.
    genera_list_raw = (args.genera_list or "").strip()
    genera_file = (args.genera_file or "").strip()
    if genera_file:
        try:
            file_items = []
            with open(genera_file, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    line = line.strip()
                    if (not line) or line.startswith("#"):
                        continue
                    file_items.append(line)
            if file_items:
                if genera_list_raw:
                    genera_list_raw = genera_list_raw + "," + ",".join(file_items)
                else:
                    genera_list_raw = ",".join(file_items)
        except FileNotFoundError:
            raise SystemExit(f"ERROR: --genera-file not found: {genera_file}")
    wanted_taxids = set()
    wanted_names = set()
    if genera_list_raw:
        # split on commas and/or whitespace
        for tok in [t for t in re.split(r"[\s,]+", genera_list_raw) if t]:
            if tok.isdigit():
                wanted_taxids.add(int(tok))
            else:
                wanted_names.add(tok.strip().strip('"').strip("'").lower())

    selected: List[Tuple[int, int]] = []  # list of (kingdom_taxid, genus_taxid)

    if wanted_taxids or wanted_names:
        for key, meta in genus_meta.items():
            _k_name, g_name, _k_rank_used = meta
            if key[1] in wanted_taxids or (g_name and g_name.lower() in wanted_names):
                selected.append(key)
    else:
        # default: top genera per kingdom by abundance
        by_kingdom: Dict[int, List[Tuple[int, int]]] = defaultdict(list)  # kingdom_taxid -> list of (genus_taxid, abundance)
        for (k_taxid, g_taxid), abund in genus_abund.items():
            by_kingdom[k_taxid].append((g_taxid, abund))
        for k_taxid, lst in by_kingdom.items():
            lst_sorted = sorted(lst, key=lambda x: x[1], reverse=True)
            top = lst_sorted[: max(0, args.top_genera)]
            selected.extend([(k_taxid, g_taxid) for g_taxid, _ in top])

    selected_set = set(selected)
    if not selected_set:
        if wanted_taxids or wanted_names:
            print("[ERROR] Nothing selected (check --genera-list values).", file=sys.stderr)
        else:
            print("[ERROR] Nothing selected (check --top-genera).", file=sys.stderr)
        sys.exit(2)


    # 4) taxon_group records by (kingdom, genus)
    taxon_group: Dict[Tuple[int, int], List[LcaRecord]] = defaultdict(list)
    for rec in records:
        k_taxid, _, _ = pick_kingdom_rank(rec.ranks)
        g_taxid, _ = get_rank(rec.ranks, "genus")
        key = (k_taxid, g_taxid)
        if key in selected_set:
            taxon_group[key].append(rec)

    # 5) sample and write outputs
    total_written = 0
    skipped_low_support = 0

    with open(args.out_fasta, "w", encoding="utf-8") as fa, open(args.out_map, "w", newline="", encoding="utf-8") as mp:
        w = csv.writer(mp, delimiter="\t")
        w.writerow([
            "query_id",
            "orig_read_id",
            "rep_index",
            "exp_kingdom_taxid",
            "exp_kingdom_name",
            "exp_genus_taxid",
            "exp_genus_name",
            "exp_family_taxid",
            "exp_family_name",
            "exp_order_taxid",
            "exp_order_name",
        ])

        for key in sorted(taxon_group.keys(), key=lambda k: (str(genus_meta.get(k, ("", "", ""))[0]), genus_meta.get(k, ("", "", ""))[1])):
            k_taxid, g_taxid = key
            entries = taxon_group[key]
            # total multiplicity mass in taxon_group
            total_reads_in_group = sum(max(0, e.count) for e in entries)
            if total_reads_in_group < args.min_reads:
                skipped_low_support += 1
                continue

            sampled = sample_reads_without_expanding(entries, args.max_reads, rng)
            if len(sampled) < args.min_reads:
                skipped_low_support += 1
                continue

            # metadata (names)
            k_name, g_name, _k_rank_used = genus_meta.get(key, ("unknown", "NA", "unknown"))

            # For expected family/order, use first sampled record (all records in genus should share those ranks typically)
            # If not, it’s still fine; we store expected per query from each read’s LCA ranks below.
            rep_counter: Counter[str] = Counter()

            for rec in sampled:
                rep_counter[rec.read_id] += 1
                rep_i = rep_counter[rec.read_id]
                query_id = f"{rec.read_id}__rep{rep_i}"

                fam_taxid, fam_name = get_rank(rec.ranks, "family")
                ord_taxid, ord_name = get_rank(rec.ranks, "order")

                # FASTA header includes expected taxids for convenience
                fa.write(f">{query_id}|exp_genus={g_taxid}|exp_family={fam_taxid}|exp_order={ord_taxid}|exp_kingdom={k_taxid}\n")
                fa.write(rec.seq + "\n")

                w.writerow([
                    query_id,
                    rec.read_id,
                    rep_i,
                    k_taxid,
                    k_name,
                    g_taxid,
                    g_name,
                    fam_taxid,
                    fam_name,
                    ord_taxid,
                    ord_name,
                ])

                total_written += 1

    print(f"[INFO] Selected genus: {len(taxon_group)}", file=sys.stderr)
    print(f"[INFO] Number of selected sequences to run through MMSeqs2: {total_written}", file=sys.stderr)
    print(f"[INFO] Skipped genus due to not enough reads: {skipped_low_support}", file=sys.stderr)


if __name__ == "__main__":
    main()
