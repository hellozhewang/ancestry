#!/home/zzwang/g25/.venv/bin/python3
"""WeGene analysis — parse wegeneV16.csv (admixture by component), build
family-vs-Han-baseline tables, and analyze the relatives geographic
distribution from the *_relatives.txt files."""
import csv
import re
from pathlib import Path
from collections import defaultdict

WEGENE = Path('/home/zzwang/wegene')

# --- Parse the wegeneV16.csv ----------------------------------------------
rows = []
with open(WEGENE / 'wegeneV16.csv') as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows.append(row)

# Identify component columns: everything after the metadata columns.
META_COLS = {'Name', 'Province', 'City (Pinyin)', 'Datasource', 'Y-DNA', 'MT-DNA'}
component_cols = [c for c in rows[0].keys() if c not in META_COLS]

# Continental / aggregate columns (the last block of WeGene's output).
AGGREGATE_COLS = {'chinese_nation', 'ne_asian', 'se_asian', 'south_asian',
                  'central_asian', 'middle_eastern', 'african', 'european',
                  'american', 'oceanian'}
ETHNIC_COLS = [c for c in component_cols if c not in AGGREGATE_COLS]

# Index by name for quick lookup
by_name = {r['Name']: r for r in rows}

# Family + reference set we care about
FAMILY = ['Zhe Wang (me)', 'Yi Bo Wang (dad)', 'Ji Hu (mom)']
REFERENCE_AVERAGES = [r['Name'] for r in rows if r['Datasource'] == 'ethnicity']
INDIVIDUALS = [r['Name'] for r in rows if r['Datasource'] == 'individual']


def f(row, col):
    """Coerce a CSV cell to float, treating '' or missing as 0."""
    v = row.get(col, '')
    return float(v) if v else 0.0


def nonzero_components(row, threshold=0.5):
    """Return [(component, value), ...] of components > threshold, sorted desc."""
    pairs = [(c, f(row, c)) for c in ETHNIC_COLS]
    pairs = [(c, v) for c, v in pairs if v >= threshold]
    pairs.sort(key=lambda p: -p[1])
    return pairs


print(f"Parsed {len(rows)} WeGene rows: {len(REFERENCE_AVERAGES)} ethnic "
      f"averages + {len(INDIVIDUALS)} individuals.")
print(f"{len(component_cols)} component columns, {len(ETHNIC_COLS)} ethnic + "
      f"{len(AGGREGATE_COLS)} continental.")
print()

# === FAMILY DECOMPOSITION ================================================
print("="*70)
print("FAMILY ADMIXTURE PROFILES (WeGene)")
print("="*70)
for name in FAMILY:
    if name not in by_name:
        print(f"  [{name} not found]")
        continue
    r = by_name[name]
    pairs = nonzero_components(r, threshold=0.3)
    han_total = sum(f(r, c) for c in ('han_northern', 'han_southern'))
    ne_total = sum(f(r, c) for c in ('mongolian', 'korean', 'japanese',
                                     'tungus', 'yakut'))
    se_total = sum(f(r, c) for c in ('dai', 'kinh', 'cambodian', 'thai',
                                     'gaoshan', 'she', 'lahu', 'miao_yao'))
    # Continental aggregates
    chinese_nation = f(r, 'chinese_nation')
    ne_asian = f(r, 'ne_asian')
    se_asian = f(r, 'se_asian')
    print(f"\n[{name}] Y={r['Y-DNA']} MT={r['MT-DNA']}  "
          f"Province={r['Province']} City={r['City (Pinyin)']}")
    print(f"  Han_North={f(r,'han_northern'):.2f}  Han_South={f(r,'han_southern'):.2f}  "
          f"(combined Han={han_total:.2f}%)")
    print(f"  NE-cluster: mongolian={f(r,'mongolian'):.2f}  korean={f(r,'korean'):.2f}  "
          f"japanese={f(r,'japanese'):.2f}  tungus={f(r,'tungus'):.2f}  "
          f"(NE total={ne_total:.2f}%)")
    print(f"  SE-cluster: dai={f(r,'dai'):.2f}  she={f(r,'she'):.2f}  "
          f"miao_yao={f(r,'miao_yao'):.2f}  kinh={f(r,'kinh'):.2f}  "
          f"gaoshan={f(r,'gaoshan'):.2f}  (SE total={se_total:.2f}%)")
    print(f"  Continental: chinese_nation={chinese_nation:.2f}%  "
          f"ne_asian={ne_asian:.2f}%  se_asian={se_asian:.2f}%")
    print(f"  All non-zero ethnic components (≥0.3%):")
    for c, v in pairs:
        print(f"     {c:>20}: {v:.2f}")

# === COMPARE FAMILY TO HAN AVERAGE ========================================
print("\n" + "="*70)
print("FAMILY vs Han_Average (deviations from baseline Han)")
print("="*70)
han_avg = by_name.get('Han Average')
if han_avg:
    for name in FAMILY:
        r = by_name[name]
        print(f"\n[{name}] component-level deltas vs Han Average:")
        deltas = [(c, f(r, c) - f(han_avg, c)) for c in ETHNIC_COLS]
        deltas.sort(key=lambda p: -abs(p[1]))
        for c, d in deltas[:10]:
            sign = '+' if d > 0 else ''
            print(f"  {c:>20}: {sign}{d:.2f}%   "
                  f"(target={f(r,c):.2f}, han_avg={f(han_avg,c):.2f})")

# === COMPARE EACH FAMILY MEMBER TO EVERY ETHNIC AVERAGE ==================
# Compute "distance" between an individual and an ethnic average using
# component-level Euclidean distance (treating each WeGene component as a
# dimension). Lower = closer.
import math

def admix_distance(row_a, row_b):
    """Euclidean distance between two rows on ethnic components."""
    s = 0.0
    for c in ETHNIC_COLS:
        s += (f(row_a, c) - f(row_b, c)) ** 2
    return math.sqrt(s)

print("\n" + "="*70)
print("FAMILY → CLOSEST ETHNIC AVERAGES (Euclidean on WeGene components)")
print("="*70)
for name in FAMILY:
    r = by_name[name]
    dists = [(eth, admix_distance(r, by_name[eth]))
             for eth in REFERENCE_AVERAGES]
    dists.sort(key=lambda p: p[1])
    print(f"\n[{name}] top 10 closest ethnic averages:")
    for eth, d in dists[:10]:
        print(f"  {d:6.2f}  {eth}")

# === FAMILY-TO-FAMILY DISTANCES ===========================================
print("\n" + "="*70)
print("FAMILY-MEMBER DISTANCES (WeGene component space)")
print("="*70)
for i, a in enumerate(FAMILY):
    for b in FAMILY[i+1:]:
        d = admix_distance(by_name[a], by_name[b])
        print(f"  {a:<25} ↔ {b:<25}  d={d:.2f}")

# Also compare to wife and Korean friend if present
for extra in ('Cecilia Wang (wife)', 'Korean friend'):
    if extra in by_name:
        print(f"\n  Distances from {extra}:")
        for fam in FAMILY:
            d = admix_distance(by_name[extra], by_name[fam])
            print(f"    ↔ {fam}: d={d:.2f}")

# === RELATIVES GEOGRAPHIC DISTRIBUTION ====================================
def parse_relatives(path):
    """Parse a relatives.txt file: alternating province / count lines.

    Format observed: blank-line separated pairs:
        浙江
        100
    """
    pairs = []
    lines = [l.strip() for l in open(path).readlines() if l.strip()]
    i = 0
    while i + 1 < len(lines):
        prov = lines[i]
        cnt = lines[i+1]
        try:
            pairs.append((prov, int(cnt)))
        except ValueError:
            i += 1
            continue
        i += 2
    return pairs


print("\n" + "="*70)
print("RELATIVES GEOGRAPHIC DISTRIBUTION (匹配亲属省份分布)")
print("="*70)
for fam_name, fname in [('zhe', 'zhe_relatives.txt'),
                         ('dad', 'dad_relatives.txt'),
                         ('mom', 'mom_relatives.txt')]:
    p = WEGENE / fname
    if not p.exists():
        continue
    rels = parse_relatives(p)
    total = sum(c for _, c in rels)
    print(f"\n[{fam_name}] {len(rels)} provinces, {total} total relatives:")
    # Group provinces into regions
    REGION_GROUP = {
        'East': ['浙江', '江苏', '上海', '安徽', '福建', '江西', '山东'],
        'North': ['北京', '天津', '河北', '山西', '内蒙古', '河南'],
        'NorthEast': ['辽宁', '吉林', '黑龙江'],
        'Central': ['湖北', '湖南'],
        'South': ['广东', '广西', '海南', '香港'],
        'SouthWest': ['四川', '贵州', '云南', '重庆', '西藏'],
        'NorthWest': ['陕西', '甘肃', '青海', '宁夏', '新疆'],
        'Other': ['台湾', '籍贯地未知'],
    }
    region_totals = defaultdict(int)
    for prov, cnt in rels:
        for region, provs in REGION_GROUP.items():
            if prov in provs:
                region_totals[region] += cnt
                break
        else:
            region_totals['?Unknown'] += cnt
    print(f"  By region:")
    for region in ['East', 'North', 'Central', 'South', 'SouthWest',
                   'NorthEast', 'NorthWest', 'Other', '?Unknown']:
        if region_totals[region]:
            pct = region_totals[region] / total * 100
            print(f"    {region:<11} {region_totals[region]:>4}  ({pct:5.1f}%)")
    print(f"  Top 8 provinces:")
    for prov, cnt in rels[:8]:
        pct = cnt / total * 100
        print(f"    {prov:<10} {cnt:>4}  ({pct:5.1f}%)")
