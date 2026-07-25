import os
import re
import numpy as np
from collections import defaultdict

KALLISTO_DIR = "/gpfs/commons/home/tmehta/proj/kallisto_lr/ds_52_sub_2000"

tx_path = os.path.join(KALLISTO_DIR, "transcripts.txt")
ec_path = os.path.join(KALLISTO_DIR, "count.ec.txt")
mtx_path = os.path.join(KALLISTO_DIR, "count.mtx")

with open(tx_path) as f:
    tx_names = [ln.strip() for ln in f if ln.strip()]

ec_to_txs = []
with open(ec_path) as f:
    for line in f:
        s = line.strip()
        if not s:
            ec_to_txs.append([])
        else:
            tokens = re.split(r"[,\s]+", s)
            ec_to_txs.append([int(x) for x in tokens if x])

from kallisto_ec_loader import read_count_mtx
nrows, ncols, nnz, rows, cols, vals = read_count_mtx(mtx_path)
ec_ids = cols - 1
counts = vals

tx_to_ecs = defaultdict(list)  # transcript -> [(ec_id, ec_size, count, other_tx_in_ec)]

for ec_id, count in zip(ec_ids, counts):
    tx_list = ec_to_txs[ec_id]
    tx_names_in_ec = [tx_names[i] for i in tx_list]
    for tx_idx in tx_list:
        tx_name = tx_names[tx_idx]
        others = [tx_names[i] for i in tx_list if i != tx_idx]
        tx_to_ecs[tx_name].append((ec_id, len(tx_list), count, others))

# Check problematic transcripts AND their EC partners
check_transcripts = [
    "GTF2IP10", "NM_001101.5",
    "CICP3", "unassigned_transcript_4819", 
    "NM_006082.3", "XM_047432837.1",
    "NM_002032.3", "NR_135154.2",
]

print("=== Full EC Membership ===\n")
for tx in check_transcripts:
    ecs = tx_to_ecs.get(tx, [])
    print(f"{tx}: appears in {len(ecs)} EC(s)")
    for ec_id, ec_size, count, others in ecs:
        print(f"  EC_{ec_id}: size={ec_size}, count={count}")
        print(f"    Partners: {others[:5]}{'...' if len(others) > 5 else ''}")
    
    # Check if any of the partners appear in other ECs
    all_partners = set()
    for _, _, _, others in ecs:
        all_partners.update(others)
    
    partners_with_other_ecs = []
    for partner in all_partners:
        partner_ecs = tx_to_ecs.get(partner, [])
        if len(partner_ecs) > 1:
            partners_with_other_ecs.append((partner, len(partner_ecs)))
    
    if partners_with_other_ecs:
        print(f"Partners appearing in OTHER ECs: {partners_with_other_ecs}")
    print()