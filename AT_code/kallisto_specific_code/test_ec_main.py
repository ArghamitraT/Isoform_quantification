from JOLI import Expec_Max
import numpy as np
import pandas as pd
import os

KALLISTO_DIR = "/gpfs/commons/home/tmehta/proj/kallisto_lr/ds52"

OUT_DIR = "/gpfs/commons/home/tmehta/proj/JOLI/Results/ec_debug"

os.makedirs(OUT_DIR, exist_ok=True)

em = Expec_Max(
    file_names=[],
    alignment_file="",
    count_file=OUT_DIR + "/",
    GD_lr=0.01,
    alpha_initial=1,
    max_em_rounds=0,
    load=0,
    load_filename="",
    experiment_num=1,
    dirichlet_builtin=0,
    EM_type="MAP",
    process='theta',
    use_kallisto_ec=True,
    kallisto_dir=KALLISTO_DIR,
    debug_max_ec=200,
)

sample_key = "sample1"

# pull arrays out of the model
ec_labels   = em.all_readName[sample_key]
tx_names    = em.all_readIsoformMap[sample_key]
compat_list = em.all_Yri[sample_key]
p_nm        = em.all_read_iso_prob[sample_key]
weights     = em.all_read_weights[sample_key]
Z_vals      = em.all_Phi_ri[sample_key]
n_vec       = em.all_n[sample_key]
theta_names = em.theta_names[sample_key]
theta_vals  = em.all_theta[sample_key]

_, ec_idx = np.unique(ec_labels, return_inverse=True)

row_df = pd.DataFrame({
    "row_idx":        np.arange(len(ec_labels)),
    "ec_label":       ec_labels,
    "ec_index":       ec_idx,
    "transcript":     tx_names,
    "compat":         compat_list,
    "p_nm":           p_nm,
    "weight":         weights,
    "Z":              Z_vals,
    "Z_times_weight": Z_vals * weights,
})

row_path = os.path.join(OUT_DIR, "ec_rows_init.tsv")
row_df.to_csv(row_path, sep="\t", index=False)
print(f"Wrote row-level debug to: {row_path}")

theta_df = pd.DataFrame({
    "transcript": theta_names,
    "theta_init": theta_vals,
    "n_init":     n_vec,
})
theta_path = os.path.join(OUT_DIR, "theta_n_init.tsv")
theta_df.to_csv(theta_path, sep="\t", index=False)
print(f"Wrote theta/n debug to: {theta_path}")