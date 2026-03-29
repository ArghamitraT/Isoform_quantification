import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv(
    "/gpfs/commons/home/tmehta/proj/JOLI/Results/ec_debug/theta_compare.tsv",
    sep="\t"
)

eps = 1e-12
x = df["theta_kallisto"].values + eps
y = df["theta_joli"].values + eps

plt.figure(figsize=(5,5))
plt.loglog(x, y, ".", alpha=0.3, markersize=4)
plt.xlabel("kallisto θ (normalized TPM)")
plt.ylabel("JOLI θ")
plt.title("JOLI vs kallisto on ds_52_sub_2000")

lims = [
    np.min([x.min(), y.min()]),
    np.max([x.max(), y.max()])
]
plt.plot(lims, lims, linestyle="--")
plt.xlim(lims)
plt.ylim(lims)

plt.tight_layout()
plt.savefig("joli_vs_kallisto_theta_scatter.png", dpi=300)
plt.show()
