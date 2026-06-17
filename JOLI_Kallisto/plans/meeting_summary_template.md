# Meeting Summary Template

This file defines the **format** for advisor meeting summaries.  
Actual summaries live at: `/gpfs/commons/home/atalukder/RNA_Splicing/files/presentations/`  
Naming convention: `YYYY_MM_DD_meetingWDavid.md` (or other advisor name).

---

## Format

Each summary file has **three main sections** (readable at a glance) followed by a **Detailed Notes** section at the bottom for everything else — full calculations, derivations, equations, and implementation specifics. The advisor reads the three main sections first; Detailed Notes are for deep-dives or reference.

```markdown
# Meeting with [Advisor Name] — YYYY-MM-DD

---

## 1. Background
*What this is about and why it matters.*

- Bullet 1: what the method/problem is and its goal
- Bullet 2: key design choices or components relevant to this meeting
- Bullet 3: what changed since last meeting, or what question motivated this work

---

## 2. Results
*What happened.*

*(To be filled after experiments.)*

---

## 3. Summary & Further Questions
*What it means and what comes next.*

- Bullet 1: interpretation — what the results tell us
- Bullet 2: leading hypothesis for unexpected findings
- Bullet 3: top priority next step
- Bullet 4: open question(s) for advisor

---

## Detailed Notes

Everything else goes here: full calculations, derivations, equations, pipeline diagrams, pseudocode, metric tables, figure references, experiment folder names, and implementation specifics.

### Background Detail
Pipeline diagram or step list. Core equations in LaTeX (E-step, M-step, MAP update, GD objective). Variable definitions.

### Experimental Detail
- Metric comparison table (Spearman, Pearson, MAE, RMSE, Bray-Curtis, TP, FP, FN)
- Figure references with full absolute paths
- Runtime
- Full experiment folder name (`exprmnt_YYYY_MM_DD__HH_MM_SS`)

### Derivations & Implementation
- Full derivations or algorithm pseudocode
- Key equations with variable definitions
- What changed in the code (file, flag, equation) and why
- Hypotheses with supporting equations
```

---

## Rules for Filling In

1. **Verdict blocks** — after any result table or analysis, add a `Verdict` block. The **entire block** (title + every bullet) must be yellow using `<span style="color:#E6A817">...</span>`. Format:
   ```
   <span style="color:#E6A817">**Verdict**</span>
   <span style="color:#E6A817">- **Finding 1**: ...</span>
   <span style="color:#E6A817">- **Finding 2**: ...</span>
   ```
   Each bullet should compare methods explicitly with numbers and percentages where available.
2. **Three sections first, always** — write the three main sections before Detailed Notes. If the advisor only reads those three sections, they should walk away with a complete picture of the meeting.
2. **Bullets, not prose** — each main section bullet should be one clear sentence.
3. **Math equations** — use LaTeX inline (`$...$`) or block (`$$...$$`). Always define variables inline on first use. All equations go in Detailed Notes, not the main sections.
4. **Figure references** — always include the full absolute path, e.g.:
   ```
   ![Caption](/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmnt_.../figures/fig.png)
   Path: /gpfs/commons/.../fig.png
   ```
5. **Metric tables** — always include Spearman, Pearson, MAE, RMSE, Bray-Curtis, TP, FP, FN. Follow these layout rules:
   - Place samples **side by side** in the same row, separated by a `\|` column (renders as a vertical divider).
   - Highlight the **best value per metric per sample** in green using `<span style="color:green">**value**</span>`.
   - Best = highest for Spearman, Pearson, TP; lowest for MAE, RMSE, Bray-Curtis, FP, FN.
   - When two methods tie, color both green.
   - Add a caption line above each table stating the direction of best (e.g. "Spearman, Pearson = highest; MAE, RMSE, BC = lowest").
   ```
   | Method | sim1 Spearman | Pearson | MAE | \| | sim2 Spearman | Pearson | MAE |
   |--------|--------------|---------|-----|---|--------------|---------|-----|
   | LK     | 0.77 | <span style="color:green">**0.97**</span> | 3.58 | \| | 0.78 | ... | ... |
   ```
6. **Experiment references** — always cite the full `exprmnt_YYYY_MM_DD__HH_MM_SS` folder name.
7. **No results before running** — leave section 2 as `*(To be filled after experiments.)*` until data exists.
