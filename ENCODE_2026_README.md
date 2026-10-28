# ENCODE 2026 replicate-aware downloader

This setup inventories and downloads the WTC11 and H1 ENCODE experiments listed
in `encode_2026_experiments.csv`. The default destination is:

```text
/gpfs/commons/home/atalukder/RNA_Splicing/data/ENCODE_2026
```

The downloader keeps LR and SR separate, then separates experiments and biological
replicates. It never assigns pooled or replicate-unresolved processed files to an
individual replicate.

## Stages

Inventory metadata without downloading data files:

```bash
python3 code/encode_2026_downloader.py inventory
```

Review the global `manifest.csv` and each experiment's
`experiment_file_manifest.csv`. Each experiment also gets an
`experiment_summary.txt` containing the ENCODE Summary and Attribution fields,
library details, platforms, tags, and its replicate table. Then download selected
files:

```bash
python3 -u code/encode_2026_downloader.py download --workers 4 \
  2>&1 | tee -a data/ENCODE_2026/download.log
```

The download mode rebuilds the inventory first. Existing files with a valid MD5
are not downloaded again, and interrupted downloads remain as `.part` files for
resume. A detailed result and overall progress block is printed after every file.
`python3 -u` makes output immediate, and `tee` both displays it and saves it in
`data/ENCODE_2026/download.log`.

Verify all selected files later:

```bash
python3 code/encode_2026_downloader.py verify
```

Show current download progress without changing any files:

```bash
python3 code/encode_2026_downloader.py status
```

For an automatically refreshed display, press `Ctrl-C` when you want to stop
watching (the downloader will continue):

```bash
watch -n 30 python3 code/encode_2026_downloader.py status
```

## File selection

LR selects released FASTQ/FASTA `reads`, TSV/text `transcript quantifications`,
and GTF/GFF `transcriptome annotations`. SR selects released FASTQ `reads` and
TSV/text `transcript quantifications`. Every associated ENCODE file—including
excluded files—is recorded in the experiment CSV with its selection reason.

The inventory should be reviewed before running `download`, especially rows marked
`pooled_replicates` or `no_biological_replicate` and assembly/annotation columns
that do not match the intended GRCh38/GENCODE V29 analysis reference.
