# QIIME2 + DADA2 Processing Workflow (Bacteria, Fungi, Eukaryote)

This document summarizes the preprocessing pipeline for converting raw FASTQ files into ASV tables using QIIME2 and DADA2. The procedure was applied to three datasets: **bacteria**, **fungi**, and **eukaryotes**.

---

## 1. Prepare Manifest Files

Create a manifest CSV file for each dataset (example: `manifest_bac.csv`, `manifest_fun.csv`, `manifest_euk.csv`) using absolute paths:

```csv
sample-id,absolute-filepath,direction
Sample01,/path/to/Sample01_R1.fastq.gz,forward
Sample01,/path/to/Sample01_R2.fastq.gz,reverse
...
```

## 2. Import FASTQ Files
```
# For bacteria
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_bac.csv \
  --output-path paired-end-demux_bac.qza \
  --input-format PairedEndFastqManifestPhred33

# For fungi
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest_fun.csv \
  --output-path paired-end-demux_fun.qza \
  --input-format PairedEndFastqManifestPhred33

```

## 3. Summarize Demultiplexed Reads
```
qiime demux summarize \
  --i-data paired-end-demux_bac.qza \
  --o-visualization demux_bac.qzv

qiime demux summarize \
  --i-data paired-end-demux_fun.qza \
  --o-visualization demux_fun.qzv

```

## 4. Subsample for Parameter Optimization (Optional)
```
qiime demux subsample-paired \
  --i-sequences paired-end-demux_bac.qza \
  --p-fraction 0.05 \
  --o-subsampled-sequences demux-subsample.qza

qiime demux summarize \
  --i-data demux-subsample.qza \
  --o-visualization demux-subsample.qzv
```
## 5. Denoising with DADA2
Bacteria (final parameters: f=260, r=180)
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux_bac.qza \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 180 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --o-table table_bac.qza \
  --o-representative-sequences rep-seqs_bac.qza \
  --o-denoising-stats stats-dada2_bac.qza \
  --p-n-threads 0
```
Fungi (final parameters: f=260, r=170)
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux_fun.qza \
  --p-trunc-len-f 260 \
  --p-trunc-len-r 170 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --o-table table_fun.qza \
  --o-representative-sequences rep-seqs_fun.qza \
  --o-denoising-stats stats-dada2_fun.qza \
  --p-n-threads 0
```

## 6. Visualization and Summary
```
# Visualize denoising stats
qiime metadata tabulate \
  --m-input-file stats-dada2_bac.qza \
  --o-visualization stats-dada2_bac.qzv

qiime metadata tabulate \
  --m-input-file stats-dada2_fun.qza \
  --o-visualization stats-dada2_fun.qzv

# Summarize feature tables
qiime feature-table summarize \
  --i-table table_bac.qza \
  --o-visualization table_bac.qzv \
  --m-sample-metadata-file metadata_bac.txt

qiime feature-table summarize \
  --i-table table_fun.qza \
  --o-visualization table_fun.qzv \
  --m-sample-metadata-file metadata_fun.txt

# Tabulate representative sequences
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_bac.qza \
  --o-visualization rep-seqs_bac.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs_fun.qza \
  --o-visualization rep-seqs_fun.qzv

```

## Note
- .qzv files can be opened using QIIME2 View
- Trimming parameters were selected based on subsampled read quality profiles.
