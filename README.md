# fastQC_summary_stats
A python script which pulls basic summary stats pre- and post-filtering from `fastQC` output files.

## Basic usage
```bash
python3 fastQC_summary_stats.py \
    --pre pre_trim \
    --post post_trim \
    --output fastqc_summary.csv
```

## Expected directory structure
```bash
project/ 
├── pre_trim/ 
│ ├── SAMPLE1_R1_fastqc.zip 
│ ├── SAMPLE1_R2_fastqc.zip 
│ ├── SAMPLE2_R1_fastqc.zip 
│ └── SAMPLE2_R2_fastqc.zip 
│ 
└── post_trim/ 
  ├── SAMPLE1_R1_fastqc.zip 
  ├── SAMPLE1_R2_fastqc.zip 
  ├── SAMPLE2_R1_fastqc.zip 
  └── SAMPLE2_R2_fastqc.zip
```

## Read-pair naming expectations

The script is deliberately NOT dependent on a particular or rigid FastQC filename convention, although it does require the `R1` and `R2` tokens to identify read-pairs.

```bash
It identifies R1/R2 using the read-pair token, e.g.:

    SAMPLE_R1_001_fastqc.zip
    SAMPLE_R2_001_fastqc.zip

    SAMPLE_R1_001__post_trim.fastqc.zip
    SAMPLE_R2_001__post_trim.fastqc.zip

    SAMPLE_R1_fastqc.zip
    SAMPLE_R2_fastqc.zip

```
The biological sample ID is everything **before** the `R1`/`R2` token.

```bash
For example:

    NMR2-COL1-HYP-24_R1_001_post_trim.fastqc.zip
                       ^
                       |
                   R1 token

becomes:

    sample = NMR2-COL1-HYP-24
    pair   = R1
```

**Note**. FastQC ZIP files are inspected directly. They do not need to be unzipped.

## Output metrics

```bash
Output:

    sample
    original_read_count
    original_quality_metric
    original_gc_percent
    post_trim_read_count
    reads_maintained_percent
    post_trim_gc_percent
    post_trim_quality_metric
    original_sequence_length
    post_trim_sequence_length
```

Things to **note** here are that the `quality_metric` is a weighted mean of the "Per sequence quality scores" distribution, and that R1 and R2 are combined at the sample level in the output `csv` file.

## Requirements

```bash
Python >= 3.9
```

**Note**. No external Python packages are required.
