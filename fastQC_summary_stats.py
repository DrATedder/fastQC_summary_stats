#!/usr/bin/env python3

"""
Robust paired-end FastQC pre/post-trimming summary.

The script is deliberately NOT dependent on a particular FastQC
filename convention.

It identifies R1/R2 using the read-pair token, e.g.:

    SAMPLE_R1_001_fastqc.zip
    SAMPLE_R2_001_fastqc.zip

    SAMPLE_R1_001__post_trim.fastqc.zip
    SAMPLE_R2_001__post_trim.fastqc.zip

    SAMPLE_R1_fastqc.zip
    SAMPLE_R2_fastqc.zip

The biological sample ID is everything BEFORE the R1/R2 token.

For example:

    NMR2-COL1-HYP-24_R1_001_post_trim.fastqc.zip
                       ^
                       |
                   R1 token

becomes:

    sample = NMR2-COL1-HYP-24
    pair   = R1

FastQC ZIP files are inspected directly. They do not need to be
unzipped.

The script also supports already-unzipped FastQC directories.

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

Quality metric:

    Weighted mean of the "Per sequence quality scores" distribution.

R1 and R2 are combined at the sample level.

Requirements:

    Python >= 3.8

No external Python packages are required.
"""

import argparse
import csv
import re
import zipfile
from pathlib import Path


# =====================================================================
# Output columns
# =====================================================================

OUTPUT_COLUMNS = [
    "sample",
    "original_read_count",
    "original_quality_metric",
    "original_gc_percent",
    "post_trim_read_count",
    "reads_maintained_percent",
    "post_trim_gc_percent",
    "post_trim_quality_metric",
    "original_sequence_length",
    "post_trim_sequence_length",
]


# =====================================================================
# Sample / R1 / R2 identification
# =====================================================================

def identify_sample_and_pair(name):
    """
    Identify the biological sample and whether the file is R1 or R2.

    The important assumption is that the filename contains a standard
    paired-end token such as:

        _R1_
        _R2_
        _R1.
        _R2.
        _R1
        _R2

    Everything BEFORE that token is treated as the biological sample.

    Examples
    --------

    NMR2-COL1-HYP-24_R1_001_fastqc.zip

        sample = NMR2-COL1-HYP-24
        pair   = R1

    NMR2-COL1-HYP-24_R2_001__post_trim.fastqc.zip

        sample = NMR2-COL1-HYP-24
        pair   = R2

    SAMPLE_A_R1_fastqc.zip

        sample = SAMPLE_A
        pair   = R1

    The suffix after R1/R2 is deliberately ignored.
    """

    name = Path(name).name

    # -------------------------------------------------------------
    # Remove common FastQC/archive endings.
    #
    # This is deliberately permissive because we don't want the
    # parser to depend on "_fastqc.zip".
    # -------------------------------------------------------------

    cleaned = re.sub(
        r"\.zip$",
        "",
        name,
        flags=re.IGNORECASE
    )

    cleaned = re.sub(
        r"\.fastqc$",
        "",
        cleaned,
        flags=re.IGNORECASE
    )

    cleaned = re.sub(
        r"_fastqc$",
        "",
        cleaned,
        flags=re.IGNORECASE
    )

    # -------------------------------------------------------------
    # Look for R1/R2.
    #
    # Require R1/R2 to be separated from the surrounding name.
    #
    # This prevents something like:
    #
    #     SAMPLE_R10
    #
    # being incorrectly interpreted as R1.
    #
    # Accepted:
    #
    #     _R1_
    #     _R1.
    #     _R1-
    #     _R1
    #
    # Same for R2.
    # -------------------------------------------------------------

    match = re.search(
        r"(?P<separator>^|[_\-.])(?P<pair>R[12])(?=$|[_\-.])",
        cleaned,
        flags=re.IGNORECASE
    )

    if match is None:
        return None, None

    pair = match.group("pair").upper()

    # Everything before the R1/R2 token is the sample.
    #
    # Example:
    #
    # NMR2-COL1-HYP-24_R1_001__post_trim
    #
    # becomes:
    #
    # NMR2-COL1-HYP-24

    sample = cleaned[:match.start("pair")]

    # Remove the separator immediately preceding R1/R2.
    sample = sample.rstrip("_-.")

    if not sample:
        return None, None

    return sample, pair


# =====================================================================
# FastQC ZIP detection
# =====================================================================

def is_fastqc_zip(path):
    """
    Determine whether a ZIP file is actually a FastQC archive.

    Rather than relying on the filename, inspect the ZIP contents
    for fastqc_data.txt.
    """

    try:

        with zipfile.ZipFile(path, "r") as archive:

            return any(
                name.endswith("fastqc_data.txt")
                for name in archive.namelist()
            )

    except (
        zipfile.BadZipFile,
        OSError,
    ):

        return False


# =====================================================================
# Find fastqc_data.txt inside ZIP
# =====================================================================

def find_fastqc_data_in_zip(archive):
    """
    Return the path to fastqc_data.txt inside a FastQC ZIP.
    """

    candidates = [
        name
        for name in archive.namelist()
        if name.endswith("fastqc_data.txt")
    ]

    if not candidates:
        return None

    return candidates[0]


# =====================================================================
# Parse FastQC text
# =====================================================================

def parse_fastqc_stream(handle):
    """
    Parse a FastQC fastqc_data.txt stream.
    """

    metrics = {
        "total_sequences": None,
        "sequence_length": None,
        "gc_percent": None,
        "quality_distribution": [],
    }

    current_module = None

    for raw_line in handle:

        if isinstance(raw_line, bytes):

            line = raw_line.decode(
                "utf-8",
                errors="replace"
            ).rstrip("\r\n")

        else:

            line = raw_line.rstrip("\r\n")

        # -------------------------------------------------------------
        # Basic FastQC summary
        # -------------------------------------------------------------

        if line.startswith("Total Sequences"):

            parts = line.split("\t")

            if len(parts) >= 2:

                try:
                    metrics["total_sequences"] = int(
                        parts[1].replace(",", "")
                    )
                except ValueError:
                    pass

        elif line.startswith("Sequence length"):

            parts = line.split("\t")

            if len(parts) >= 2:
                metrics["sequence_length"] = parts[1]

        elif line.startswith("%GC"):

            parts = line.split("\t")

            if len(parts) >= 2:

                try:
                    metrics["gc_percent"] = float(
                        parts[1]
                    )
                except ValueError:
                    pass

        # -------------------------------------------------------------
        # FastQC module start
        # -------------------------------------------------------------

        elif line.startswith(">>"):

            if line.startswith(">>END_MODULE"):

                current_module = None

            else:

                current_module = line[2:].split("\t")[0]

        # -------------------------------------------------------------
        # Per sequence quality scores
        # -------------------------------------------------------------

        elif (
            current_module == "Per sequence quality scores"
            and line
            and not line.startswith("#")
        ):

            parts = line.split("\t")

            if len(parts) < 2:
                continue

            try:

                quality = float(parts[0])
                count = float(parts[1])

            except ValueError:

                continue

            metrics["quality_distribution"].append(
                (quality, count)
            )

    return metrics


# =====================================================================
# Parse a FastQC ZIP
# =====================================================================

def parse_fastqc_zip(path):

    with zipfile.ZipFile(path, "r") as archive:

        fastqc_data = find_fastqc_data_in_zip(
            archive
        )

        if fastqc_data is None:

            raise RuntimeError(
                f"No fastqc_data.txt found in {path}"
            )

        with archive.open(
            fastqc_data,
            "r"
        ) as handle:

            return parse_fastqc_stream(
                handle
            )


# =====================================================================
# Parse an unzipped FastQC directory
# =====================================================================

def parse_fastqc_directory(path):

    fastqc_data = path / "fastqc_data.txt"

    if not fastqc_data.exists():

        raise RuntimeError(
            f"Could not find fastqc_data.txt in {path}"
        )

    with open(
        fastqc_data,
        "r",
        encoding="utf-8",
        errors="replace"
    ) as handle:

        return parse_fastqc_stream(
            handle
        )


# =====================================================================
# Find FastQC reports
# =====================================================================

def find_fastqc_reports(directory):
    """
    Search recursively for FastQC reports.

    ZIP files are identified by CONTENT rather than filename.

    Unzipped reports are identified by the presence of:

        fastqc_data.txt
    """

    directory = Path(directory)

    if not directory.exists():

        raise FileNotFoundError(
            f"Directory does not exist:\n{directory}"
        )

    reports = []

    # -----------------------------------------------------------------
    # ZIP files
    # -----------------------------------------------------------------

    print("  Searching ZIP files...")

    zip_files = sorted(
        directory.rglob("*.zip")
    )

    for path in zip_files:

        if not is_fastqc_zip(path):
            continue

        sample, pair = identify_sample_and_pair(
            path.name
        )

        if sample is None:

            print(
                f"  WARNING: FastQC ZIP found but "
                f"could not identify R1/R2:\n"
                f"           {path}"
            )

            continue

        reports.append(
            {
                "sample": sample,
                "pair": pair,
                "path": path,
                "type": "zip",
            }
        )

    # -----------------------------------------------------------------
    # Unzipped FastQC directories
    # -----------------------------------------------------------------

    print("  Searching unzipped FastQC reports...")

    fastqc_data_files = sorted(
        directory.rglob("fastqc_data.txt")
    )

    for fastqc_data in fastqc_data_files:

        fastqc_dir = fastqc_data.parent

        sample, pair = identify_sample_and_pair(
            fastqc_dir.name
        )

        if sample is None:

            print(
                f"  WARNING: fastqc_data.txt found but "
                f"could not identify R1/R2:\n"
                f"           {fastqc_dir}"
            )

            continue

        reports.append(
            {
                "sample": sample,
                "pair": pair,
                "path": fastqc_dir,
                "type": "directory",
            }
        )

    return reports


# =====================================================================
# Calculate weighted mean quality
# =====================================================================

def calculate_mean_quality(distribution):

    if not distribution:
        return None

    total_reads = sum(
        count
        for _, count in distribution
    )

    if total_reads <= 0:
        return None

    weighted_sum = sum(
        quality * count
        for quality, count in distribution
    )

    return weighted_sum / total_reads


# =====================================================================
# Parse all FastQC reports for one stage
# =====================================================================

def summarise_directory(directory):

    print()
    print(
        f"Scanning:\n  {directory}"
    )

    reports = find_fastqc_reports(
        directory
    )

    if not reports:

        raise RuntimeError(
            f"\nNo FastQC reports could be identified in:\n"
            f"{directory}\n\n"
            f"Check that the directory contains FastQC ZIP "
            f"files or unzipped FastQC directories."
        )

    print(
        f"  Identified {len(reports)} FastQC reports."
    )

    samples = {}

    # -----------------------------------------------------------------
    # Parse individual R1/R2 reports
    # -----------------------------------------------------------------

    for report in reports:

        sample = report["sample"]
        pair = report["pair"]

        if sample not in samples:
            samples[sample] = {}

        if pair in samples[sample]:

            print(
                f"  WARNING: Duplicate {pair} report for "
                f"{sample}:"
            )

            print(
                f"           {report['path']}"
            )

            continue

        # -------------------------------------------------------------
        # Parse report
        # -------------------------------------------------------------

        try:

            if report["type"] == "zip":

                metrics = parse_fastqc_zip(
                    report["path"]
                )

            else:

                metrics = parse_fastqc_directory(
                    report["path"]
                )

        except Exception as error:

            print(
                f"  ERROR reading {report['path']}:"
            )

            print(
                f"        {error}"
            )

            continue

        quality = calculate_mean_quality(
            metrics["quality_distribution"]
        )

        samples[sample][pair] = {
            "read_count":
                metrics["total_sequences"],

            "quality":
                quality,

            "gc":
                metrics["gc_percent"],

            "length":
                metrics["sequence_length"],
        }

        print(
            f"  {sample:<35} {pair}"
        )

    # -----------------------------------------------------------------
    # Combine R1 and R2
    # -----------------------------------------------------------------

    combined = {}

    for sample in sorted(samples):

        reads = samples[sample]

        r1 = reads.get("R1")
        r2 = reads.get("R2")

        if r1 is None:

            print(
                f"  WARNING: {sample} is missing R1"
            )

        if r2 is None:

            print(
                f"  WARNING: {sample} is missing R2"
            )

        # -------------------------------------------------------------
        # Read counts
        # -------------------------------------------------------------

        r1_count = (
            r1["read_count"]
            if (
                r1 is not None
                and r1["read_count"] is not None
            )
            else 0
        )

        r2_count = (
            r2["read_count"]
            if (
                r2 is not None
                and r2["read_count"] is not None
            )
            else 0
        )

        total_reads = (
            r1_count + r2_count
        )

        # -------------------------------------------------------------
        # Quality
        # -------------------------------------------------------------

        quality_values = []

        if (
            r1 is not None
            and r1["quality"] is not None
        ):

            quality_values.append(
                (
                    r1["quality"],
                    r1_count
                )
            )

        if (
            r2 is not None
            and r2["quality"] is not None
        ):

            quality_values.append(
                (
                    r2["quality"],
                    r2_count
                )
            )

        if quality_values:

            denominator = sum(
                count
                for _, count
                in quality_values
            )

            combined_quality = (
                sum(
                    quality * count
                    for quality, count
                    in quality_values
                )
                / denominator
            )

        else:

            combined_quality = None

        # -------------------------------------------------------------
        # GC
        # -------------------------------------------------------------

        gc_values = []

        if (
            r1 is not None
            and r1["gc"] is not None
        ):

            gc_values.append(
                (
                    r1["gc"],
                    r1_count
                )
            )

        if (
            r2 is not None
            and r2["gc"] is not None
        ):

            gc_values.append(
                (
                    r2["gc"],
                    r2_count
                )
            )

        if gc_values:

            denominator = sum(
                count
                for _, count
                in gc_values
            )

            combined_gc = (
                sum(
                    gc * count
                    for gc, count
                    in gc_values
                )
                / denominator
            )

        else:

            combined_gc = None

        # -------------------------------------------------------------
        # Store
        # -------------------------------------------------------------

        combined[sample] = {
            "read_count":
                total_reads,

            "quality":
                combined_quality,

            "gc":
                combined_gc,

            "r1_length":
                r1["length"]
                if r1 is not None
                else None,

            "r2_length":
                r2["length"]
                if r2 is not None
                else None,
        }

    return combined


# =====================================================================
# Format sequence lengths
# =====================================================================

def format_lengths(r1, r2):

    if r1 is None and r2 is None:
        return None

    r1 = (
        "NA"
        if r1 is None
        else str(r1)
    )

    r2 = (
        "NA"
        if r2 is None
        else str(r2)
    )

    return f"{r1}/{r2}"


# =====================================================================
# Combine pre- and post-trimming
# =====================================================================

def combine_results(pre, post):

    all_samples = sorted(
        set(pre.keys()) |
        set(post.keys())
    )

    results = []

    for sample in all_samples:

        before = pre.get(
            sample,
            {}
        )

        after = post.get(
            sample,
            {}
        )

        original_reads = before.get(
            "read_count"
        )

        post_trim_reads = after.get(
            "read_count"
        )

        # -------------------------------------------------------------
        # Reads maintained
        # -------------------------------------------------------------

        if (
            original_reads is not None
            and post_trim_reads is not None
            and original_reads > 0
        ):

            reads_maintained = (
                post_trim_reads
                / original_reads
                * 100
            )

        else:

            reads_maintained = None

        # -------------------------------------------------------------
        # Lengths
        # -------------------------------------------------------------

        original_length = format_lengths(
            before.get("r1_length"),
            before.get("r2_length")
        )

        post_trim_length = format_lengths(
            after.get("r1_length"),
            after.get("r2_length")
        )

        # -------------------------------------------------------------
        # Row
        # -------------------------------------------------------------

        row = {

            "sample":
                sample,

            "original_read_count":
                original_reads,

            "original_quality_metric":
                round(
                    before["quality"],
                    3
                )
                if before.get("quality") is not None
                else None,

            "original_gc_percent":
                round(
                    before["gc"],
                    2
                )
                if before.get("gc") is not None
                else None,

            "post_trim_read_count":
                post_trim_reads,

            "reads_maintained_percent":
                round(
                    reads_maintained,
                    2
                )
                if reads_maintained is not None
                else None,

            "post_trim_gc_percent":
                round(
                    after["gc"],
                    2
                )
                if after.get("gc") is not None
                else None,

            "post_trim_quality_metric":
                round(
                    after["quality"],
                    3
                )
                if after.get("quality") is not None
                else None,

            "original_sequence_length":
                original_length,

            "post_trim_sequence_length":
                post_trim_length,
        }

        results.append(row)

    return results


# =====================================================================
# Write CSV
# =====================================================================

def write_csv(records, output_file):

    with open(
        output_file,
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS
        )

        writer.writeheader()
        writer.writerows(records)


# =====================================================================
# Main
# =====================================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Generate a robust paired-end FastQC "
            "pre/post-trimming summary."
        )
    )

    parser.add_argument(
        "--pre",
        required=True,
        help="Pre-trimming FastQC directory"
    )

    parser.add_argument(
        "--post",
        required=True,
        help="Post-trimming FastQC directory"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file"
    )

    args = parser.parse_args()

    print("=" * 70)
    print("FASTQC PRE/POST-TRIMMING SUMMARY")
    print("=" * 70)

    # -----------------------------------------------------------------
    # Pre-trimming
    # -----------------------------------------------------------------

    print("\n[1/4] PRE-TRIMMING")

    pre = summarise_directory(
        args.pre
    )

    print(
        f"\nPre-trimming samples identified: "
        f"{len(pre)}"
    )

    # -----------------------------------------------------------------
    # Post-trimming
    # -----------------------------------------------------------------

    print("\n[2/4] POST-TRIMMING")

    post = summarise_directory(
        args.post
    )

    print(
        f"\nPost-trimming samples identified: "
        f"{len(post)}"
    )

    # -----------------------------------------------------------------
    # Compare sample sets
    # -----------------------------------------------------------------

    print("\n[3/4] COMPARING SAMPLE SETS")

    pre_only = sorted(
        set(pre) - set(post)
    )

    post_only = sorted(
        set(post) - set(pre)
    )

    common = sorted(
        set(pre) & set(post)
    )

    print(
        f"  Samples in both:       {len(common)}"
    )

    print(
        f"  Pre-trim only:          {len(pre_only)}"
    )

    print(
        f"  Post-trim only:         {len(post_only)}"
    )

    if pre_only:

        print("\n  WARNING: samples only in pre-trimming:")

        for sample in pre_only:
            print(
                f"    {sample}"
            )

    if post_only:

        print("\n  WARNING: samples only in post-trimming:")

        for sample in post_only:
            print(
                f"    {sample}"
            )

    # -----------------------------------------------------------------
    # Combine
    # -----------------------------------------------------------------

    results = combine_results(
        pre,
        post
    )

    # -----------------------------------------------------------------
    # Write
    # -----------------------------------------------------------------

    print("\n[4/4] WRITING OUTPUT")

    write_csv(
        results,
        args.output
    )

    print(
        f"\nOutput:"
        f"\n  {args.output}"
    )

    print(
        f"\nRows written: {len(results)}"
    )

    print("\nDone.")
    print("=" * 70)


if __name__ == "__main__":
    main()
