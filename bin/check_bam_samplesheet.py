#!/usr/bin/env python


"""Provide a command line tool to validate and transform BAM samplesheets."""
"""Modified from check_samplesheet.py by Bilal Sharif <bilal.bioinfo@gmail.com>"""

import argparse
import csv
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


class BamRowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    def __init__(
        self,
        sample_id_col="sample_id",
        bam_col="bam",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_id_col (str): sample name (default "sample_id").
            bam_col (str): BAM file path (default "bam").
        """
        super().__init__(**kwargs)
        self._sample_id_col = sample_id_col
        self._bam_col = bam_col
        self._seen_samples = set()
        self._seen_bams = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample_id(row)
        self._validate_bam(row)
        self.modified.append(row)

    def _validate_sample_id(self, row):
        """Assert that the sample_id exists, is unique and convert spaces and underscores to dashes."""
        if len(row[self._sample_id_col]) <= 0:
            raise AssertionError("A sample_id is required.")
        # Sanitize samples slightly.
        row[self._sample_id_col] = row[self._sample_id_col].replace(" ", "-")
        row[self._sample_id_col] = row[self._sample_id_col].replace("_", "-")
        if row[self._sample_id_col] in self._seen_samples:
            raise AssertionError(f"Duplicate sample_id: {row[self._sample_id_col]}. Each sample should have only one BAM file.")
        self._seen_samples.add(row[self._sample_id_col])

    def _validate_bam(self, row):
        """Assert that the BAM entry is non-empty, has the right format and a unique file name."""
        if len(row[self._bam_col]) <= 0:
            raise AssertionError("A BAM file is required.")
        if not row[self._bam_col].endswith(".bam"):
            raise AssertionError(f"The BAM file has an unrecognized extension: {row[self._bam_col]}\nIt should be: .bam")
        # BAM files are staged together for joint variant calling, so file names must be unique
        bam_name = Path(row[self._bam_col]).name
        if bam_name in self._seen_bams:
            raise AssertionError(f"Duplicate BAM file name: {bam_name}. BAM file names must be unique, also when located in different directories.")
        self._seen_bams.add(bam_name)


def read_head(handle, num_lines=1):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """Detect the tabular format. The read position is expected to be at the beginning (index 0)."""
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_bam_samplesheet(file_in, file_out):
    """
    Check that the tabular BAM samplesheet has the structure expected by DNAharvester.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        sample_id,bam
        Sample_1,/path/to/data/Sample_1.bam
        Sample_2,/path/to/data/Sample_2.bam

    """
    required_columns = {"sample_id", "bam"}
    # See https://docs.python.org/3.9/library_id/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        diallect = sniff_format(in_handle)
        reader = csv.DictReader(in_handle, dialect=diallect)
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The BAM sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = BamRowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library_id/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular BAM samplesheet.",
        epilog="Example: python check_bam_samplesheet.py bam_samplesheet.csv bam_samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input BAM samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output BAM samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_bam_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
