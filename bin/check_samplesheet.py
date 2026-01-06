#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""
"""From https://github.com/nf-core/tools/blob/e5ce6ce20304835bd40f102f038b7e1aadc888b2/nf_core/pipeline-template/bin/check_samplesheet.py"""
"""Modified by Bilal Sharif <bilal.bioinfo@gmail.com> and Verena Kutschera <verena.kutschera@scilifelab.se>"""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
        ".fq",
        ".fastq",
    )

    def __init__(
        self,
        sample_id_col="sample_id",
        library_id_col="library_id",
        lane_col="lane",
        sample_type_col="sample_type",
        library_type_col="library_type",
        first_col="fastq_1",
        second_col="fastq_2",
        single_col="single_end",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_id_col (str): sample name (default "sample_id").
            library_id_col (str): library_id (default "library_id").
            lane_col (str): lane number (default "lane").
            sample_type_col (str): type of the sample - ancient or modern (default "sample_type").
            library_type_col (str): type of the library - single-stranded or double-stranded (default "library_type").
            first_col (str): first (or only) FASTQ file path (default "fastq_1").
            second_col (str): second (if any) FASTQ file path (default "fastq_2").
            single_col (str): new column that will be inserted and records whether the sample contains single- or paired-end sequencing reads (default "single_end").
        """
        super().__init__(**kwargs)
        self._sample_id_col = sample_id_col
        self._library_id_col = library_id_col
        self._lane_col = lane_col
        self._sample_type_col = sample_type_col
        self._library_type_col = library_type_col
        self._first_col = first_col
        self._second_col = second_col
        self._single_col = single_col
        self._seen = set()
        self.modified = []
        self._multiple_fastq = set()

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        multiple_fastq = (row[self._sample_id_col], row[self._library_id_col], row[self._lane_col])
        if multiple_fastq in self._multiple_fastq:
            raise AssertionError(f"Duplicate FASTQ entries for the same sample/library_id/lane combination: {multiple_fastq}. it is recommended to add a unique identifier to the lane column in such scenarios.")
        self._multiple_fastq.add(multiple_fastq)

        self._validate_sample_id(row)
        self._validate_library_id(row)
        self._validate_lane(row)
        self._validate_sample_type(row)
        self._validate_library_type(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_pair(row)
        self._seen.add((row[self._sample_id_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample_id(self, row):
        """Assert that the sample_id exists and convert spaces and underscores to dashes."""
        if len(row[self._sample_id_col]) <= 0:
            raise AssertionError("A sample_id is required.")
        # Sanitize samples slightly.
        row[self._sample_id_col] = row[self._sample_id_col].replace(" ", "-")
        row[self._sample_id_col] = row[self._sample_id_col].replace("_", "-")

    def _validate_library_id(self, row):
        """Assert that the library_id ID exists."""
        if len(row[self._library_id_col]) <= 0:
            raise AssertionError("A library_id is required.")
        # Sanitize library_id IDs slightly.
        row[self._library_id_col] = row[self._library_id_col].replace(" ", "-")
        row[self._library_id_col] = row[self._library_id_col].replace("_", "-")

    def _validate_lane(self, row):
        """Assert that the lane number exists."""
        if len(row[self._lane_col]) <= 0:
            raise AssertionError("A lane number is required.")
        # Sanitize lane numbers slightly.
        row[self._lane_col] = row[self._lane_col].replace(" ", "-")
        row[self._lane_col] = row[self._lane_col].replace("_", "-")

    def _validate_sample_type(self, row):
        """Assert that the sample type exists and it only contains one of the following values: 'ancient', 'modern'."""
        if len(row[self._sample_type_col]) <= 0:
            raise AssertionError("A sample_type is required.")
        if row[self._sample_type_col] not in ["ancient", "modern"]:
            raise AssertionError("sample_type must be either 'ancient' or 'modern'.")

    def _validate_library_type(self, row):
        """Assert that the library type exists and it only contains one of the following values: 'single-stranded', 'double-stranded'."""
        if len(row[self._library_type_col]) <= 0:
            raise AssertionError("A library_type is required.")
        if row[self._library_type_col] not in ["single-stranded", "double-stranded"]:
            raise AssertionError("library_type must be either 'single-stranded' or 'double-stranded'.")

    def _validate_first(self, row):
        """Assert that the first FASTQ entry is non-empty and has the right format."""
        if len(row[self._first_col]) <= 0:
            raise AssertionError("At least the first FASTQ file is required.")
        self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the second FASTQ entry has the right format if it exists."""
        second = (row.get(self._second_col) or "").strip()
        if len(second) > 0:
            self._validate_fastq_format(second)

    def _validate_pair(self, row):
        """Assert that read pairs have the same file extension. Report pair status."""
        if row[self._first_col] and row[self._second_col]:
            row[self._single_col] = False
            first_col_suffix = Path(row[self._first_col]).suffixes[-2:]
            second_col_suffix = Path(row[self._second_col]).suffixes[-2:]
            if first_col_suffix != second_col_suffix:
                raise AssertionError("FASTQ pairs must have the same file extensions.")
        else:
            row[self._single_col] = True

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FORMATS):
            raise AssertionError(
                f"The FASTQ file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FORMATS)}"
            )

    def validate_unique_samples(self):
        """
        Assert that the combination of sample name and FASTQ filename is unique.
        """
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and FASTQ filename must be unique.")


def read_head(handle, num_lines=1):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by DNAharvester.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        Sample_1,library_1,L001,single-stranded,XYZ0123456,Illumina,/path/to/data/AB1_R1.fq.gz,/path/to/data/AB1_R2.fq.gz
        Sample_1,library_2,L001,single-stranded,XYZ0123456,Illumina,/path/to/data/AB1_R1.fq.gz,/path/to/data/AB1_R2.fq.gz
        Sample_2,library_1,L001,double-stranded,XYZ0123456,Illumina,/path/to/data/AB2_R1.fq.gz,/path/to/data/AB2_R2.fq.gz
        Sample_2,library_1,L002,double-stranded,XYZ0123456,Illumina,/path/to/data/AB2_R1.fq.gz,/path/to/data/AB2_R2.fq.gz
        Sample_3,library_1,L001,double-stranded,XYZ0123456,Illumina,/path/to/data/AB3.fq.gz,

    """
    required_columns = {"sample_id", "library_id", "lane", "sample_type", "library_type", "fastq_1", "fastq_2"}
    # See https://docs.python.org/3.9/library_id/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        diallect = sniff_format(in_handle)
        reader = csv.DictReader(in_handle, dialect=diallect)
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    header.insert(1, "single_end")
    # See https://docs.python.org/3.9/library_id/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
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
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())