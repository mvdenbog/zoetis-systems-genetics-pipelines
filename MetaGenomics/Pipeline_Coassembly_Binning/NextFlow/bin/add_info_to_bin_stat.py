#!/usr/bin/env python3

"""----------------------------------------------------------------------------
  Script Name: add_info_to_bin_stat.py
  Description: Classifies metagenomic bins into quality categories (High-quality,
               Medium-quality, Low-quality, High-contamination) based on 
               completeness and contamination thresholds from CheckM/CheckM2.
               Calculates the size of unbinned contigs by subtracting the sum 
               of binned sizes from the total assembly length (from QUAST), 
               and appends a 'Not_binned' entry to the output table.
  Input files: 
    - Bins statistics TSV (with columns: genome, Size, completeness, contamination)
    - QUAST report TSV (with assembly metrics including total length)
  Output: TSV file with original bin stats plus 'Quality' column and 'Not_binned' row
  Created By:  Mathias Vandenbogaert
  Date:        2026-05-06
-------------------------------------------------------------------------------
"""

# Metadata
__author__ = 'Mathias Vandenbogaert'
__version__ = '0.1'
__email__ = 'mathias.vandenbogaert@lizard.bio'
__status__ = 'dev'


import logging
import sys
import argparse
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from pathlib import Path
from typing import Final

import pandas as pd


# ============================================================================
# CONSTANTS
# ============================================================================
DEFAULT_OUTPUT: Final[str] = "bins_stat_and_quality.tsv"
DEFAULT_MAX_CONTAMINATION: Final[int] = 10
DEFAULT_MAX_HIGH_QUAL_CONTAMINATION: Final[int] = 5
DEFAULT_MEDIUM_QUAL_COMPLETENESS: Final[int] = 50
DEFAULT_HIGH_QUAL_COMPLETENESS: Final[int] = 90

# Quality category labels
QUAL_HIGH_CONTAM: Final[str] = 'High-contamination'
QUAL_LOW: Final[str] = 'Low-quality'
QUAL_MEDIUM: Final[str] = 'Medium-quality'
QUAL_HIGH: Final[str] = 'High-quality'
QUAL_NOT_BINNED: Final[str] = 'Not-binned'

# Column names
COL_GENOME: Final[str] = 'genome'
COL_SIZE: Final[str] = 'Size'
COL_COMPLETENESS: Final[str] = 'completeness'
COL_CONTAMINATION: Final[str] = 'contamination'
COL_QUALITY: Final[str] = 'Quality'
COL_ASSEMBLY: Final[str] = 'Assembly'
COL_TOTAL_LENGTH: Final[str] = 'Total length (>= 0 bp)'


# ============================================================================
# FUNCTIONS
# ============================================================================

def parse_arguments() -> argparse.Namespace:
    """Parse and validate script command-line arguments."""
    parser = ArgumentParser(
        description="Classify metagenomic bins by quality and add unbinned contig size.",
        formatter_class=ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        '-s', '--bins_stat',
        required=True,
        type=Path,
        help="Path to bins stat file (TSV format with genome, Size, completeness, contamination)."
    )
    parser.add_argument(
        '-q', '--quast_report',
        required=True,
        type=Path,
        help="Path to QUAST report file (TSV format with assembly metrics)."
    )
    parser.add_argument(
        '-o', '--output_file',
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output table path with quality classifications and unbinned entry."
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help="Increase output verbosity."
    )

    args = parser.parse_args()
    
    # Validate input files exist
    for file_path, file_desc in [(args.bins_stat, "bins stat"), 
                                  (args.quast_report, "QUAST report")]:
        if not file_path.exists():
            logging.error(f"{file_desc} file not found: {file_path}")
            sys.exit(1)
    
    return args


def configure_logging(verbose: bool) -> None:
    """Configure logging based on verbosity flag."""
    level = logging.DEBUG if verbose else logging.WARNING
    log_format = "%(levelname)s: %(message)s"
    logging.basicConfig(level=level, format=log_format)
    if verbose:
        logging.debug('Verbose mode enabled')


def load_bins_stat(file_path: Path) -> pd.DataFrame:
    """Load and preprocess bins statistics file."""
    df = pd.read_csv(file_path, sep='\t')
    if COL_GENOME not in df.columns:
        logging.error(f"Required column '{COL_GENOME}' not found in bins stat file")
        sys.exit(1)
    df.set_index(COL_GENOME, drop=False, inplace=True)
    return df


def load_quast_assembly_size(file_path: Path) -> int:
    """Extract total assembly length from QUAST report."""
    df = pd.read_csv(file_path, sep='\t', index_col=COL_ASSEMBLY)
    if df.empty:
        logging.error("QUAST report is empty or malformed")
        sys.exit(1)
    
    first_column = df.columns[0]
    try:
        total_length = df.loc[COL_TOTAL_LENGTH, first_column]
        return int(total_length)
    except KeyError:
        logging.error(f"Required metric '{COL_TOTAL_LENGTH}' not found in QUAST report")
        sys.exit(1)


def classify_bin_quality(
    df: pd.DataFrame,
    max_contamination: int = DEFAULT_MAX_CONTAMINATION,
    max_high_qual_contamination: int = DEFAULT_MAX_HIGH_QUAL_CONTAMINATION,
    medium_qual_completeness: int = DEFAULT_MEDIUM_QUAL_COMPLETENESS,
    high_qual_completeness: int = DEFAULT_HIGH_QUAL_COMPLETENESS
) -> pd.DataFrame:
    """
    Classify metagenomic bins into quality categories based on completeness 
    and contamination thresholds (MIMAG-inspired criteria).
    
    Classification rules (applied in priority order):
    1. High-contamination: contamination >= max_contamination (default: 10%)
    2. High-quality: completeness > high_qual_completeness (default: 90%) AND 
                     contamination < max_high_qual_contamination (default: 5%)
    3. Medium-quality: completeness >= medium_qual_completeness (default: 50%) AND 
                       contamination < max_contamination
    4. Low-quality: completeness < medium_qual_completeness AND 
                    contamination < max_contamination
    """
    # Initialize Quality column
    df[COL_QUALITY] = pd.NA
    
    # Define boolean masks for each category
    is_high_contam = df[COL_CONTAMINATION] >= max_contamination
    is_low_qual = (df[COL_COMPLETENESS] < medium_qual_completeness) & ~is_high_contam
    is_medium_qual = (df[COL_COMPLETENESS] >= medium_qual_completeness) & ~is_high_contam
    is_high_qual = (df[COL_COMPLETENESS] > high_qual_completeness) & \
                   (df[COL_CONTAMINATION] < max_high_qual_contamination)
    
    # Apply classifications (order ensures high-contamination takes priority)
    df.loc[is_high_contam, COL_QUALITY] = QUAL_HIGH_CONTAM
    df.loc[is_high_qual, COL_QUALITY] = QUAL_HIGH
    df.loc[is_medium_qual, COL_QUALITY] = QUAL_MEDIUM
    df.loc[is_low_qual, COL_QUALITY] = QUAL_LOW
    
    return df


def add_not_binned_entry(
    df: pd.DataFrame,
    total_assembly_size: int,
    binned_size_sum: int
) -> pd.DataFrame:
    """
    Calculate and append a row representing unbinned contigs.
    
    The unbinned size is computed as: total_assembly_size - sum(binned_sizes).
    If the result is negative (due to rounding or input inconsistencies), 
    it is clamped to 0 with a warning.
    """
    not_binned_size = total_assembly_size - binned_size_sum
    
    if not_binned_size < 0:
        logging.warning(
            f"Calculated unbinned size is negative ({not_binned_size}). "
            "Setting to 0. Verify input files for consistency."
        )
        not_binned_size = 0
    
    not_binned_row = pd.DataFrame({
        COL_GENOME: ['Not_binned'],
        COL_SIZE: [not_binned_size],
        COL_QUALITY: [QUAL_NOT_BINNED]
    }).set_index(COL_GENOME)
    
    return pd.concat([df, not_binned_row], ignore_index=False)


def write_output(df: pd.DataFrame, output_path: Path) -> None:
    """Write the processed dataframe to TSV file."""
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, sep='\t', index=False)
        logging.info(f"Output written to: {output_path}")
    except OSError as e:
        logging.error(f"Failed to write output file: {e}")
        sys.exit(1)


def main() -> None:
    """Main execution function."""
    args = parse_arguments()
    configure_logging(args.verbose)
    
    # Load input data
    logging.debug(f"Loading bins stat from: {args.bins_stat}")
    df_bins = load_bins_stat(args.bins_stat)
    
    logging.debug(f"Loading assembly size from: {args.quast_report}")
    total_assembly_size = load_quast_assembly_size(args.quast_report)
    
    # Calculate total binned size before classification
    binned_size_sum = df_bins[COL_SIZE].sum()
    
    # Classify bin quality
    logging.debug("Classifying bins by quality thresholds")
    df_classified = classify_bin_quality(df_bins)
    
    # Add unbinned entry
    logging.debug("Adding unbinned contigs entry")
    df_final = add_not_binned_entry(df_classified, total_assembly_size, binned_size_sum)
    
    # Write output
    write_output(df_final, args.output_file)
    
    # Summary logging
    quality_counts = df_final[COL_QUALITY].value_counts()
    logging.info(f"Quality distribution: {dict(quality_counts)}")
    logging.info("Processing complete")


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == '__main__':
    main()
