"""Tests for FileIO.pyx — parse_flag, generate_extra_info."""

import pytest
import pandas as pd
import numpy as np

# Flag constants matching cluster_utils.h
CLT_PASS = 1
CLT_ASSEMBLED = 2
CLT_LEFT_FLANK_MAP = 4
CLT_RIGHT_FLANK_MAP = 8
CLT_DIFF_FLANK_MAP = 16
CLT_SAME_FLANK_MAP = 32
CLT_TE_MAP = 64
CLT_POLYA = 128
CLT_TSD = 256
CLT_5P_FULL = 512
CLT_3P_FULL = 1024
CLT_5P_UNKNOWN = 2048
CLT_3P_UNKNOWN = 4096
CLT_SINGLE_TE = 8192
CLT_SELF_TO_SELF = 16384
CLT_LINE = 32768
CLT_SINE = 65536
CLT_RETROPOSON = 131072
CLT_LTR = 262144
CLT_DNA = 524288
CLT_SOLO_LTR = 8388608


def test_parse_flag_full_te():
    """Test parse_flag with a full-length LINE insertion (both flank mapped)."""
    from LOCATE.FileIO import parse_flag

    # DIFF_FLANK_MAP | SAME_FLANK_MAP yields "both_end"
    flag = CLT_PASS | CLT_ASSEMBLED | CLT_POLYA | CLT_TSD | \
           CLT_5P_FULL | CLT_3P_FULL | CLT_LINE | \
           CLT_DIFF_FLANK_MAP | CLT_SAME_FLANK_MAP
    df = pd.DataFrame({"flag": [flag], "frequency": [0.95]})
    parse_flag(df, genotyper='threshold')

    assert df["passed"].iloc[0] == True
    assert df["assembled"].iloc[0] == True
    assert df["has_polya"].iloc[0] == True
    assert df["has_tsd"].iloc[0] == True
    assert df["reconstructed_ends"].iloc[0] == "both_end"
    assert df["truncation"].iloc[0] == "full"
    assert df["te_class"].iloc[0] == "LINE"
    assert df["genotype"].iloc[0] == "1/1"


def test_parse_flag_truncated():
    """Test parse_flag with a 5' truncated, 3' full insertion."""
    from LOCATE.FileIO import parse_flag

    flag = CLT_PASS | CLT_3P_FULL | CLT_DNA | CLT_LEFT_FLANK_MAP
    df = pd.DataFrame({"flag": [flag], "frequency": [0.1]})
    parse_flag(df, genotyper='threshold')

    assert df["truncation"].iloc[0] == "5p_truncated"
    assert df["te_class"].iloc[0] == "DNA"
    assert df["genotype"].iloc[0] == "0/0"
    assert df["reconstructed_ends"].iloc[0] == "only_left"


def test_parse_flag_solo_ltr():
    """Test parse_flag with a solo LTR insertion."""
    from LOCATE.FileIO import parse_flag

    # only DIFF_FLANK_MAP -> "both_end"; LEFT_FLANK_MAP would take priority
    flag = CLT_PASS | CLT_ASSEMBLED | CLT_LTR | CLT_SOLO_LTR | \
           CLT_SINGLE_TE | CLT_DIFF_FLANK_MAP
    df = pd.DataFrame({"flag": [flag], "frequency": [0.5]})
    parse_flag(df, genotyper='threshold')

    assert df["solo_ltr"].iloc[0] == True
    assert df["singleton"].iloc[0] == True
    assert df["te_class"].iloc[0] == "LTR"
    assert df["genotype"].iloc[0] == "0/1"
    assert df["reconstructed_ends"].iloc[0] == "both_end"


def test_parse_flag_empty():
    """Test parse_flag with an unannotated flag."""
    from LOCATE.FileIO import parse_flag

    df = pd.DataFrame({"flag": [0], "frequency": [0.0]})
    parse_flag(df, genotyper='threshold')

    assert df["passed"].iloc[0] == False
    assert df["assembled"].iloc[0] == False
    assert df["reconstructed_ends"].iloc[0] == "unknown"
    assert df["truncation"].iloc[0] == "5p3p_truncated"
    assert df["te_class"].iloc[0] == "unknown"
    assert df["genotype"].iloc[0] == "0/0"


def test_generate_extra_info():
    """Test generate_extra_info format string with fully parsed row."""
    from LOCATE.FileIO import parse_flag, generate_extra_info

    flag = CLT_PASS | CLT_ASSEMBLED | CLT_LTR
    df = pd.DataFrame({
        "flag": [flag], "frequency": [0.95],
        "insertion_id": ["0-1"],
    })
    parse_flag(df, genotyper='threshold')
    # Add required columns for generate_extra_info
    df["leftclip_reads"] = 10
    df["spanning_reads"] = 5
    df["rightclip_reads"] = 8
    df["probability"] = 0.95

    row = df.iloc[0]
    result = generate_extra_info(row)
    assert "reconstructedEnds=unknown" in result
    assert "teClass=LTR" in result


def test_parse_flag_bayesian_hom_ref():
    """Bayesian genotyper: all ref reads → 0/0 with quality score."""
    from LOCATE.FileIO import parse_flag

    flag = CLT_PASS | CLT_ASSEMBLED
    df = pd.DataFrame({
        "flag": [flag], "frequency": [0.0],
        "leftclip_reads": [0], "spanning_reads": [0],
        "rightclip_reads": [0], "num_ref": [50],
    })
    parse_flag(df)

    assert df["genotype"].iloc[0] == "0/0"
    assert df["genotype_quality"].iloc[0] >= 30


def test_parse_flag_bayesian_hom_alt():
    """Bayesian genotyper: mostly alt reads → 1/1 with quality score."""
    from LOCATE.FileIO import parse_flag

    flag = CLT_PASS | CLT_ASSEMBLED
    df = pd.DataFrame({
        "flag": [flag], "frequency": [0.95],
        "leftclip_reads": [40], "spanning_reads": [5],
        "rightclip_reads": [30], "num_ref": [1],
    })
    parse_flag(df)

    assert df["genotype"].iloc[0] == "1/1"
    assert df["genotype_quality"].iloc[0] >= 20


def test_parse_flag_bayesian_het():
    """Bayesian genotyper: balanced alt/ref → 0/1 with quality score."""
    from LOCATE.FileIO import parse_flag

    flag = CLT_PASS | CLT_ASSEMBLED
    df = pd.DataFrame({
        "flag": [flag], "frequency": [0.4],
        "leftclip_reads": [10], "spanning_reads": [5],
        "rightclip_reads": [8], "num_ref": [30],
    })
    parse_flag(df)

    assert df["genotype"].iloc[0] == "0/1"
    assert "genotype_quality" in df.columns


def test_parse_flag_bayesian_genotype_quality_column():
    """Bayesian genotyper should add genotype_quality column."""
    from LOCATE.FileIO import parse_flag

    flag = CLT_PASS
    df = pd.DataFrame({
        "flag": [flag], "frequency": [0.5],
        "leftclip_reads": [5], "spanning_reads": [2],
        "rightclip_reads": [5], "num_ref": [10],
    })
    parse_flag(df)
    assert "genotype_quality" in df.columns
    assert df["genotype_quality"].iloc[0] >= 0


def test_parse_flag_threshold_no_genotype_quality():
    """Threshold genotyper should NOT add genotype_quality column."""
    from LOCATE.FileIO import parse_flag

    flag = CLT_PASS
    df = pd.DataFrame({"flag": [flag], "frequency": [0.5]})
    parse_flag(df, genotyper='threshold')
    assert "genotype_quality" not in df.columns
