import numpy as np
import pytest
from scipy.stats import spearmanr

from plm_interp.baseline.scoring import masked_logprobs, score_mutants, window_for
from plm_interp.data.proteingym import dms_info, load_dms, parse_mutant, published_benchmark, scoring_region, shift_mutants


def test_parse_mutant():
    assert parse_mutant("A12G") == [("A", 11, "G")]
    assert parse_mutant("A12G:L30P") == [("A", 11, "G"), ("L", 29, "P")]
    with pytest.raises(ValueError):
        parse_mutant("12G")


def test_window_short_sequence_is_whole():
    assert window_for(5, 100) == (0, 100)


def test_window_long_sequence_contains_position():
    for pos in (0, 600, 1999):
        s, e = window_for(pos, 2000)
        assert e - s == 1022 and s <= pos < e


def test_logprobs_are_normalised_over_20_aa(bundle, seq):
    lp = masked_logprobs(bundle, seq, [0, 10, 40])
    assert lp.shape == (3, 20)
    # 20 amino acids carry nearly all the mass (special tokens get almost none)
    assert np.all(np.exp(lp).sum(axis=1).between(0.95, 1.0001))


def test_scores_are_additive(bundle, seq):
    lp = masked_logprobs(bundle, seq, [4, 9])
    a, b, ab = score_mutants([f"{seq[4]}5A", f"{seq[9]}10W", f"{seq[4]}5A:{seq[9]}10W"], lp, seq)
    assert ab == pytest.approx(a + b)


def test_wrong_wildtype_raises(bundle, seq):
    lp = masked_logprobs(bundle, seq, [0])
    with pytest.raises(ValueError):
        score_mutants(["A1G"], lp, seq)  # ubiquitin starts with M


def test_reproduces_published_proteingym_spearman(bundle):
    dms_id = "RASH_HUMAN_Bandaru_2017"
    seq = dms_info(dms_id).target_seq
    df = load_dms(dms_id)
    lp = masked_logprobs(bundle, seq)
    rho = spearmanr(score_mutants(df.mutant, lp, seq), df.DMS_score)[0]
    published = published_benchmark().set_index("DMS_id").loc[dms_id, "ESM2 (8M)"]
    assert rho == pytest.approx(published, abs=1e-3)


def test_shift_mutants_renumbers():
    assert shift_mutants(["A12G", "A12G:L30P"], 10) == ["A2G", "A2G:L20P"]
    assert shift_mutants(["A12G"], 0) == ["A12G"]


def test_scoring_region_follows_proteingym():
    # Short protein: full sequence
    seq, off = scoring_region("RASH_HUMAN_Bandaru_2017")
    assert off == 0 and seq == dms_info("RASH_HUMAN_Bandaru_2017").target_seq
    # Long protein with a partial MSA: scored on the MSA region only, and the mutated residues line up
    dms_id = "KCNH2_HUMAN_Kozek_2020"
    seq, off = scoring_region(dms_id)
    assert off == 534 and len(seq) == 31
    for w, p, _ in (t for m in shift_mutants(load_dms(dms_id).mutant, off) for t in parse_mutant(m)):
        assert seq[p] == w
