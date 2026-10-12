from fastapi.testclient import TestClient

from plm_interp.app.server import app

client = TestClient(app)
SEQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def test_index_and_protein_list():
    assert client.get("/").status_code == 200
    proteins = client.get("/api/proteins").json()
    assert {p["key"] for p in proteins} >= {"hras", "p53"}


def test_protein_detail_has_structure_and_lab_data():
    d = client.get("/api/protein/hras").json()
    assert len(d["sequence"]) == 189 and d["pdb"] and len(d["lab"]) > 3000


def test_predict_returns_20_probabilities():
    r = client.post("/api/predict", json={"sequence": SEQ, "position": 3, "model": "8M"}).json()
    assert len(r["probs"]) == 20 and r["wt"] == "F"
    assert abs(sum(p["p"] for p in r["probs"]) - 1) < 0.05
    assert r["tokens"][4]["token"] == "<mask>"


def test_bad_input_gives_readable_error():
    r = client.post("/api/predict", json={"sequence": "MKVXX" * 4, "position": 1, "model": "8M"})
    assert r.status_code == 400 and "standard amino acids" in r.json()["error"]
    r = client.post("/api/predict", json={"sequence": SEQ, "position": 500, "model": "8M"})
    assert r.status_code == 400


def test_inside_shapes():
    r = client.post("/api/inside", json={"sequence": SEQ, "position": 10, "mutant": "W", "model": "8M"}).json()
    assert len(r["layers"]) == 7 and len(r["heads"]) == 6 and len(r["heads"][0]) == 20


def test_batched_head_effects_match_one_at_a_time():
    import numpy as np
    import torch

    from plm_interp.app.engine import _head_effects, bundle
    from plm_interp.explain.ablation import head_means
    from plm_interp.models.hooks import ablate_heads

    b = bundle("8M")
    pos, wt, mt = 10, SEQ[10], "W"
    batched = _head_effects(b, SEQ, pos, wt, mt, chunk=7)
    means = head_means(b, SEQ)
    ids = b.encode(SEQ)
    ids[0, pos + 1] = b.mask_id

    def llr():
        with torch.no_grad():
            lp = b.model(input_ids=ids).logits[0, pos + 1].log_softmax(-1)
        return (lp[b.aa_ids["ACDEFGHIKLMNPQRSTVWY".index(mt)]] - lp[b.aa_ids["ACDEFGHIKLMNPQRSTVWY".index(wt)]]).item()

    clean = llr()
    for layer, head in [(0, 0), (2, 7), (5, 19)]:
        with ablate_heads(b, {layer: [head]}, means):
            assert np.isclose(llr() - clean, batched[layer, head], atol=1e-4)
