"""Download ESM-2 8M and 35M into models/<name>/ and verify they load and predict.

Usage (from the project root, with the plm-interp env active):
    python scripts/01_download_models.py
"""

from pathlib import Path

import torch
from huggingface_hub import snapshot_download
from transformers import AutoTokenizer, EsmForMaskedLM

ROOT = Path(__file__).resolve().parents[1]
MODELS_DIR = ROOT / "models"

MODELS = {
    "esm2_t6_8M_UR50D": "facebook/esm2_t6_8M_UR50D",
    "esm2_t12_35M_UR50D": "facebook/esm2_t12_35M_UR50D",
}

# Config + tokenizer files + safetensors only (skips the duplicate .bin and TF .h5 weights)
ALLOW = ["*.json", "*.txt", "*.safetensors"]

# Ubiquitin, used only as a smoke test
TEST_SEQ = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def dir_size_mb(path: Path) -> float:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / 1e6


def download(name: str, repo_id: str) -> Path:
    local_dir = MODELS_DIR / name
    local_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n[{name}] downloading {repo_id} -> {local_dir.relative_to(ROOT)}")
    snapshot_download(repo_id=repo_id, local_dir=local_dir, allow_patterns=ALLOW)

    if not list(local_dir.glob("*.safetensors")):
        print(f"[{name}] no safetensors in repo, falling back to pytorch_model.bin")
        snapshot_download(repo_id=repo_id, local_dir=local_dir, allow_patterns=["pytorch_model.bin"])
    return local_dir


@torch.no_grad()
def verify(name: str, local_dir: Path) -> None:
    tok = AutoTokenizer.from_pretrained(local_dir, local_files_only=True)
    model = EsmForMaskedLM.from_pretrained(local_dir, local_files_only=True).eval()

    n_params = sum(p.numel() for p in model.parameters())
    cfg = model.config
    print(
        f"[{name}] loaded: {n_params / 1e6:.1f}M params, {cfg.num_hidden_layers} layers, "
        f"{cfg.num_attention_heads} heads, d_model={cfg.hidden_size}"
    )

    # Mask one residue and check the model recovers it (index +1 for the BOS/<cls> token)
    pos = 10
    inputs = tok(TEST_SEQ, return_tensors="pt")
    true_id = inputs["input_ids"][0, pos + 1].item()
    inputs["input_ids"][0, pos + 1] = tok.mask_token_id
    probs = model(**inputs).logits[0, pos + 1].softmax(-1)
    top_id = probs.argmax().item()
    print(
        f"[{name}] smoke test at position {pos}: true={tok.convert_ids_to_tokens(true_id)}, "
        f"predicted={tok.convert_ids_to_tokens(top_id)} (p={probs[top_id]:.2f}), "
        f"p(true)={probs[true_id]:.2f}"
    )
    print(f"[{name}] on disk: {dir_size_mb(local_dir):.0f} MB")


def main() -> None:
    torch.set_num_threads(max(1, torch.get_num_threads()))
    for name, repo_id in MODELS.items():
        local_dir = download(name, repo_id)
        verify(name, local_dir)
    print("\nDone. Load models with EsmForMaskedLM.from_pretrained('models/<name>', local_files_only=True)")


if __name__ == "__main__":
    main()