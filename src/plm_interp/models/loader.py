"""Load local ESM-2 checkpoints and expose the pieces the analyses need."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import torch
from transformers import AutoTokenizer, EsmConfig, EsmForMaskedLM

from plm_interp.paths import MODELS_DIR

MODEL_NAMES = {
    "8M": "esm2_t6_8M_UR50D",
    "35M": "esm2_t12_35M_UR50D",
}

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# ESM-2 positional limit: 1026 positions minus <cls> and <eos>, minus 2 for safety as in fair-esm
MAX_RESIDUES = 1022


@dataclass
class ModelBundle:
    """A model, its tokenizer, and cached index helpers."""

    name: str
    model: EsmForMaskedLM
    tokenizer: AutoTokenizer

    @property
    def n_layers(self) -> int:
        return self.model.config.num_hidden_layers

    @property
    def n_heads(self) -> int:
        return self.model.config.num_attention_heads

    @property
    def d_model(self) -> int:
        return self.model.config.hidden_size

    @property
    def d_head(self) -> int:
        return self.d_model // self.n_heads

    @property
    def mask_id(self) -> int:
        return self.tokenizer.mask_token_id

    @cached_property
    def aa_ids(self) -> torch.Tensor:
        """Vocabulary ids of the 20 standard amino acids, in AMINO_ACIDS order."""
        ids = self.tokenizer.convert_tokens_to_ids(list(AMINO_ACIDS))
        return torch.tensor(ids, dtype=torch.long)

    def encode(self, seq: str) -> torch.Tensor:
        """Token ids (1, L+2) with <cls> at 0 and <eos> at L+1."""
        return self.tokenizer(seq, return_tensors="pt")["input_ids"]

    @property
    def layers(self):
        return self.model.esm.encoder.layer

    def unembed(self, resid: torch.Tensor) -> torch.Tensor:
        """Final layer norm + LM head applied to a residual-stream tensor (logit lens)."""
        return self.model.lm_head(self.model.esm.encoder.emb_layer_norm_after(resid))


def load_model(size: str = "35M", attn: str = "sdpa", random_init: bool = False, seed: int = 0) -> ModelBundle:
    """Load ESM-2 from models/.

    attn="eager" is needed whenever attention weights are read; "sdpa" is faster otherwise.
    random_init=True returns the same architecture with freshly initialised weights, used as a
    control: any "biology" a probe or attribution finds there comes from the method, not training.
    """
    path = MODELS_DIR / MODEL_NAMES[size]
    tokenizer = AutoTokenizer.from_pretrained(path, local_files_only=True)
    if random_init:
        torch.manual_seed(seed)
        config = EsmConfig.from_pretrained(path, local_files_only=True)
        config._attn_implementation = attn
        model = EsmForMaskedLM(config)
    else:
        model = EsmForMaskedLM.from_pretrained(path, local_files_only=True, attn_implementation=attn)
    model.eval()
    for p in model.parameters():
        p.requires_grad_(False)
    label = f"{size}-random" if random_init else size
    return ModelBundle(name=label, model=model, tokenizer=tokenizer)
