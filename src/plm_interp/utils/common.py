"""Small shared helpers: config, logging, seeding, threads."""

from __future__ import annotations

import logging
import random
import sys
from pathlib import Path

import numpy as np
import torch
import yaml

from plm_interp.paths import CONFIGS_DIR, ROOT


def load_config(name: str = "default") -> dict:
    with open(CONFIGS_DIR / f"{name}.yaml") as f:
        return yaml.safe_load(f)


def seed_everything(seed: int = 0) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def set_threads(n: int | None) -> None:
    if n:
        torch.set_num_threads(n)


def get_logger(name: str, logfile: str | None = None) -> logging.Logger:
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%H:%M:%S")
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    if logfile:
        path = ROOT / "logs" / logfile
        path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(path)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path
