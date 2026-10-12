"""FastAPI server for the local demo.

Run:
    python -m plm_interp.app              # opens http://127.0.0.1:8000 in the browser
    python -m plm_interp.app --warm       # precompute the showcase examples first (recommended before presenting)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import threading
import time
import webbrowser
from pathlib import Path

import torch
import uvicorn
from fastapi import FastAPI, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from plm_interp.app import engine
from plm_interp.paths import FIGURES_DIR, INTERIM_DIR

STATIC = Path(__file__).parent / "static"
CACHE = INTERIM_DIR / "app_cache" / "results"

app = FastAPI(title="PLM Interpreter")
app.mount("/static", StaticFiles(directory=STATIC), name="static")
if FIGURES_DIR.exists():
    app.mount("/figures", StaticFiles(directory=FIGURES_DIR), name="figures")


class SeqRequest(BaseModel):
    sequence: str
    model: str = "35M"


class SiteRequest(SeqRequest):
    position: int


class VariantRequest(SiteRequest):
    mutant: str
    protein: str | None = None


@app.exception_handler(engine.DemoError)
async def demo_error(_request: Request, exc: engine.DemoError):
    return JSONResponse(status_code=400, content={"error": str(exc)})


def _cached(kind: str, payload: dict, compute):
    """Disk cache for the slower calls, so a repeated demo step answers instantly."""
    CACHE.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha1(json.dumps([kind, payload], sort_keys=True).encode()).hexdigest()
    path = CACHE / f"{kind}_{key}.json"
    if path.exists():
        return json.loads(path.read_text())
    result = compute()
    path.write_text(json.dumps(result))
    return result


@app.get("/")
def index():
    return FileResponse(STATIC / "index.html")


@app.get("/api/proteins")
def proteins():
    return engine.protein_list()


@app.get("/api/protein/{key}")
def protein(key: str):
    return engine.protein_detail(key)


@app.post("/api/tokens")
def tokens(req: SeqRequest):
    return engine.tokens(engine.clean_sequence(req.sequence))


@app.post("/api/predict")
def predict(req: SiteRequest):
    return engine.predict(req.sequence, req.position, req.model)


@app.post("/api/landscape")
def landscape(req: SeqRequest):
    return engine.landscape(req.sequence, req.model)


@app.post("/api/score")
def score(req: VariantRequest):
    return engine.score(req.sequence, req.position, req.mutant, req.model)


@app.post("/api/explain")
def explain(req: VariantRequest):
    seq = engine.clean_sequence(req.sequence)
    payload = dict(seq=seq, pos=req.position, mut=req.mutant, model=req.model, protein=req.protein)
    return _cached("explain", payload, lambda: engine.explain(seq, req.position, req.mutant, req.model, req.protein))


@app.post("/api/inside")
def inside(req: VariantRequest):
    seq = engine.clean_sequence(req.sequence)
    payload = dict(seq=seq, pos=req.position, mut=req.mutant, model=req.model)
    return _cached("inside", payload, lambda: engine.inside(seq, req.position, req.mutant, req.model))


def warm(models: list[str]) -> None:
    """Precompute maps, explanations and model internals for every suggested showcase mutation."""
    for size in models:
        for p in engine.protein_list():
            d = engine.protein_detail(p["key"])
            t = time.time()
            engine.landscape(d["sequence"], size)
            for s in d["suggest"]:
                pos = s["position"] - 1
                req = VariantRequest(sequence=d["sequence"], position=pos, mutant=s["mutant"], model=size, protein=d["key"])
                explain(req)
                inside(req)
            print(f"  warmed {size} {p['name']} in {time.time() - t:.0f}s", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description="PLM Interpreter local demo")
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--threads", type=int, default=None, help="CPU threads for the model (default: all)")
    ap.add_argument("--warm", action="store_true", help="precompute the showcase examples before starting")
    ap.add_argument("--warm-models", nargs="+", default=["35M", "8M"])
    ap.add_argument("--no-browser", action="store_true")
    args = ap.parse_args()
    if args.threads:
        torch.set_num_threads(args.threads)

    print("Loading ESM-2 models...", flush=True)
    engine.bundle("35M")
    engine.bundle("8M")
    if args.warm:
        print("Precomputing showcase examples (one-off; results are cached on disk)...", flush=True)
        warm(args.warm_models)

    url = f"http://127.0.0.1:{args.port}"
    if not args.no_browser:
        threading.Timer(1.5, lambda: webbrowser.open(url)).start()
    print(f"PLM Interpreter is running at {url}  (Ctrl+C to stop)", flush=True)
    uvicorn.run(app, host="127.0.0.1", port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
