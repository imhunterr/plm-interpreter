/* PLM Interpreter demo front end. Plain JS, no build step. */
"use strict";

const AA = "ACDEFGHIKLMNPQRSTVWY";
const MAP_ORDER = "AVILMFWCGPSTYNQHKRDE"; // grouped by chemistry for the mutation maps
const AA_NAME = { A: "Alanine", C: "Cysteine", D: "Aspartate", E: "Glutamate", F: "Phenylalanine", G: "Glycine", H: "Histidine", I: "Isoleucine", K: "Lysine", L: "Leucine", M: "Methionine", N: "Asparagine", P: "Proline", Q: "Glutamine", R: "Arginine", S: "Serine", T: "Threonine", V: "Valine", W: "Tryptophan", Y: "Tyrosine" };
const AA_CLASS = { A: "hydrophobic", V: "hydrophobic", I: "hydrophobic", L: "hydrophobic", M: "hydrophobic", F: "hydrophobic", W: "hydrophobic", C: "hydrophobic", G: "special", P: "special", S: "polar", T: "polar", Y: "polar", N: "polar", Q: "polar", H: "positive", K: "positive", R: "positive", D: "negative", E: "negative" };
const CLASS_LABEL = { hydrophobic: "water-avoiding", polar: "polar", positive: "positive charge", negative: "negative charge", special: "glycine / proline" };
const SS_NAME = { H: "helix", E: "strand", "-": "loop" };

const state = {
  proteins: [], protein: null, tokens: null, model: "35M",
  pos: null, mut: null, predict: null, score: null, explain: null, inside: null, land: null,
  colorMode: "class", showIds: false, req: 0,
};

/* ------------------------------------------------------------ helpers */
const $ = (sel, root = document) => root.querySelector(sel);
function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") node.className = v;
    else if (k === "html") node.innerHTML = v;
    else if (k.startsWith("on")) node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v === true ? "" : v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined && c !== false) node.append(c instanceof Node ? c : document.createTextNode(c));
  return node;
}
function svg(tag, attrs = {}, ...children) {
  const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) node.setAttribute(k, v);
  for (const c of children.flat()) if (c) node.append(c instanceof Node ? c : document.createTextNode(c));
  return node;
}
const fmt = (x, d = 2) => (x > 0 ? "+" : x < 0 ? "−" : "") + Math.abs(x).toFixed(d);
const pct = (x) => `${Math.round(x)}%`;
const label = (pos, wt, mt) => `${wt}${pos + 1}${mt}`;

async function api(path, body) {
  const opts = body ? { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) } : {};
  const res = await fetch(path, opts);
  const data = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(data.error || data.detail?.[0]?.msg || `Request failed (${res.status})`);
  return data;
}

function cssVar(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
function hexToRgb(h) { h = h.replace("#", ""); if (h.length === 3) h = h.split("").map((c) => c + c).join(""); const n = parseInt(h, 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
function mix(a, b, t) { const A = hexToRgb(a), B = hexToRgb(b); return "#" + A.map((v, i) => Math.round(v + (B[i] - v) * t).toString(16).padStart(2, "0")).join(""); }
function diverging(v, vmax) {
  const t = Math.max(-1, Math.min(1, vmax > 0 ? v / vmax : 0));
  const s = Math.pow(Math.abs(t), 0.75);
  return t < 0 ? mix(cssVar("--mid"), cssVar("--down"), s) : mix(cssVar("--mid"), cssVar("--up"), s);
}
function ramp(t) { return mix(cssVar("--ramp-lo"), cssVar("--ramp-hi"), Math.max(0, Math.min(1, t))); }

function status(id, text, kind = "busy") {
  const node = $(id);
  node.className = "status" + (kind === "error" ? " error" : "");
  node.replaceChildren(kind === "busy" && text ? el("span", { class: "spinner", "aria-hidden": "true" }) : "", text || "");
}

/* tooltip */
const tip = $("#tooltip");
function showTip(html, x, y) {
  tip.innerHTML = html; tip.hidden = false;
  const r = tip.getBoundingClientRect();
  tip.style.left = Math.min(window.innerWidth - r.width - 8, x + 14) + "px";
  tip.style.top = Math.min(window.innerHeight - r.height - 8, y + 14) + "px";
}
function hideTip() { tip.hidden = true; }

/* ------------------------------------------------------------ theme, tabs, model */
const THEMES = ["auto", "light", "dark"];
let themeIdx = 0;
try { themeIdx = Math.max(0, THEMES.indexOf(localStorage.getItem("plm-theme") || "auto")); } catch (e) { /* storage unavailable */ }
function applyTheme() {
  const t = THEMES[themeIdx];
  if (t === "auto") delete document.documentElement.dataset.theme; else document.documentElement.dataset.theme = t;
  $("#theme-toggle").textContent = `Theme: ${t}`;
  try { localStorage.setItem("plm-theme", t); } catch (e) { /* ignore */ }
  repaint();
}
$("#theme-toggle").addEventListener("click", () => { themeIdx = (themeIdx + 1) % THEMES.length; applyTheme(); });
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", repaint);

document.querySelectorAll(".tab").forEach((tab) => tab.addEventListener("click", () => {
  document.querySelectorAll(".tab").forEach((t) => t.setAttribute("aria-selected", String(t === tab)));
  document.querySelectorAll(".view").forEach((v) => { v.hidden = v.id !== `view-${tab.dataset.view}`; });
  window.scrollTo({ top: 0 });
  if (tab.dataset.view === "demo" && viewer) setTimeout(() => { viewer.resize(); viewer.render(); }, 30);
}));

document.querySelectorAll("[data-model]").forEach((b) => b.addEventListener("click", async () => {
  if (state.model === b.dataset.model) return;
  state.model = b.dataset.model;
  document.querySelectorAll("[data-model]").forEach((x) => x.setAttribute("aria-pressed", String(x === b)));
  state.land = null; $("#map-out").replaceChildren(); status("#map-status", "");
  if (state.pos !== null) { const mut = state.mut; await selectPosition(state.pos); if (mut) await chooseMutant(mut); }
}));

/* ------------------------------------------------------------ proteins */
async function init() {
  applyTheme();
  initViewer();
  renderFindings();
  try {
    state.proteins = await api("/api/proteins");
  } catch (e) {
    $("#protein-head").replaceChildren(el("p", { class: "status error" }, `Could not reach the model server: ${e.message}`));
    return;
  }
  const picker = $("#picker");
  for (const p of state.proteins) {
    picker.append(el("button", { class: "chip", "data-key": p.key, "aria-pressed": "false", onclick: () => loadProtein(p.key) }, p.name, el("small", {}, `${p.length} aa`)));
  }
  picker.append(el("button", { class: "chip ghost", id: "chip-custom", "aria-pressed": "false", onclick: openCustom }, "Your own sequence"));
  await loadProtein(state.proteins[0].key);
}

async function loadProtein(key) {
  document.querySelectorAll(".chip").forEach((c) => c.setAttribute("aria-pressed", String(c.dataset.key === key)));
  $("#protein-head").replaceChildren(el("p", { class: "status" }, el("span", { class: "spinner" }), "Loading protein…"));
  const d = await api(`/api/protein/${key}`);
  await setProtein(d);
}

async function setProtein(d) {
  state.protein = d;
  state.tokens = await api("/api/tokens", { sequence: d.sequence });
  resetSelection();
  state.land = null; $("#map-out").replaceChildren(); status("#map-status", "");
  renderHead();
  renderSequence();
  renderSuggestions();
  showStructure();
  renderSiteInfo();
}

function renderHead() {
  const d = state.protein;
  const meta = [el("span", {}, el("b", {}, d.sequence.length), " amino acids")];
  if (d.gene) meta.push(el("span", {}, "Gene ", el("b", {}, d.gene)));
  if (d.assay) meta.push(el("span", {}, "Lab data: ", el("b", {}, `${d.assay}`), ` (${d.assay_type.toLowerCase()}, ${Object.keys(d.lab).length.toLocaleString()} mutations)`));
  $("#protein-head").replaceChildren(el("h1", {}, d.name), el("p", {}, d.blurb || ""), el("div", { class: "meta" }, meta));
}

function openCustom() { $("#custom-error").textContent = ""; $("#custom-dialog").showModal(); }
$("#custom-cancel").addEventListener("click", () => $("#custom-dialog").close());
$("#custom-form").addEventListener("submit", async (ev) => {
  ev.preventDefault();
  try {
    const toks = await api("/api/tokens", { sequence: $("#custom-seq").value });
    const seq = toks.slice(1, -1).map((t) => t.token).join("");
    $("#custom-dialog").close();
    document.querySelectorAll(".chip").forEach((c) => c.setAttribute("aria-pressed", String(c.id === "chip-custom")));
    await setProtein({ key: null, name: "Your sequence", gene: null, blurb: "A sequence you pasted. The model works the same way; there is no structure or lab data to compare with.", sequence: seq, residues: null, pdb: null, lab: {}, suggest: [] });
  } catch (e) { $("#custom-error").textContent = e.message; }
});

/* ------------------------------------------------------------ step 1: sequence */
function renderSequence() {
  const d = state.protein;
  const seqEl = $("#seq");
  seqEl.className = "seq" + (state.showIds ? " show-ids" : "") + (state.colorMode === "importance" && state.explain ? " importance" : "");
  const blocks = [];
  if (state.showIds) blocks.push(el("span", { class: "special-tok" }, `<cls> ${state.tokens[0].id}`));
  for (let start = 0; start < d.sequence.length; start += 10) {
    const tiles = [];
    for (let i = start; i < Math.min(d.sequence.length, start + 10); i++) {
      const aa = d.sequence[i];
      const t = el("button", {
        class: "tile", type: "button", "data-pos": i, "aria-label": `${AA_NAME[aa]} ${i + 1}`,
        style: `--cls: var(--cls-${AA_CLASS[aa]})`,
        onclick: () => selectPosition(i),
        onmouseenter: (e) => showTip(tileTip(i), e.clientX, e.clientY), onmouseleave: hideTip,
      }, el("span", {}, aa), el("span", { class: "tid" }, state.tokens[i + 1].id));
      tiles.push(t);
    }
    blocks.push(el("div", { class: "block" }, el("span", { class: "block-num" }, start + 1), el("div", { class: "block-tiles" }, tiles)));
  }
  if (state.showIds) blocks.push(el("span", { class: "special-tok" }, `<eos> ${state.tokens[state.tokens.length - 1].id}`));
  seqEl.replaceChildren(...blocks);
  paintTiles();
  const legend = $("#seq-legend");
  if (state.colorMode === "importance" && state.explain) {
    legend.replaceChildren(
      el("span", {}, el("i", { class: "sw", style: `background:${cssVar("--down")}` }), "made the mutation look worse"),
      el("span", {}, el("i", { class: "sw", style: `background:${cssVar("--up")}` }), "made it look more acceptable"),
      el("button", { class: "icon-btn", type: "button", onclick: () => { state.colorMode = "class"; renderSequence(); styleStructure(); } }, "Show amino-acid colours"),
    );
  } else {
    legend.replaceChildren(...Object.entries(CLASS_LABEL).map(([k, v]) => el("span", {}, el("i", { style: `background: var(--cls-${k})` }), v)));
    if (state.explain) legend.append(el("button", { class: "icon-btn", type: "button", onclick: () => { state.colorMode = "importance"; renderSequence(); styleStructure(); } }, "Show importance colours"));
  }
}
$("#show-ids").addEventListener("change", (e) => { state.showIds = e.target.checked; renderSequence(); });

function paintTiles() {
  const imp = state.colorMode === "importance" && state.explain && state.explain.position === state.pos ? state.explain.attribution : null;
  const vmax = imp ? Math.max(...imp.map(Math.abs)) : 1;
  document.querySelectorAll("#seq .tile").forEach((t) => {
    const i = +t.dataset.pos;
    t.classList.toggle("selected", i === state.pos);
    t.style.background = imp && i !== state.pos ? diverging(imp[i], vmax) : "";
    t.firstChild.textContent = i === state.pos ? "?" : state.protein.sequence[i];
    t.lastChild.textContent = i === state.pos ? 32 : state.tokens[i + 1].id;
  });
}

function tileTip(i) {
  const d = state.protein, aa = d.sequence[i];
  let html = `<b>${aa}${i + 1}</b> · ${AA_NAME[aa]} (${CLASS_LABEL[AA_CLASS[aa]]})`;
  const r = d.residues?.[i];
  if (r) html += `<br>${SS_NAME[r.ss] || "loop"}, ${r.rsa < 0.25 ? "buried" : "exposed"}${r.conservation !== null ? `, conservation ${r.conservation.toFixed(2)}` : ""}`;
  if (state.explain && state.colorMode === "importance" && i !== state.pos) html += `<br>effect on score: ${fmt(state.explain.attribution[i], 3)}`;
  if (state.showIds) html += `<br>token number ${state.tokens[i + 1].id}`;
  return html;
}

function renderSuggestions() {
  const d = state.protein;
  const box = $("#suggest");
  if (!d.suggest?.length) { box.replaceChildren(); return; }
  box.replaceChildren(el("span", { class: "muted" }, "Try:"), ...d.suggest.map((s) => {
    const wt = d.sequence[s.position - 1];
    return el("button", { type: "button", onclick: async () => { await selectPosition(s.position - 1); await chooseMutant(s.mutant); } }, `${wt}${s.position}${s.mutant}`);
  }));
}

/* ------------------------------------------------------------ steps 2-3: mask and predict */
function resetSelection() {
  state.pos = null; state.mut = null; state.predict = null; state.score = null; state.explain = null; state.inside = null; state.colorMode = "class";
  for (const id of ["#step-mask", "#step-predict", "#step-score", "#step-why", "#step-inside"]) $(id).classList.add("locked");
  for (const id of ["#tokstrip", "#predict-out", "#score-out", "#why-out", "#inside-out"]) $(id).replaceChildren();
  status("#why-status", ""); status("#inside-status", "");
}

async function selectPosition(pos) {
  const d = state.protein;
  const req = ++state.req;
  state.pos = pos; state.mut = null; state.predict = null; state.score = null; state.explain = null; state.inside = null;
  state.colorMode = "class";
  renderSequence();
  for (const id of ["#step-score", "#step-why", "#step-inside"]) $(id).classList.add("locked");
  for (const id of ["#score-out", "#why-out", "#inside-out"]) $(id).replaceChildren();
  status("#why-status", ""); status("#inside-status", "");
  $("#step-mask").classList.remove("locked");
  $("#step-predict").classList.remove("locked");
  renderTokStrip();
  styleStructure(); focusStructure();
  renderSiteInfo();
  drawMaps();
  $("#predict-out").replaceChildren(el("p", { class: "status" }, el("span", { class: "spinner" }), `ESM-2 ${state.model} is guessing the hidden residue…`));
  try {
    const p = await api("/api/predict", { sequence: d.sequence, position: pos, model: state.model });
    if (req !== state.req) return;
    state.predict = p;
    renderPredict();
  } catch (e) {
    $("#predict-out").replaceChildren(el("p", { class: "status error" }, e.message));
  }
}

function renderTokStrip() {
  const d = state.protein, pos = state.pos;
  const lo = Math.max(0, pos - 5), hi = Math.min(d.sequence.length - 1, pos + 5);
  const items = [];
  if (lo > 0) items.push(el("div", { class: "tok" }, el("b", {}, "<cls>"), el("span", {}, state.tokens[0].id)), el("div", { class: "tok gap" }, "…"));
  else items.push(el("div", { class: "tok" }, el("b", {}, "<cls>"), el("span", {}, state.tokens[0].id)));
  for (let i = lo; i <= hi; i++) {
    if (i === pos) items.push(el("div", { class: "tok mask" }, el("b", {}, "<mask>"), el("span", {}, "32")));
    else items.push(el("div", { class: "tok" }, el("b", {}, d.sequence[i]), el("span", {}, state.tokens[i + 1].id)));
  }
  if (hi < d.sequence.length - 1) items.push(el("div", { class: "tok gap" }, "…"));
  items.push(el("div", { class: "tok" }, el("b", {}, "<eos>"), el("span", {}, state.tokens[state.tokens.length - 1].id)));
  $("#tokstrip").replaceChildren(...items);
  const note = $("#step-mask .arrow-note") || el("p", { class: "arrow-note" });
  note.innerHTML = `Position <b>${pos + 1}</b> was <b>${d.sequence[pos]}</b> (${AA_NAME[d.sequence[pos]]}). The model only sees token numbers; <code>${d.sequence[pos]}</code> = ${state.tokens[pos + 1].id} became <code>&lt;mask&gt;</code> = 32.`;
  $("#step-mask .step-body").append(note);
}

function labFor(pos, aa) {
  const d = state.protein;
  return d.lab?.[label(pos, d.sequence[pos], aa)] || null;
}
function labVerdict(p) { return p < 20 ? ["damaging", "bad"] : p < 50 ? ["reduced", "mid"] : ["works", "good"]; }

function renderPredict() {
  const p = state.predict, d = state.protein, pos = state.pos;
  const maxP = p.probs[0].p;
  const hasLab = Object.keys(d.lab || {}).length > 0;
  const rows = [el("span", { class: "h" }, "amino acid"), el("span", { class: "h" }, "model probability"), el("span", { class: "h", style: "text-align:right" }, ""), el("span", { class: "h" }, hasLab ? "lab result" : "")];
  for (const r of p.probs) {
    const isWt = r.aa === p.wt, isMut = r.aa === state.mut;
    const lab = isWt ? null : labFor(pos, r.aa);
    let labCell = el("span", { class: "lab-cell muted" }, isWt ? "" : hasLab ? "not tested" : "");
    if (lab) {
      const [txt] = labVerdict(lab[1]);
      labCell = el("span", { class: "lab-cell" }, el("span", { class: "lab-dot", style: `background:${ramp(1 - lab[1] / 100)}` }), txt);
    }
    rows.push(el("div", {
      class: "prob-row" + (isWt ? " is-wt" : "") + (isMut ? " is-mut" : ""), role: "button", tabindex: isWt ? -1 : 0,
      onclick: () => { if (!isWt) chooseMutant(r.aa); },
      onkeydown: (e) => { if ((e.key === "Enter" || e.key === " ") && !isWt) { e.preventDefault(); chooseMutant(r.aa); } },
      onmousemove: (e) => showTip(`<b>${r.aa}</b> ${AA_NAME[r.aa]}<br>model probability ${(r.p * 100).toFixed(2)}%${lab ? `<br>lab: score ${lab[0].toFixed(3)}, better than ${pct(lab[1])} of tested mutations` : ""}${isWt ? "<br>this is the original residue" : "<br>click to score this mutation"}`, e.clientX, e.clientY),
      onmouseleave: hideTip,
    },
    el("span", { class: "aa-label" }, r.aa, el("small", {}, AA_NAME[r.aa].slice(0, 3))),
    el("div", { class: "bar-track" }, el("div", { class: "bar", style: `width:${Math.max(0.3, (r.p / maxP) * 100)}%` })),
    el("span", { class: "pval" }, r.p >= 0.001 ? `${(r.p * 100).toFixed(1)}%` : "<0.1%"),
    labCell));
  }
  const wtName = AA_NAME[p.wt];
  const verdict = p.wt_rank === 1
    ? `The model's best guess is <b>${p.top}</b>, the original residue (${(p.wt_prob * 100).toFixed(0)}% probability). It has recovered the hidden letter from context alone.`
    : `The model's best guess is <b>${p.top}</b> (${AA_NAME[p.top]}). The original, <b>${p.wt}</b> (${wtName}), is its #${p.wt_rank} choice at ${(p.wt_prob * 100).toFixed(1)}%.`;
  const certainty = p.uncertainty < 0.35 ? "The model is confident here: only a few amino acids fit." : p.uncertainty < 0.6 ? "Several amino acids could fit here." : "The model is unsure: many amino acids could fit here.";
  $("#predict-out").replaceChildren(
    el("div", { class: "callout", html: `${verdict} ${certainty}` }),
    el("div", { class: "prob-table" }, rows),
  );
}

/* ------------------------------------------------------------ step 4: score */
async function chooseMutant(aa) {
  const d = state.protein, pos = state.pos;
  if (pos === null || aa === d.sequence[pos]) return;
  const req = ++state.req;
  state.mut = aa; state.score = null; state.explain = null; state.inside = null;
  if (state.colorMode === "importance") { state.colorMode = "class"; renderSequence(); styleStructure(); }
  if (state.predict) renderPredict();
  for (const id of ["#step-score", "#step-why", "#step-inside"]) $(id).classList.remove("locked");
  $("#why-out").replaceChildren(); $("#inside-out").replaceChildren(); status("#why-status", ""); status("#inside-status", "");
  $("#btn-why").textContent = `Explain ${label(pos, d.sequence[pos], aa)}`;
  $("#btn-inside").textContent = `Look inside for ${label(pos, d.sequence[pos], aa)}`;
  $("#score-out").replaceChildren(el("p", { class: "status" }, el("span", { class: "spinner" }), state.land ? "Scoring…" : "Scoring (the first time, the model scores every possible mutation in this protein for comparison)…"));
  try {
    const s = await api("/api/score", { sequence: d.sequence, position: pos, mutant: aa, model: state.model });
    if (req !== state.req) return;
    state.score = s;
    renderScore();
    drawMaps();
  } catch (e) { $("#score-out").replaceChildren(el("p", { class: "status error" }, e.message)); }
}

function modelVerdict(p) { return p < 15 ? ["likely damaging", "bad"] : p < 40 ? ["unusual", "mid"] : ["looks tolerated", "good"]; }

function renderScore() {
  const s = state.score, d = state.protein;
  const name = label(s.position, s.wt, s.mutant);
  const ratio = Math.exp(Math.abs(s.llr));
  const ratioText = s.llr < 0
    ? `The model finds <b>${s.mutant}</b> about <b>${ratio < 10 ? ratio.toFixed(1) : Math.round(ratio)}×</b> less likely than <b>${s.wt}</b> at this position.`
    : `The model finds <b>${s.mutant}</b> about <b>${ratio.toFixed(1)}×</b> more likely than <b>${s.wt}</b> here: it looks at least as natural as the original.`;
  const [mv, mcls] = modelVerdict(s.percentile);
  const modelBox = el("div", { class: "score-box" },
    el("span", { class: "k" }, `ESM-2 ${s.model} score for ${name}`),
    el("span", { class: "big" }, fmt(s.llr, 2)),
    el("span", { class: `verdict ${mcls}` }, mv),
    el("p", { class: "small", html: ratioText }),
    el("div", { class: "gauge" }, el("span", { class: "marker", style: `left:${s.percentile}%` })),
    el("div", { class: "gauge-labels" }, el("span", {}, "most unusual"), el("span", {}, "most natural")),
    el("p", { class: "small muted" }, `${pct(s.percentile)} of all ${(d.sequence.length * 19).toLocaleString()} possible mutations in this protein look even more unusual.`));
  const lab = labFor(s.position, s.mutant);
  let labBox, agree = null;
  if (lab) {
    const [lv, lcls] = labVerdict(lab[1]);
    labBox = el("div", { class: "score-box" },
      el("span", { class: "k" }, `Lab measurement (${d.assay})`),
      el("span", { class: "big" }, lab[0].toFixed(2)),
      el("span", { class: `verdict ${lcls}` }, lv),
      el("p", { class: "small" }, `In the experiment, ${name} did better than ${pct(lab[1])} of the ${Object.keys(d.lab).length.toLocaleString()} tested mutations.`),
      el("div", { class: "gauge" }, el("span", { class: "marker", style: `left:${lab[1]}%` })),
      el("div", { class: "gauge-labels" }, el("span", {}, "most damaging"), el("span", {}, "works best")));
    const modelBad = s.percentile < 40, labBad = lab[1] < 50;
    if (modelBad === labBad) agree = el("div", { class: "agree yes", html: `<b>The model and the lab agree.</b> ${labBad ? "Both rate this mutation as harmful." : "Both rate this mutation as acceptable."}` });
    else if (modelBad) agree = el("div", { class: "agree no", html: "<b>The model is harsher than the lab.</b> ESM-2 scores how natural a change looks to evolution, not whether the protein still works. Changes that are rare in nature get flagged even when the protein copes, or even works better. That is the case for activating cancer mutations such as HRas G12V." });
    else agree = el("div", { class: "agree no", html: "<b>The lab finds this damaging, but the model barely objects.</b> Some harmful changes still look natural to evolution, for example when the original residue has a job specific to this protein." });
  } else {
    labBox = el("div", { class: "score-box" }, el("span", { class: "k" }, "Lab measurement"), el("p", { class: "muted" }, Object.keys(d.lab || {}).length ? "This mutation was not tested in the experiment." : "No lab data for this sequence."));
  }
  $("#score-out").replaceChildren(el("div", { class: "score-grid" }, modelBox, labBox), agree || "");
}

/* ------------------------------------------------------------ step 5: why */
$("#btn-why").addEventListener("click", async () => {
  const d = state.protein, pos = state.pos, mut = state.mut;
  if (pos === null || !mut) return;
  const req = state.req;
  $("#btn-why").disabled = true;
  status("#why-status", `Computing integrated gradients with ESM-2 ${state.model} (a few seconds)…`);
  try {
    const x = await api("/api/explain", { sequence: d.sequence, position: pos, mutant: mut, model: state.model, protein: d.key });
    if (req !== state.req) return;
    state.explain = x; state.colorMode = "importance";
    status("#why-status", "");
    renderWhy(); renderSequence(); styleStructure(); renderSiteInfo();
  } catch (e) { status("#why-status", e.message, "error"); }
  finally { $("#btn-why").disabled = false; }
});

function renderWhy() {
  const x = state.explain, d = state.protein;
  const vmax = Math.max(...x.attribution.map(Math.abs));
  const rows = x.top.map((r, k) => el("tr", { onclick: () => highlightResidue(r.position), onmouseenter: () => highlightResidue(r.position), onmouseleave: () => highlightResidue(null) },
    el("td", {}, k + 1),
    el("td", { class: "mono" }, `${r.aa}${r.position + 1}`),
    el("td", {}, el("span", { class: "effect-chip", style: `background:${diverging(r.value, vmax)}` }), fmt(r.value, 3)),
    el("td", {}, r.seq_separation),
    x.has_structure ? el("td", {}, `${r.distance.toFixed(1)} Å`) : null,
    x.has_structure ? el("td", {}, r.contact ? el("span", { class: "yes-pill" }, "yes") : "no") : null));
  const head = el("tr", {}, el("th", {}, "#"), el("th", {}, "residue"), el("th", {}, "effect"), el("th", {}, "positions away"), x.has_structure ? el("th", {}, "3D distance") : null, x.has_structure ? el("th", {}, "touching?") : null);
  const nContact = x.top.filter((r) => r.contact).length;
  const nLocal = x.top.filter((r) => r.seq_separation <= 4).length;
  const nFar = x.top.filter((r) => r.contact && r.seq_separation >= 6).length;
  let summary = `<b>${nLocal} of the top 10</b> residues sit within 4 positions of the site in the sequence.`;
  if (x.has_structure) summary += ` <b>${nContact} of 10</b> physically touch it in 3D (within 8 Å)${nFar ? `, including <b>${nFar}</b> that are far apart in the sequence: the model has picked up the fold` : ""}.`;
  const down = x.top.filter((r) => r.value < 0).length;
  summary += ` ${down} of them push the score down, making ${label(x.position, x.wt, x.mutant)} look worse.`;
  const note = x.converged ? "" : el("p", { class: "small muted" }, `Approximate: the importance values add up to the score change within ${Math.abs(x.completeness_gap).toFixed(2)} (the integral did not fully converge in ${x.steps} steps).`);
  $("#why-out").replaceChildren(
    el("div", { class: "callout", html: summary }),
    el("table", { class: "top-table" }, el("thead", {}, head), el("tbody", {}, rows)),
    note,
  );
}

/* ------------------------------------------------------------ step 6: inside */
$("#btn-inside").addEventListener("click", async () => {
  const d = state.protein, pos = state.pos, mut = state.mut;
  if (pos === null || !mut) return;
  const req = state.req;
  $("#btn-inside").disabled = true;
  status("#inside-status", `Reading every layer and switching off ${state.model === "35M" ? 240 : 120} heads one at a time…`);
  try {
    const x = await api("/api/inside", { sequence: d.sequence, position: pos, mutant: mut, model: state.model });
    if (req !== state.req) return;
    state.inside = x;
    status("#inside-status", "");
    renderInside();
  } catch (e) { status("#inside-status", e.message, "error"); }
  finally { $("#btn-inside").disabled = false; }
});

function renderInside() {
  const x = state.inside;
  const final = x.layers[x.layers.length - 1];
  const settle = x.layers.findIndex((v) => Math.abs(v) >= 0.5 * Math.abs(final) && Math.sign(v) === Math.sign(final));
  const contrib = x.heads.map((row) => row.map((v) => -v)); // head's contribution = clean - ablated
  const flat = contrib.flatMap((row, l) => row.map((v, h) => ({ l, h, v })));
  flat.sort((a, b) => Math.abs(b.v) - Math.abs(a.v));
  const top3 = flat.slice(0, 3).map((t) => `L${t.l}H${t.h} (${fmt(t.v, 2)})`).join(", ");
  const known = x.known_heads;
  const contactSet = new Set(known.contact.map(([l, h]) => `${l}-${h}`));
  const topContact = flat.slice(0, 10).filter((t) => contactSet.has(`${t.l}-${t.h}`)).length;
  const summary = `The score is ${fmt(final, 2)} at the end. It ${settle > 0 ? `only takes shape around <b>layer ${settle}</b> of ${x.n_layers}` : "is present from the start"}. The heads that matter most for this mutation are <b>${top3}</b>; ${topContact ? `${topContact} of the 10 most important are 3D-contact heads.` : "none of the 10 most important are 3D-contact heads."}`;
  $("#inside-out").replaceChildren(
    el("div", { class: "callout", html: summary }),
    el("div", { class: "two-col" },
      el("div", { class: "chart-box" }, el("h3", {}, "Score after each layer"), el("p", { class: "small muted" }, "The model's own output layer applied early (logit lens). 0 = before any layer."), layerChart(x)),
      el("div", { class: "chart-box" }, el("h3", {}, "What each attention head contributes"), headLegend(), headChart(x, contrib))),
  );
}

function layerChart(x) {
  const W = 420, H = 230, m = { l: 40, r: 14, t: 18, b: 34 };
  const n = x.layers.length;
  const lo = Math.min(0, ...x.layers), hi = Math.max(0, ...x.layers);
  const pad = (hi - lo) * 0.1 || 1;
  const y0 = lo - pad, y1 = hi + pad;
  const sx = (i) => m.l + (i / (n - 1)) * (W - m.l - m.r);
  const sy = (v) => m.t + (1 - (v - y0) / (y1 - y0)) * (H - m.t - m.b);
  const g = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: "100%", role: "img", "aria-label": "Score after each layer" });
  const ticks = niceTicks(y0, y1, 4);
  for (const t of ticks) {
    g.append(svg("line", { x1: m.l, x2: W - m.r, y1: sy(t), y2: sy(t), class: "axis", "stroke-width": 0.6 }));
    g.append(svg("text", { x: m.l - 6, y: sy(t) + 4, "text-anchor": "end" }, t.toFixed(Math.abs(t) < 1 && t !== 0 ? 1 : 0)));
  }
  g.append(svg("line", { x1: m.l, x2: W - m.r, y1: sy(0), y2: sy(0), class: "zero" }));
  for (let i = 0; i < n; i++) g.append(svg("text", { x: sx(i), y: H - 14, "text-anchor": "middle" }, i));
  g.append(svg("text", { x: W - m.r, y: H - 1, "text-anchor": "end" }, "layer →"));
  g.append(svg("text", { x: 4, y: m.t - 2, "text-anchor": "start" }, "score"));
  const path = x.layers.map((v, i) => `${i ? "L" : "M"}${sx(i).toFixed(1)},${sy(v).toFixed(1)}`).join("");
  g.append(svg("path", { d: path, fill: "none", stroke: cssVar("--accent"), "stroke-width": 2 }));
  x.layers.forEach((v, i) => {
    const last = i === n - 1;
    const c = svg("circle", { cx: sx(i), cy: sy(v), r: last ? 5 : 3.5, fill: last ? cssVar("--mask") : cssVar("--accent"), stroke: cssVar("--surface"), "stroke-width": 1.5 });
    const hit = svg("rect", { x: sx(i) - 12, y: m.t, width: 24, height: H - m.t - m.b, fill: "transparent" });
    hit.addEventListener("mousemove", (e) => showTip(`<b>${i === 0 ? "Before layer 1 (embeddings)" : `After layer ${i}`}</b><br>score ${fmt(v, 2)}<br>model's top guess here: ${x.layer_top_guess[i]}`, e.clientX, e.clientY));
    hit.addEventListener("mouseleave", hideTip);
    g.append(c, hit);
  });
  g.append(svg("text", { x: sx(n - 1) - 6, y: sy(x.layers[n - 1]) + (x.layers[n - 1] < 0 ? 18 : -10), "text-anchor": "end", style: `fill:${cssVar("--fg")}` }, `final ${fmt(x.layers[n - 1], 2)}`));
  return g;
}

function niceTicks(a, b, n) {
  const span = b - a, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const step = [1, 2, 5, 10].map((s) => s * mag).find((s) => span / s <= n + 1);
  const out = [];
  for (let t = Math.ceil(a / step) * step; t <= b + 1e-9; t += step) out.push(Math.abs(t) < 1e-9 ? 0 : t);
  return out;
}

function headLegend() {
  return el("div", { class: "legend" },
    el("span", {}, el("i", { class: "sw", style: `background:${cssVar("--down")}` }), "pushes the score down"),
    el("span", {}, el("i", { class: "sw", style: `background:${cssVar("--up")}` }), "pushes it up"),
    el("span", {}, el("i", { class: "sw", style: `background:transparent; outline:2px solid ${cssVar("--contact")}` }), "3D-contact head"),
    el("span", {}, el("i", { class: "sw", style: `background:${cssVar("--fg")}; width:6px; height:6px; border-radius:50%` }), "matters most across proteins"));
}

function headChart(x, contrib) {
  const cell = x.n_heads >= 20 ? 19 : 22, m = { l: 34, t: 8, b: 26 };
  const W = m.l + x.n_heads * cell + 4, H = m.t + x.n_layers * cell + m.b;
  const vmax = Math.max(0.05, ...contrib.flat().map(Math.abs));
  const g = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: "100%", role: "img", "aria-label": "Effect of switching off each attention head" });
  const contact = new Set(x.known_heads.contact.map(([l, h]) => `${l}-${h}`));
  const important = new Set(x.known_heads.mutation.map(([l, h]) => `${l}-${h}`));
  for (let l = 0; l < x.n_layers; l++) {
    g.append(svg("text", { x: m.l - 6, y: m.t + l * cell + cell * 0.68, "text-anchor": "end" }, `L${l}`));
    for (let h = 0; h < x.n_heads; h++) {
      const v = contrib[l][h];
      const rect = svg("rect", { x: m.l + h * cell + 1, y: m.t + l * cell + 1, width: cell - 2, height: cell - 2, rx: 2, fill: diverging(v, vmax) });
      if (contact.has(`${l}-${h}`)) { rect.setAttribute("stroke", cssVar("--contact")); rect.setAttribute("stroke-width", 2); }
      rect.addEventListener("mousemove", (e) => showTip(`<b>Layer ${l}, head ${h}</b><br>contribution to the score: ${fmt(v, 3)}<br>(switching it off moves the score to ${fmt(state.inside.layers[state.inside.layers.length - 1] - v, 2)})${contact.has(`${l}-${h}`) ? "<br>one of the 3D-contact heads" : ""}${important.has(`${l}-${h}`) ? "<br>one of the heads that matter most across proteins" : ""}`, e.clientX, e.clientY));
      rect.addEventListener("mouseleave", hideTip);
      g.append(rect);
      if (important.has(`${l}-${h}`)) g.append(svg("circle", { cx: m.l + h * cell + cell / 2, cy: m.t + l * cell + cell / 2, r: 2.6, fill: cssVar("--fg"), "pointer-events": "none" }));
    }
  }
  for (let h = 0; h < x.n_heads; h += 5) g.append(svg("text", { x: m.l + h * cell + cell / 2, y: H - 10, "text-anchor": "middle" }, `H${h}`));
  return g;
}

/* ------------------------------------------------------------ step 7: maps */
$("#btn-map").addEventListener("click", async () => {
  const d = state.protein;
  $("#btn-map").disabled = true;
  status("#map-status", `ESM-2 ${state.model} is scoring all ${(d.sequence.length * 19).toLocaleString()} possible mutations…`);
  try {
    const land = await api("/api/landscape", { sequence: d.sequence, model: state.model });
    if (state.protein !== d) return;
    state.land = land;
    status("#map-status", "");
    drawMaps();
  } catch (e) { status("#map-status", e.message, "error"); }
  finally { $("#btn-map").disabled = false; }
});

function spearman(a, b) {
  const rank = (v) => { const idx = v.map((x, i) => [x, i]).sort((p, q) => p[0] - q[0]); const r = new Array(v.length); let i = 0; while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++; for (let k = i; k <= j; k++) r[idx[k][1]] = (i + j) / 2; i = j + 1; } return r; };
  const ra = rank(a), rb = rank(b), n = a.length, ma = (n - 1) / 2;
  let num = 0, da = 0, db = 0;
  for (let i = 0; i < n; i++) { num += (ra[i] - ma) * (rb[i] - ma); da += (ra[i] - ma) ** 2; db += (rb[i] - ma) ** 2; }
  return num / Math.sqrt(da * db);
}

function drawMaps() {
  if (!state.land) return;
  const d = state.protein, L = d.sequence.length, llr = state.land.llr;
  const cw = L > 300 ? 5 : L > 120 ? 6 : 9, ch = 10;
  // model damage percentile among all non-wild-type cells
  const vals = [];
  for (let i = 0; i < L; i++) for (let a = 0; a < 20; a++) if (AA[a] !== d.sequence[i]) vals.push(llr[i][a]);
  const sorted = [...vals].sort((p, q) => p - q);
  const pctOf = (v) => { let lo = 0, hi = sorted.length; while (lo < hi) { const mid = (lo + hi) >> 1; if (sorted[mid] < v) lo = mid + 1; else hi = mid; } return (lo / sorted.length) * 100; };
  const hasLab = Object.keys(d.lab || {}).length > 0;

  const out = [];
  const makeMap = (title, cellColor) => {
    const canvas = el("canvas", { class: "map", "aria-label": title });
    const dpr = window.devicePixelRatio || 1;
    canvas.width = L * cw * dpr; canvas.height = (20 * ch + 14) * dpr;
    canvas.style.width = `${L * cw}px`; canvas.style.height = `${20 * ch + 14}px`;
    const ctx = canvas.getContext("2d"); ctx.scale(dpr, dpr);
    ctx.fillStyle = cssVar("--muted"); ctx.font = `9px ${cssVar("--f-mono")}`; ctx.textAlign = "center";
    for (let i = 9; i < L; i += 10) ctx.fillText(String(i + 1), i * cw + cw / 2, 9);
    for (let i = 0; i < L; i++) for (let r = 0; r < 20; r++) {
      const aa = MAP_ORDER[r];
      ctx.fillStyle = aa === d.sequence[i] ? cssVar("--fg") : cellColor(i, aa);
      if (aa === d.sequence[i]) { ctx.globalAlpha = 0.85; ctx.fillRect(i * cw + cw / 2 - 1.5, 14 + r * ch + ch / 2 - 1.5, 3, 3); ctx.globalAlpha = 1; }
      else ctx.fillRect(i * cw, 14 + r * ch, cw - (cw > 5 ? 1 : 0), ch - 1);
    }
    if (state.pos !== null) { ctx.strokeStyle = cssVar("--mask"); ctx.lineWidth = 2; ctx.strokeRect(state.pos * cw - 1, 13, cw + 1, 20 * ch + 1); }
    canvas.addEventListener("mousemove", (e) => {
      const rect = canvas.getBoundingClientRect();
      const i = Math.floor((e.clientX - rect.left) / cw), r = Math.floor((e.clientY - rect.top - 14) / ch);
      if (i < 0 || i >= L || r < 0 || r >= 20) { hideTip(); return; }
      const aa = MAP_ORDER[r], wt = d.sequence[i];
      if (aa === wt) { showTip(`<b>${wt}${i + 1}</b> original residue`, e.clientX, e.clientY); return; }
      const v = llr[i][AA.indexOf(aa)], lab = d.lab?.[label(i, wt, aa)];
      showTip(`<b>${label(i, wt, aa)}</b><br>model score ${fmt(v, 2)} (more unusual than ${pct(100 - pctOf(v))} of mutations)${lab ? `<br>lab: better than ${pct(lab[1])} of tested mutations` : hasLab ? "<br>lab: not tested" : ""}<br>click to explore`, e.clientX, e.clientY);
    });
    canvas.addEventListener("mouseleave", hideTip);
    canvas.addEventListener("click", async (e) => {
      const rect = canvas.getBoundingClientRect();
      const i = Math.floor((e.clientX - rect.left) / cw), r = Math.floor((e.clientY - rect.top - 14) / ch);
      if (i < 0 || i >= L || r < 0 || r >= 20 || MAP_ORDER[r] === d.sequence[i]) return;
      hideTip();
      await selectPosition(i); await chooseMutant(MAP_ORDER[r]);
      document.getElementById("step-score").scrollIntoView({ behavior: "smooth", block: "start" });
    });
    const labels = el("div", { class: "map-labels", style: `padding-top:14px` }, ...MAP_ORDER.split("").map((a) => el("span", { style: `height:${ch}px; line-height:${ch}px` }, a)));
    return el("div", {}, el("div", { class: "map-title" }, title), el("div", { class: "map-row" }, labels, canvas));
  };
  out.push(makeMap(`Model: ESM-2 ${state.model}`, (i, aa) => ramp(1 - pctOf(llr[i][AA.indexOf(aa)]) / 100)));
  let agreement = "";
  if (hasLab) {
    out.push(makeMap(`Lab: ${d.assay} (${d.assay_type.toLowerCase()})`, (i, aa) => {
      const lab = d.lab[label(i, d.sequence[i], aa)];
      return lab ? ramp(1 - lab[1] / 100) : cssVar("--surface-2");
    }));
    const m = [], l = [];
    for (const [mutant, [score]] of Object.entries(d.lab)) {
      const i = parseInt(mutant.slice(1, -1), 10) - 1, a = AA.indexOf(mutant.slice(-1));
      if (i >= 0 && i < L && a >= 0) { m.push(llr[i][a]); l.push(score); }
    }
    const rho = spearman(m, l);
    agreement = el("div", { class: "callout", html: `<b>Agreement for this protein:</b> Spearman ρ = ${rho.toFixed(2)} over ${m.length.toLocaleString()} tested mutations. 0 would mean no relation, 1 a perfect match. Across all 217 experiments, ESM-2 35M averages 0.32.` });
  }
  $("#map-out").replaceChildren(
    el("div", { class: "legend" }, "harmless", el("span", { class: "ramp" }), "damaging", el("span", {}, "· dot = original residue · grey = not tested")),
    el("div", { class: "map-scroll" }, out),
    agreement,
  );
}

/* ------------------------------------------------------------ 3D structure */
let viewer = null;
let highlighted = null;
function initViewer() {
  if (!window.$3Dmol) { $("#viewer").replaceChildren(el("div", { class: "viewer-empty" }, "3D viewer could not load.")); return; }
  viewer = $3Dmol.createViewer($("#viewer"), { backgroundColor: cssVar("--viewer-bg") || "#ffffff", antialias: true });
  let mx = 0, my = 0;
  $("#viewer").addEventListener("mousemove", (e) => { mx = e.clientX; my = e.clientY; });
  viewer._hoverPos = () => [mx, my];
  window.addEventListener("resize", () => { if (viewer) { viewer.resize(); viewer.render(); } });
}

function showStructure() {
  if (!viewer) return;
  const d = state.protein;
  viewer.clear();
  $("#viewer-note").textContent = d.pdb ? `AlphaFold2 model of ${d.name}` : "no structure";
  if (!d.pdb) {
    viewer.render();
    $("#viewer-legend").replaceChildren(el("span", {}, "There is no 3D structure for a pasted sequence."));
    return;
  }
  viewer.addModel(d.pdb, "pdb");
  viewer.setHoverable({}, true, (atom) => {
    if (!atom) return;
    const [x, y] = viewer._hoverPos();
    showTip(tileTip(atom.resi - 1), x, y);
  }, () => hideTip());
  viewer.setClickable({}, true, (atom) => { hideTip(); selectPosition(atom.resi - 1); });
  styleStructure();
  viewer.zoomTo();
  viewer.render();
}

function styleStructure() {
  if (!viewer || !state.protein?.pdb) return;
  const base = mix(cssVar("--viewer-bg"), cssVar("--muted"), 0.55);
  const imp = state.colorMode === "importance" && state.explain ? state.explain.attribution : null;
  const vmax = imp ? Math.max(...imp.map(Math.abs)) : 1;
  viewer.setStyle({}, { cartoon: imp ? { colorfunc: (a) => (a.resi - 1 === state.pos ? cssVar("--mask") : diverging(imp[a.resi - 1] || 0, vmax)) } : { color: base } });
  if (imp) {
    const top = state.explain.top.map((r) => r.position + 1);
    viewer.addStyle({ resi: top }, { stick: { colorfunc: (a) => diverging(imp[a.resi - 1], vmax), radius: 0.22 } });
  }
  if (state.pos !== null) viewer.addStyle({ resi: state.pos + 1 }, { stick: { color: cssVar("--mask"), radius: 0.3 }, sphere: { color: cssVar("--mask"), radius: 1.4, opacity: 0.6 } });
  if (highlighted !== null) viewer.addStyle({ resi: highlighted + 1 }, { stick: { color: cssVar("--fg"), radius: 0.32 } });
  viewer.setBackgroundColor(cssVar("--viewer-bg"));
  viewer.render();
  $("#viewer-legend").replaceChildren(...(imp
    ? [el("span", {}, el("i", { class: "sw", style: `background:${cssVar("--mask")}` }), "hidden site"), el("span", {}, el("i", { class: "sw", style: `background:${cssVar("--down")}` }), "made it look worse"), el("span", {}, el("i", { class: "sw", style: `background:${cssVar("--up")}` }), "more acceptable"), el("span", {}, "sticks = top 10 residues")]
    : [el("span", {}, el("i", { class: "sw", style: `background:${cssVar("--mask")}` }), state.pos !== null ? "hidden site" : "click a residue to hide it"), el("span", {}, "drag to rotate · scroll to zoom")]));
}

function focusStructure() {
  if (!viewer || !state.protein?.pdb || state.pos === null) return;
  viewer.zoomTo({ within: { distance: 22, sel: { resi: state.pos + 1 } } }, 600);
}

function highlightResidue(pos) {
  highlighted = pos;
  styleStructure();
  document.querySelectorAll("#seq .tile").forEach((t) => t.classList.toggle("neighbour", +t.dataset.pos === pos));
}

function renderSiteInfo() {
  const d = state.protein, box = $("#site-info");
  if (!d || state.pos === null) { box.replaceChildren(); return; }
  const i = state.pos, aa = d.sequence[i], r = d.residues?.[i];
  const items = [["Hidden residue", `${aa}${i + 1} · ${AA_NAME[aa]}`]];
  if (r) {
    items.push(["Shape", SS_NAME[r.ss] || "loop"]);
    items.push(["Position", `${r.rsa < 0.25 ? "buried inside the protein" : "on the surface"} (exposure ${r.rsa.toFixed(2)})`]);
    if (r.conservation !== null) items.push(["Conservation", `${r.conservation.toFixed(2)} (${r.conservation > 0.6 ? "highly conserved across species" : r.conservation > 0.35 ? "moderately conserved" : "variable across species"})`]);
    items.push(["AlphaFold confidence", `${r.plddt.toFixed(0)} / 100`]);
  }
  box.replaceChildren(...items.flatMap(([k, v]) => [el("dt", {}, k), el("dd", {}, v)]));
}

/* ------------------------------------------------------------ repaint on theme change */
function repaint() {
  if (!state.protein) return;
  renderSequence();
  styleStructure();
  if (state.predict) renderPredict();
  if (state.score) renderScore();
  if (state.explain) renderWhy();
  if (state.inside) renderInside();
  drawMaps();
}

/* ------------------------------------------------------------ findings */
const FINDINGS = [
  { fig: "fig1_baseline", tag: "solid", title: "Our scores reproduce the benchmark exactly", stat: "<b>217 / 217</b> experiments match ProteinGym", text: "Before explaining the model we checked we compute its scores correctly. Mean Spearman with lab data: 0.32 for ESM-2 35M and 0.21 for 8M. The model is best on stability experiments, where structure matters most. On 62,727 clinical variants, 35M separates disease-causing from harmless with AUROC 0.74." },
  { fig: "fig2_probing", tag: "solid", title: "The model learned protein structure on its own", stat: "Helix / strand / loop: <b>0.80</b> vs 0.45 from the letter alone", text: "Simple probes read secondary structure, burial, solvent exposure, contact count and conservation from the model's internal numbers on proteins it never saw. The signal grows with depth and is absent in an untrained copy, so it came from training on sequences only." },
  { fig: "fig3_attention_contacts", tag: "solid", title: "A few late attention heads detect 3D contacts", stat: "Head L11H13: <b>11×</b> chance (0.39 vs 0.035)", text: "Some heads link residues that touch in the folded protein. They sit in the last layers. No head of an untrained model does this." },
  { fig: "fig4_attribution", tag: "care", title: "Explanations point at the site's neighbours", stat: "AUROC <b>0.74</b> for 3D contacts (untrained model 0.46–0.52)", text: "Integrated gradients finds residues that physically surround a mutation, and hiding the top one moves the score 10× more than hiding a random one. But simply taking the nearest residues in the sequence works just as well: the model leans heavily on local context." },
  { fig: "fig5_ablation", tag: "care", title: "The score depends on a few heads", stat: "Top 16 heads: accuracy drops <b>0.16</b>; 16 random heads: 0.03", text: "Heads were ranked on half of the positions and tested on the other half. Switching off the 3D-contact heads barely changes mutation scores: seeing the fold and judging a mutation are done by different parts of the model." },
  { fig: "fig6_logit_lens", tag: "solid", title: "The verdict forms in the last layers", stat: "Near zero until layer 8 of 12, then <b>0.34</b>", text: "Reading the model's guess after every layer shows that the early and middle layers build features, and the last third turns them into a judgement about the mutation." },
  { fig: "fig7_evolution_vs_function", tag: "solid", title: "The model knows what breaks a protein, not what improves it", stat: "Ranking within the worst quarter <b>0.18</b> vs best quarter 0.04", text: "Across 209 experiments the model separates damaging mutations well but barely ranks the beneficial ones (74% of experiments, 88% of stability ones). Evolution shows what is forbidden, not what is better. HRas cancer mutations are the extreme case: they make the protein over-active, yet the model calls them unusual." },
  { fig: "fig8_head_roles", tag: "solid", title: "The heads that judge mutations look nearby", stat: "<b>0</b> of the top 10 mutation heads are contact heads", text: "From the hidden site, the heads that matter for the score put about 40% of their attention on residues 1–4 positions away; contact heads look far across the chain. This explains why sequence neighbours explain the score as well as our attributions do." },
];
function renderFindings() {
  $("#findings").replaceChildren(...FINDINGS.map((f, k) => el("article", { class: "finding" },
    el("span", { class: "num" }, `Finding ${k + 1}`),
    el("h2", {}, f.title),
    el("span", { class: `tag ${f.tag}` }, f.tag === "solid" ? "strong result" : "good, with a catch"),
    el("p", { class: "stat", html: f.stat }),
    el("p", {}, f.text),
    el("div", { class: "plate" }, el("img", { src: `/figures/${f.fig}.png`, alt: f.title, loading: "lazy" })))));
}

init();
