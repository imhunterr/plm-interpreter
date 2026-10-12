/* PLM Interpreter: strand-and-decoder demo. Plain JS, no build step. */
"use strict";

const AA = "ACDEFGHIKLMNPQRSTVWY";
const WHEEL_ORDER = "AVILMFWCGPSTYNQHKRDE"; // grouped by chemistry
const AA_NAME = { A: "Alanine", C: "Cysteine", D: "Aspartate", E: "Glutamate", F: "Phenylalanine", G: "Glycine", H: "Histidine", I: "Isoleucine", K: "Lysine", L: "Leucine", M: "Methionine", N: "Asparagine", P: "Proline", Q: "Glutamine", R: "Arginine", S: "Serine", T: "Threonine", V: "Valine", W: "Tryptophan", Y: "Tyrosine" };
const AA_CLASS = { A: "hydrophobic", V: "hydrophobic", I: "hydrophobic", L: "hydrophobic", M: "hydrophobic", F: "hydrophobic", W: "hydrophobic", C: "hydrophobic", G: "special", P: "special", S: "polar", T: "polar", Y: "polar", N: "polar", Q: "polar", H: "positive", K: "positive", R: "positive", D: "negative", E: "negative" };
const CLASS_LABEL = { hydrophobic: "water-avoiding", polar: "polar", positive: "positive", negative: "negative", special: "glycine / proline" };
const SS_NAME = { H: "helix", E: "strand", "-": "loop" };
const STRIDE = 30; // bead width + gap in px

const STAGES = [
  { key: "read", label: "Read" },
  { key: "hide", label: "Hide" },
  { key: "guess", label: "Guess" },
  { key: "score", label: "Score" },
  { key: "why", label: "Why?" },
  { key: "inside", label: "Inside" },
  { key: "map", label: "All mutations" },
];

const state = {
  proteins: [], protein: null, tokens: null, model: "35M", stage: 0,
  pos: null, mut: null, predict: null, score: null, explain: null, inside: null, land: null,
  showIds: false, req: 0, busy: {},
};

/* ================================================================ helpers */
const $ = (s, r = document) => r.querySelector(s);
function el(tag, attrs = {}, ...kids) {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") n.className = v;
    else if (k === "html") n.innerHTML = v;
    else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
    else n.setAttribute(k, v === true ? "" : v);
  }
  for (const c of kids.flat()) if (c !== null && c !== undefined && c !== false) n.append(c instanceof Node ? c : document.createTextNode(c));
  return n;
}
function svg(tag, attrs = {}, ...kids) {
  const n = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) n.setAttribute(k, v);
  for (const c of kids.flat()) if (c !== null && c !== undefined) n.append(c instanceof Node ? c : document.createTextNode(c));
  return n;
}
const fmt = (x, d = 2) => (x > 0 ? "+" : x < 0 ? "−" : "") + Math.abs(x).toFixed(d);
const pct = (x) => `${Math.round(x)}%`;
const mname = (pos, wt, mt) => `${wt}${pos + 1}${mt}`;
const cssVar = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
function rgb(h) { h = h.replace("#", ""); if (h.length === 3) h = h.split("").map((c) => c + c).join(""); const n = parseInt(h, 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255]; }
function mix(a, b, t) { const A = rgb(a), B = rgb(b); return "#" + A.map((v, i) => Math.round(v + (B[i] - v) * t).toString(16).padStart(2, "0")).join(""); }
function diverging(v, vmax) {
  const t = Math.max(-1, Math.min(1, vmax > 0 ? v / vmax : 0)), s = Math.pow(Math.abs(t), 0.7);
  return mix(cssVar("--mid"), t < 0 ? cssVar("--down") : cssVar("--up"), s);
}
const ramp = (t) => mix(cssVar("--ramp-lo"), cssVar("--ramp-hi"), Math.max(0, Math.min(1, t)));

async function api(path, body) {
  const opts = body ? { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) } : {};
  const res = await fetch(path, opts);
  const data = await res.json().catch(() => ({}));
  if (!res.ok) throw new Error(data.error || data.detail?.[0]?.msg || `Request failed (${res.status})`);
  return data;
}

const tip = $("#tooltip");
function showTip(html, x, y) {
  tip.innerHTML = html; tip.hidden = false;
  const r = tip.getBoundingClientRect();
  tip.style.left = Math.min(innerWidth - r.width - 10, x + 16) + "px";
  tip.style.top = Math.min(innerHeight - r.height - 10, y + 16) + "px";
}
const hideTip = () => { tip.hidden = true; };
const busy = (text) => el("div", { class: "busy" }, el("span", { class: "spinner", "aria-hidden": "true" }), text);

/* ================================================================ theme, views, model */
let theme = "dark";
try { theme = localStorage.getItem("plm-theme") || "dark"; } catch (e) { /* storage unavailable */ }
function applyTheme() {
  document.documentElement.dataset.theme = theme;
  try { localStorage.setItem("plm-theme", theme); } catch (e) { /* ignore */ }
  if (state.protein) { renderStage(); updateStrand(); styleStructure(); }
}
$("#theme-toggle").addEventListener("click", () => { theme = theme === "dark" ? "light" : "dark"; applyTheme(); });

document.querySelectorAll(".views button").forEach((b) => b.addEventListener("click", () => {
  document.querySelectorAll(".views button").forEach((x) => x.setAttribute("aria-selected", String(x === b)));
  document.querySelectorAll(".view").forEach((v) => { v.hidden = v.id !== `view-${b.dataset.view}`; });
  window.scrollTo({ top: 0 });
  if (b.dataset.view === "demo") resizeViewer();
}));

document.querySelectorAll("[data-model]").forEach((b) => b.addEventListener("click", async () => {
  if (state.model === b.dataset.model) return;
  state.model = b.dataset.model;
  document.querySelectorAll("[data-model]").forEach((x) => x.setAttribute("aria-pressed", String(x === b)));
  const { pos, mut, stage } = state;
  state.land = null;
  if (pos !== null) { await selectPosition(pos, { stay: true }); if (mut) await chooseMutant(mut, { stay: true }); }
  go(stage);
}));

/* ================================================================ proteins */
const menu = $("#protein-list"), menuBtn = $("#protein-btn");
function toggleMenu(open) { menu.hidden = !open; menuBtn.setAttribute("aria-expanded", String(open)); }
menuBtn.addEventListener("click", (e) => { e.stopPropagation(); toggleMenu(menu.hidden); });
document.addEventListener("click", (e) => { if (!menu.hidden && !menu.contains(e.target)) toggleMenu(false); });

async function init() {
  applyTheme();
  initViewer();
  renderFindings();
  try { state.proteins = await api("/api/proteins"); }
  catch (e) { $("#stage-body").replaceChildren(el("p", { class: "error-text" }, `Could not reach the model server: ${e.message}`)); return; }
  menu.replaceChildren(...state.proteins.map((p) => el("button", { role: "option", "data-key": p.key, onclick: () => { toggleMenu(false); loadProtein(p.key); } },
    el("b", {}, p.name), el("span", { class: "len" }, `${p.length} aa`), el("small", {}, p.blurb))),
    el("button", { role: "option", onclick: () => { toggleMenu(false); openCustom(); } }, el("b", {}, "Your own sequence"), el("span", { class: "len" }, "paste"), el("small", {}, "Paste any protein sequence and run the same analysis.")));
  await loadProtein(state.proteins[0].key);
}

async function loadProtein(key) {
  $("#protein-name").textContent = "Loading…";
  const d = await api(`/api/protein/${key}`);
  await setProtein(d);
}

async function setProtein(d) {
  state.protein = d;
  state.tokens = await api("/api/tokens", { sequence: d.sequence });
  Object.assign(state, { pos: null, mut: null, predict: null, score: null, explain: null, inside: null, land: null });
  state.req++;
  $("#protein-name").textContent = d.name;
  menu.querySelectorAll("button").forEach((b) => b.setAttribute("aria-selected", String(b.dataset.key === d.key)));
  renderStrand();
  showStructure();
  go(0);
}

function openCustom() { $("#custom-error").textContent = ""; $("#custom-dialog").showModal(); }
$("#custom-cancel").addEventListener("click", () => $("#custom-dialog").close());
$("#custom-form").addEventListener("submit", async (ev) => {
  ev.preventDefault();
  try {
    const toks = await api("/api/tokens", { sequence: $("#custom-seq").value });
    const seq = toks.slice(1, -1).map((t) => t.token).join("");
    $("#custom-dialog").close();
    await setProtein({ key: null, name: "Your sequence", gene: null, blurb: "A sequence you pasted.", sequence: seq, residues: null, pdb: null, lab: {}, suggest: [] });
  } catch (e) { $("#custom-error").textContent = e.message; }
});

/* ================================================================ stage machine */
function unlocked(i) {
  if (i === 0 || i === 6) return true;
  if (i <= 2) return state.pos !== null;
  return state.pos !== null && state.mut !== null;
}

function go(i) {
  i = Math.max(0, Math.min(STAGES.length - 1, i));
  if (!unlocked(i)) return;
  state.stage = i;
  renderRail();
  renderStage();
  updateStrand();
  styleStructure();
  const wide = ["inside", "map"].includes(STAGES[i].key);
  $("#workspace").classList.toggle("wide", wide);
  $("#view-demo").classList.toggle("no-strand", wide);
  $("#strand").classList.toggle("compact", STAGES[i].key !== "why");
  resizeViewer();
  if (STAGES[i].key === "why" && !state.explain) runExplain();
  if (STAGES[i].key === "inside" && !state.inside) runInside();
  if (STAGES[i].key === "map" && !state.land) runLandscape();
}

function renderRail() {
  $("#rail").replaceChildren(...STAGES.map((s, i) => el("button", {
    type: "button", disabled: !unlocked(i), "aria-current": i === state.stage ? "step" : null,
    class: i < state.stage && unlocked(i) ? "done" : null, onclick: () => go(i),
  }, el("span", { class: "n" }, i + 1), el("span", { class: "t" }, s.label))));
  const prev = state.stage - 1, next = state.stage + 1;
  $("#prev").disabled = prev < 0;
  $("#prev-label").textContent = prev >= 0 ? STAGES[prev].label : "Back";
  const canNext = next < STAGES.length && unlocked(next);
  $("#next").disabled = !canNext;
  $("#next-label").textContent = next >= STAGES.length ? "Done" : canNext ? STAGES[next].label : state.pos === null ? "Pick a residue first" : "Pick a replacement first";
}
$("#prev").addEventListener("click", () => go(state.stage - 1));
$("#next").addEventListener("click", () => go(state.stage + 1));
document.addEventListener("keydown", (e) => {
  if ($("#view-demo").hidden || e.target.closest("textarea, input, dialog")) return;
  if (e.key === "ArrowRight") { e.preventDefault(); go(state.stage + 1); }
  if (e.key === "ArrowLeft") { e.preventDefault(); go(state.stage - 1); }
  if (e.key === "Escape") toggleMenu(false);
});

function head(title, text) {
  const s = STAGES[state.stage];
  $("#stage-kicker").textContent = `Step ${state.stage + 1} of ${STAGES.length} · ${s.label}`;
  $("#stage-title").textContent = title;
  $("#stage-text").innerHTML = text;
}

function renderStage() {
  const body = $("#stage-body");
  const key = STAGES[state.stage].key;
  const content = { read: stageRead, hide: stageHide, guess: stageGuess, score: stageScore, why: stageWhy, inside: stageInside, map: stageMap }[key]();
  body.replaceChildren(...[content].flat());
  body.style.animation = "none"; void body.offsetWidth; body.style.animation = "";
  renderRail();
}

/* ================================================================ stage 1: read */
function stageRead() {
  const d = state.protein;
  head("A protein is a sentence written in 20 letters",
    `ESM-2 has read millions of natural proteins and learned one skill: <b>hide a letter, then guess it back</b> from everything around it. This is <b>${d.name}</b>${d.gene ? ` (${d.gene})` : ""}, shown as a strand of beads below. ${d.blurb || ""}`);
  const nLab = Object.keys(d.lab || {}).length;
  const facts = el("div", { class: "facts" },
    el("div", { class: "fact" }, el("b", {}, d.sequence.length), el("span", {}, "amino acids")),
    el("div", { class: "fact" }, el("b", {}, nLab ? nLab.toLocaleString() : "—"), el("span", {}, nLab ? `mutations measured (${d.assay})` : "no lab data")),
    el("div", { class: "fact" }, el("b", {}, d.pdb ? "3D" : "—"), el("span", {}, d.pdb ? "AlphaFold2 structure" : "no structure")),
    el("div", { class: "fact" }, el("b", {}, state.model), el("span", {}, `ESM-2 model (${state.model === "35M" ? "12 layers, 240 heads" : "6 layers, 120 heads"})`)));
  const out = [facts, el("p", { class: "lead" }, "Click any bead on the strand to hide that residue, or start from a famous mutation:")];
  if (d.suggest?.length) {
    out.push(el("div", { class: "try" }, ...d.suggest.map((s) => {
      const wt = d.sequence[s.position - 1];
      return el("button", { type: "button", onclick: async () => { await selectPosition(s.position - 1, { stay: true }); await chooseMutant(s.mutant); } }, `${wt}${s.position}${s.mutant}`);
    })));
  }
  out.push(el("div", { class: "callout", html: "Bead colours show each amino acid's chemistry. Tick <b>token numbers</b> to see what the model actually receives: every letter becomes a number." }));
  return out;
}

/* ================================================================ stage 2: hide */
function stageHide() {
  const d = state.protein, pos = state.pos, wt = d.sequence[pos];
  head(`Hide position ${pos + 1}`,
    `The <b>${AA_NAME[wt]}</b> at position ${pos + 1} is swapped for a special <code>&lt;mask&gt;</code> token. The model sees every other residue, and nothing about what was here.`);
  const lo = Math.max(0, pos - 6), hi = Math.min(d.sequence.length - 1, pos + 6);
  const toks = [el("div", { class: "tok" }, el("b", {}, "<cls>"), el("span", {}, state.tokens[0].id))];
  if (lo > 0) toks.push(el("div", { class: "tok gap" }, "…"));
  for (let i = lo; i <= hi; i++) {
    toks.push(i === pos ? el("div", { class: "tok mask" }, el("b", {}, "<mask>"), el("span", {}, "32"))
      : el("div", { class: "tok" }, el("b", {}, d.sequence[i]), el("span", {}, state.tokens[i + 1].id)));
  }
  if (hi < d.sequence.length - 1) toks.push(el("div", { class: "tok gap" }, "…"));
  toks.push(el("div", { class: "tok" }, el("b", {}, "<eos>"), el("span", {}, state.tokens[state.tokens.length - 1].id)));
  return [
    el("div", { class: "tokens" }, toks),
    el("div", { class: "callout", html: `The model never sees letters, only numbers. <code>${wt}</code> was token <b>${state.tokens[pos + 1].id}</b>; now it is <code>&lt;mask&gt;</code>, token <b>32</b>. Two extra tokens mark the start and end of the protein.` }),
    el("p", { class: "muted" }, "Press → to see what the model guesses."),
  ];
}

/* ================================================================ stage 3: guess (decoder wheel) */
function stageGuess() {
  const d = state.protein, pos = state.pos;
  head("The model fills in the blank",
    `ESM-2 gives every amino acid a probability for position ${pos + 1}. Longer petals mean more likely. Dots on the rim show what the lab measured for that replacement. <b>Tap a petal</b> to score that mutation.`);
  if (!state.predict) return busy(`ESM-2 ${state.model} is reading the protein…`);
  const p = state.predict;
  const verdict = p.wt_rank === 1
    ? `Its top guess is <b>${p.top}</b>, exactly the original residue, at ${(p.wt_prob * 100).toFixed(0)}%. It recovered the hidden letter from context alone.`
    : `Its top guess is <b>${p.top}</b> (${AA_NAME[p.top]}). The original <b>${p.wt}</b> is its #${p.wt_rank} choice at ${(p.wt_prob * 100).toFixed(1)}%.`;
  const certainty = p.uncertainty < 0.35 ? "Very few amino acids fit here." : p.uncertainty < 0.6 ? "A handful of amino acids could fit here." : "Many amino acids could fit here: the model is unsure.";
  const top = p.probs.slice(0, 6);
  const list = el("div", { class: "toplist" }, el("span", { class: "kicker" }, "Top guesses"), ...top.map((r) => el("div", { class: "tl-row" + (r.aa === p.wt ? " wt" : "") },
    el("span", { class: "aa" }, r.aa), el("div", { class: "tl-bar" }, el("i", { style: `width:${(r.p / top[0].p) * 100}%` })), el("span", { class: "p" }, `${(r.p * 100).toFixed(1)}%`))));
  const hasLab = Object.keys(d.lab || {}).length > 0;
  const legend = el("div", { class: "legend" },
    el("span", {}, el("i", { style: `background:${cssVar("--primary")}` }), "model probability"),
    el("span", {}, el("i", { style: `background:${cssVar("--faint")}` }), "original residue"),
    hasLab ? el("span", {}, el("i", { style: `background:${ramp(0.9)}` }), "lab: damaging") : null,
    hasLab ? el("span", {}, el("i", { style: `background:${ramp(0.1)}; box-shadow: inset 0 0 0 1px ${cssVar("--line")}` }), "lab: works") : null);
  return el("div", { class: "decoder" }, wheel(p), el("div", { class: "toplist-col", style: "display:flex;flex-direction:column;gap:16px" },
    el("div", { class: "callout", html: `${verdict} ${certainty}` }), list, legend,
    el("p", { class: "muted", style: "font-size:.82rem" }, "Petal length is proportional to the square root of the probability, so small values stay visible.")));
}

function wheel(p) {
  const d = state.protein, pos = state.pos, C = 220, rs = 86, rMax = 166, rLab = 186, rDot = 205;
  const prob = Object.fromEntries(p.probs.map((r) => [r.aa, r.p]));
  const pmax = Math.max(...p.probs.map((r) => r.p));
  const g = svg("svg", { class: "wheel", viewBox: "0 0 440 440", role: "img", "aria-label": "Probability of each amino acid at the hidden position" });
  const ang = (k) => ((-90 + (k + 0.5) * 18) * Math.PI) / 180;
  const pt = (r, a) => [C + r * Math.cos(a), C + r * Math.sin(a)];
  // chemistry groups as arcs around the hub
  const groups = {};
  [...WHEEL_ORDER].forEach((a, k) => { (groups[AA_CLASS[a]] ||= []).push(k); });
  for (const [cls, ks] of Object.entries(groups)) {
    const a0 = ang(ks[0]) - 0.13, a1 = ang(ks[ks.length - 1]) + 0.13, r = 76;
    const [x0, y0] = pt(r, a0), [x1, y1] = pt(r, a1);
    g.append(svg("path", { d: `M${x0},${y0} A${r},${r} 0 ${a1 - a0 > Math.PI ? 1 : 0} 1 ${x1},${y1}`, fill: "none", stroke: cssVar(`--cls-${cls}`), "stroke-width": 4, "stroke-linecap": "round", opacity: 0.9 }));
  }
  for (let k = 0; k < 20; k++) {
    const aa = WHEEL_ORDER[k], a = ang(k), isWt = aa === p.wt, isMut = aa === state.mut;
    const len = Math.max(5, Math.sqrt(prob[aa] / pmax) * (rMax - rs));
    const [x0, y0] = pt(rs, a), [x1, y1] = pt(rs + len, a), [lx, ly] = pt(rLab, a), [hx, hy] = pt(rDot + 8, a);
    const color = isMut ? cssVar("--mask") : isWt ? cssVar("--faint") : cssVar("--primary");
    const grp = svg("g", { class: "petal-g" });
    const [tx, ty] = pt(rMax, a);
    grp.append(svg("line", { x1: x0, y1: y0, x2: tx, y2: ty, stroke: cssVar("--surface-3"), "stroke-width": 13, "stroke-linecap": "round", opacity: 0.6 }));
    grp.append(svg("line", { x1: x0, y1: y0, x2: x1, y2: y1, stroke: color, "stroke-width": 13, "stroke-linecap": "round", class: "petal" }));
    grp.append(svg("text", { x: lx, y: ly + 5, "text-anchor": "middle", "font-size": 15, "font-weight": isWt || isMut ? 700 : 500, style: `fill:${isMut ? cssVar("--mask") : cssVar("--fg")}` }, aa));
    const lab = isWt ? null : d.lab?.[mname(pos, p.wt, aa)];
    if (lab) { const [dx, dy] = pt(rDot, a); grp.append(svg("circle", { cx: dx, cy: dy, r: 5.5, fill: ramp(1 - lab[1] / 100), stroke: cssVar("--line"), "stroke-width": 1 })); }
    if (isWt) { const [dx, dy] = pt(rDot, a); grp.append(svg("text", { x: dx, y: dy + 4, "text-anchor": "middle", "font-size": 10 }, "orig")); }
    const hit = svg("line", { x1: x0, y1: y0, x2: hx, y2: hy, stroke: "transparent", "stroke-width": 26, style: `cursor:${isWt ? "default" : "pointer"}` });
    hit.addEventListener("mousemove", (e) => showTip(`<b>${aa} · ${AA_NAME[aa]}</b><br>probability ${(prob[aa] * 100).toFixed(2)}%${lab ? `<br>lab: did better than ${pct(lab[1])} of tested mutations` : ""}${isWt ? "<br>the original residue" : "<br>click to score " + mname(pos, p.wt, aa)}`, e.clientX, e.clientY));
    hit.addEventListener("mouseleave", hideTip);
    if (!isWt) hit.addEventListener("click", () => { hideTip(); chooseMutant(aa); });
    grp.append(hit);
    g.append(grp);
  }
  // hub with the masked bead
  g.append(svg("circle", { cx: C, cy: C, r: 64, fill: cssVar("--surface-2") }));
  g.append(svg("circle", { cx: C, cy: C - 8, r: 30, fill: cssVar("--mask") }));
  g.append(svg("text", { x: C, y: C + 2, "text-anchor": "middle", "font-size": 30, "font-weight": 700, style: "fill:#fff" }, "?"));
  g.append(svg("text", { x: C, y: C + 44, "text-anchor": "middle", "font-size": 12 }, `position ${pos + 1}`));
  return g;
}

/* ================================================================ stage 4: score */
function stageScore() {
  const d = state.protein, pos = state.pos, wt = d.sequence[pos], mt = state.mut;
  head(`How bad is ${mname(pos, wt, mt)}?`,
    `The score compares the new amino acid with the original: <code>log P(${mt}) − log P(${wt})</code>. Zero means "just as natural". Below zero means less natural, and very negative usually means damaging.`);
  if (!state.score) return busy(state.land ? "Scoring…" : "Scoring every possible mutation in this protein for comparison (first time only)…");
  const s = state.score, ratio = Math.exp(Math.abs(s.llr));
  const ratioText = s.llr < 0
    ? `The model finds <b>${mt}</b> about <b>${ratio < 10 ? ratio.toFixed(1) : Math.round(ratio)}×</b> less likely than <b>${wt}</b> here.`
    : `The model finds <b>${mt}</b> about <b>${ratio.toFixed(1)}×</b> more likely than <b>${wt}</b> here.`;
  const [mv, mc] = s.percentile < 15 ? ["likely damaging", "bad"] : s.percentile < 40 ? ["unusual", "mid"] : ["looks tolerated", "good"];
  const lab = d.lab?.[mname(pos, wt, mt)];
  const hero = el("div", { class: "score-hero" },
    el("div", { style: "display:flex;flex-direction:column;gap:8px" },
      el("span", { class: "mutname", html: `<b>${mname(pos, wt, mt)}</b> · ${AA_NAME[wt]} → ${AA_NAME[mt]}` }),
      el("span", { class: "big" }, fmt(s.llr, 2)),
      el("span", { class: `verdict ${mc}` }, `ESM-2: ${mv}`)),
    el("p", { class: "lead", style: "flex:1;min-width:240px", html: `${ratioText}<br><span class="muted">${pct(s.percentile)} of all ${(d.sequence.length * 19).toLocaleString()} possible mutations in this protein look even more unusual.</span>` }));
  const meters = el("div", { class: "cols" },
    meter("Model", `ESM-2 ${s.model}`, s.percentile, "most unusual", "most natural"),
    lab ? meter("Lab", `${d.assay}, score ${lab[0].toFixed(2)}`, lab[1], "most damaging", "works best")
      : el("div", { class: "card" }, el("h3", {}, "Lab"), el("p", { class: "muted" }, Object.keys(d.lab || {}).length ? "This mutation was not tested in the experiment." : "No lab data for this sequence.")));
  let agree = null;
  if (lab) {
    const modelBad = s.percentile < 40, labBad = lab[1] < 50;
    agree = modelBad === labBad
      ? el("div", { class: "agree yes", html: `<b>The model and the lab agree:</b> both rate this mutation as ${labBad ? "harmful" : "acceptable"}.` })
      : modelBad
        ? el("div", { class: "agree no", html: "<b>The model is harsher than the lab.</b> ESM-2 scores how natural a change looks to evolution, not whether the protein still works. Rare-in-nature changes get flagged even when the protein copes, or works even better, as with activating cancer mutations like HRas G12V." })
        : el("div", { class: "agree no", html: "<b>The lab finds this damaging, but the model barely objects.</b> Some harmful changes still look natural to evolution, for example when the original residue does a job specific to this protein." });
  }
  return [hero, meters, agree, el("p", { class: "muted" }, "Press → to ask why the model gave this score.")];
}

function meter(title, sub, value, left, right) {
  return el("div", { class: "card meter" },
    el("div", { class: "meter-top" }, el("b", {}, title), el("span", { class: "muted", style: "font-size:.82rem" }, sub)),
    el("div", { class: "track" }, el("span", { class: "dot", style: `left:${value}%` })),
    el("div", { class: "track-labels" }, el("span", {}, left), el("span", {}, right)));
}

/* ================================================================ stage 5: why */
async function runExplain() {
  const d = state.protein, req = state.req;
  try {
    const x = await api("/api/explain", { sequence: d.sequence, position: state.pos, mutant: state.mut, model: state.model, protein: d.key });
    if (req !== state.req) return;
    state.explain = x;
  } catch (e) { state.explain = { error: e.message }; }
  if (STAGES[state.stage].key === "why") { renderStage(); updateStrand(); styleStructure(); }
}

function stageWhy() {
  const d = state.protein, pos = state.pos, wt = d.sequence[pos], mt = state.mut;
  head(`Why did the model score ${mname(pos, wt, mt)} that way?`,
    `Integrated gradients traces the score back to every other residue. <b>Arcs on the strand</b> link the hidden site to the 10 residues that mattered most. Red pushed the score down (made ${mt} look worse), blue pushed it up. The same colours are painted on the 3D structure.`);
  const x = state.explain;
  if (!x) return busy(`Tracing the score back through ESM-2 ${state.model} (a few seconds)…`);
  if (x.error) return el("p", { class: "error-text" }, x.error);
  const vmax = Math.max(...x.attribution.map(Math.abs));
  const nLocal = x.top.filter((r) => r.seq_separation <= 4).length;
  const nTouch = x.top.filter((r) => r.contact).length;
  const nFar = x.top.filter((r) => r.contact && r.seq_separation >= 6).length;
  const nums = el("div", { class: "facts compact" },
    el("div", { class: "fact" }, el("b", {}, `${nLocal}/10`), el("span", {}, "within 4 positions in the sequence")),
    x.has_structure ? el("div", { class: "fact" }, el("b", {}, `${nTouch}/10`), el("span", {}, "touch the site in 3D (< 8 Å)")) : null,
    x.has_structure ? el("div", { class: "fact" }, el("b", {}, nFar), el("span", {}, "far in sequence but touching in 3D")) : null);
  const story = x.has_structure
    ? `${nTouch >= 4 ? "Most of what drives the score sits right around the site in the folded protein" : "The score draws on residues both near and far"}. ${nFar ? `<b>${nFar}</b> of them are far apart in the sequence yet touch the site in 3D: the model has picked up how the chain folds.` : "Most influence comes from immediate sequence neighbours."}`
    : "Without a structure we can only see sequence distance, but the pattern is the same: nearby residues dominate.";
  const rows = x.top.map((r, k) => {
    const w = (Math.abs(r.value) / vmax) * 50;
    return el("div", { class: "rank-row", onmouseenter: () => hoverResidue(r.position), onmouseleave: () => hoverResidue(null), onclick: () => hoverResidue(r.position) },
      el("span", { class: "i" }, k + 1), el("span", { class: "res" }, `${r.aa}${r.position + 1}`),
      el("div", { class: "effect" }, el("i", { style: `width:${w}%; ${r.value < 0 ? `right:50%` : `left:50%`}; background:${diverging(r.value, vmax)}` })),
      el("span", { class: "v" }, fmt(r.value, 3)),
      x.has_structure ? el("span", { class: "d" + (r.contact ? " touch" : "") }, r.contact ? `${r.distance.toFixed(1)} Å · touching` : `${r.distance.toFixed(1)} Å`) : el("span", { class: "d" }, `${r.seq_separation} away`));
  });
  const note = x.converged ? null : el("p", { class: "muted", style: "font-size:.82rem" }, `Approximate: the integral was stopped at ${x.steps} steps, with a gap of ${Math.abs(x.completeness_gap).toFixed(2)}.`);
  return el("div", { class: "why-grid" },
    el("div", { style: "display:flex;flex-direction:column;gap:14px" }, nums, el("div", { class: "callout", html: story }), note),
    el("div", { class: "rank" }, el("span", { class: "kicker" }, "The 10 residues that mattered most"), rows));
}

/* ================================================================ stage 6: inside */
async function runInside() {
  const d = state.protein, req = state.req;
  try {
    const x = await api("/api/inside", { sequence: d.sequence, position: state.pos, mutant: state.mut, model: state.model });
    if (req !== state.req) return;
    state.inside = x;
  } catch (e) { state.inside = { error: e.message }; }
  if (STAGES[state.stage].key === "inside") renderStage();
}

function stageInside() {
  const d = state.protein, pos = state.pos, wt = d.sequence[pos], mt = state.mut;
  head("Inside the network",
    `Two views of where the decision happens. <b>Left:</b> the score read out after every layer (the "logit lens"). <b>Right:</b> what each of the ${state.model === "35M" ? 240 : 120} attention heads contributes, measured by switching it off.`);
  const x = state.inside;
  if (!x) return busy(`Reading every layer and switching off ${state.model === "35M" ? 240 : 120} heads one at a time…`);
  if (x.error) return el("p", { class: "error-text" }, x.error);
  const final = x.layers[x.layers.length - 1];
  const settle = x.layers.findIndex((v, i) => i > 0 && Math.abs(v) >= 0.5 * Math.abs(final) && Math.sign(v) === Math.sign(final));
  const contrib = x.heads.map((row) => row.map((v) => -v));
  const flat = contrib.flatMap((row, l) => row.map((v, h) => ({ l, h, v }))).sort((a, b) => Math.abs(b.v) - Math.abs(a.v));
  const contactSet = new Set(x.known_heads.contact.map(([l, h]) => `${l}-${h}`));
  const nContact = flat.slice(0, 10).filter((t) => contactSet.has(`${t.l}-${t.h}`)).length;
  const summary = `For ${mname(pos, wt, mt)} the score is near zero for most of the network and ${settle > 0 ? `takes shape at <b>layer ${settle}</b> of ${x.n_layers}` : "is there from the start"}. The heads that matter most are <b>${flat.slice(0, 3).map((t) => `L${t.l}H${t.h}`).join(", ")}</b>. ${nContact ? `${nContact} of the top 10` : "None of the top 10"} are 3D-contact heads.`;
  return [
    el("div", { class: "callout", html: summary }),
    el("div", { class: "inside-grid" },
      el("div", { class: "card" }, el("h3", {}, "Score after each layer"), layerChart(x)),
      el("div", { class: "card" }, el("h3", {}, "What each attention head contributes"), headLegend(), headChart(x, contrib))),
  ];
}

function niceTicks(a, b, n) {
  const span = b - a, mag = Math.pow(10, Math.floor(Math.log10(span / n)));
  const step = [1, 2, 5, 10].map((s) => s * mag).find((s) => span / s <= n + 1);
  const out = [];
  for (let t = Math.ceil(a / step) * step; t <= b + 1e-9; t += step) out.push(Math.abs(t) < 1e-9 ? 0 : t);
  return out;
}

function layerChart(x) {
  const W = 440, H = 250, m = { l: 40, r: 16, t: 22, b: 36 }, n = x.layers.length;
  const lo = Math.min(0, ...x.layers), hi = Math.max(0, ...x.layers), pad = (hi - lo) * 0.12 || 1;
  const y0 = lo - pad, y1 = hi + pad;
  const sx = (i) => m.l + (i / (n - 1)) * (W - m.l - m.r), sy = (v) => m.t + (1 - (v - y0) / (y1 - y0)) * (H - m.t - m.b);
  const g = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: "100%", role: "img", "aria-label": "Score after each layer" });
  for (const t of niceTicks(y0, y1, 4)) {
    g.append(svg("line", { x1: m.l, x2: W - m.r, y1: sy(t), y2: sy(t), class: "grid", "stroke-width": 0.7 }));
    g.append(svg("text", { x: m.l - 8, y: sy(t) + 4, "text-anchor": "end" }, Number.isInteger(t) ? t : t.toFixed(1)));
  }
  g.append(svg("line", { x1: m.l, x2: W - m.r, y1: sy(0), y2: sy(0), class: "zero" }));
  for (let i = 0; i < n; i++) g.append(svg("text", { x: sx(i), y: H - 16, "text-anchor": "middle" }, i));
  g.append(svg("text", { x: W - m.r, y: H - 2, "text-anchor": "end" }, "layer →"));
  g.append(svg("text", { x: 2, y: 12 }, "score"));
  const line = x.layers.map((v, i) => `${i ? "L" : "M"}${sx(i).toFixed(1)},${sy(v).toFixed(1)}`).join("");
  const area = `${line}L${sx(n - 1)},${sy(0)}L${sx(0)},${sy(0)}Z`;
  g.append(svg("path", { d: area, fill: cssVar("--primary"), opacity: 0.12 }));
  g.append(svg("path", { d: line, fill: "none", stroke: cssVar("--primary"), "stroke-width": 2.5, "stroke-linejoin": "round" }));
  x.layers.forEach((v, i) => {
    const last = i === n - 1;
    g.append(svg("circle", { cx: sx(i), cy: sy(v), r: last ? 6 : 4, fill: last ? cssVar("--mask") : cssVar("--primary"), stroke: cssVar("--surface-2"), "stroke-width": 2 }));
    const hit = svg("rect", { x: sx(i) - 14, y: m.t, width: 28, height: H - m.t - m.b, fill: "transparent" });
    hit.addEventListener("mousemove", (e) => showTip(`<b>${i === 0 ? "Before layer 1" : `After layer ${i}`}</b><br>score ${fmt(v, 2)}<br>top guess at this depth: ${x.layer_top_guess[i]}`, e.clientX, e.clientY));
    hit.addEventListener("mouseleave", hideTip);
    g.append(hit);
  });
  g.append(svg("text", { x: sx(n - 1) - 10, y: sy(x.layers[n - 1]) + (x.layers[n - 1] < 0 ? 20 : -12), "text-anchor": "end", style: `fill:${cssVar("--fg")}` }, `final ${fmt(x.layers[n - 1], 2)}`));
  return g;
}

function headLegend() {
  return el("div", { class: "legend" },
    el("span", {}, el("i", { style: `background:${cssVar("--down")}` }), "pushes the score down"),
    el("span", {}, el("i", { style: `background:${cssVar("--up")}` }), "pushes it up"),
    el("span", {}, el("i", { style: `background:transparent; box-shadow: inset 0 0 0 2px ${cssVar("--contact")}` }), "3D-contact head"),
    el("span", {}, el("i", { style: `background:${cssVar("--fg")}; width:6px; height:6px` }), "key head across proteins"));
}

function headChart(x, contrib) {
  const cell = x.n_heads >= 20 ? 20 : 24, m = { l: 34, t: 6, b: 24 };
  const W = m.l + x.n_heads * cell + 4, H = m.t + x.n_layers * cell + m.b;
  const vmax = Math.max(0.05, ...contrib.flat().map(Math.abs));
  const contact = new Set(x.known_heads.contact.map(([l, h]) => `${l}-${h}`));
  const key = new Set(x.known_heads.mutation.map(([l, h]) => `${l}-${h}`));
  const g = svg("svg", { viewBox: `0 0 ${W} ${H}`, width: "100%", role: "img", "aria-label": "Contribution of each attention head" });
  for (let l = 0; l < x.n_layers; l++) {
    g.append(svg("text", { x: m.l - 6, y: m.t + l * cell + cell * 0.66, "text-anchor": "end" }, `L${l}`));
    for (let h = 0; h < x.n_heads; h++) {
      const v = contrib[l][h], id = `${l}-${h}`;
      const r = svg("rect", { x: m.l + h * cell + 1.5, y: m.t + l * cell + 1.5, width: cell - 3, height: cell - 3, rx: 5, fill: diverging(v, vmax) });
      if (contact.has(id)) { r.setAttribute("stroke", cssVar("--contact")); r.setAttribute("stroke-width", 2.2); }
      r.addEventListener("mousemove", (e) => showTip(`<b>Layer ${l}, head ${h}</b><br>contribution ${fmt(v, 3)}<br>switched off, the score becomes ${fmt(x.layers[x.layers.length - 1] - v, 2)}${contact.has(id) ? "<br>a 3D-contact head" : ""}${key.has(id) ? "<br>one of the key heads across proteins" : ""}`, e.clientX, e.clientY));
      r.addEventListener("mouseleave", hideTip);
      g.append(r);
      if (key.has(id)) g.append(svg("circle", { cx: m.l + h * cell + cell / 2, cy: m.t + l * cell + cell / 2, r: 2.6, fill: cssVar("--fg"), "pointer-events": "none" }));
    }
  }
  for (let h = 0; h < x.n_heads; h += 5) g.append(svg("text", { x: m.l + h * cell + cell / 2, y: H - 8, "text-anchor": "middle" }, `H${h}`));
  return g;
}

/* ================================================================ stage 7: map */
async function runLandscape() {
  const d = state.protein;
  try { state.land = await api("/api/landscape", { sequence: d.sequence, model: state.model }); }
  catch (e) { state.land = { error: e.message }; }
  if (state.protein === d && STAGES[state.stage].key === "map") renderStage();
}

function spearman(a, b) {
  const rank = (v) => { const idx = v.map((x, i) => [x, i]).sort((p, q) => p[0] - q[0]); const r = new Array(v.length); let i = 0; while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++; for (let k = i; k <= j; k++) r[idx[k][1]] = (i + j) / 2; i = j + 1; } return r; };
  const ra = rank(a), rb = rank(b), n = a.length, mu = (n - 1) / 2;
  let num = 0, da = 0, db = 0;
  for (let i = 0; i < n; i++) { num += (ra[i] - mu) * (rb[i] - mu); da += (ra[i] - mu) ** 2; db += (rb[i] - mu) ** 2; }
  return num / Math.sqrt(da * db);
}

function stageMap() {
  const d = state.protein;
  head("Every possible mutation at once",
    `ESM-2 scores all 19 replacements at every position. Below it, the lab's measurements for the same protein. <b>Darker means more damaging.</b> The more the two maps look alike, the better the model. Click any cell to explore that mutation.`);
  if (!state.land) return busy(`ESM-2 ${state.model} is scoring all ${(d.sequence.length * 19).toLocaleString()} possible mutations…`);
  if (state.land.error) return el("p", { class: "error-text" }, state.land.error);
  const L = d.sequence.length, llr = state.land.llr, cw = L > 300 ? 5 : L > 150 ? 7 : 10, ch = 9;
  const vals = [];
  for (let i = 0; i < L; i++) for (let a = 0; a < 20; a++) if (AA[a] !== d.sequence[i]) vals.push(llr[i][a]);
  const sorted = [...vals].sort((p, q) => p - q);
  const pctOf = (v) => { let lo = 0, hi = sorted.length; while (lo < hi) { const md = (lo + hi) >> 1; if (sorted[md] < v) lo = md + 1; else hi = md; } return (lo / sorted.length) * 100; };
  const hasLab = Object.keys(d.lab || {}).length > 0;
  const makeMap = (title, color) => {
    const c = el("canvas", { class: "map", "aria-label": title }), dpr = devicePixelRatio || 1;
    c.width = L * cw * dpr; c.height = (20 * ch + 14) * dpr; c.style.width = `${L * cw}px`; c.style.height = `${20 * ch + 14}px`;
    const ctx = c.getContext("2d"); ctx.scale(dpr, dpr);
    ctx.fillStyle = cssVar("--muted"); ctx.font = `9px ${cssVar("--f-code")}`; ctx.textAlign = "center";
    for (let i = 9; i < L; i += 10) ctx.fillText(String(i + 1), i * cw + cw / 2, 9);
    for (let i = 0; i < L; i++) for (let r = 0; r < 20; r++) {
      const aa = WHEEL_ORDER[r];
      if (aa === d.sequence[i]) { ctx.fillStyle = cssVar("--fg"); ctx.beginPath(); ctx.arc(i * cw + cw / 2, 14 + r * ch + ch / 2, Math.min(cw, ch) * 0.22, 0, 7); ctx.fill(); continue; }
      ctx.fillStyle = color(i, aa);
      ctx.fillRect(i * cw, 14 + r * ch, cw - (cw > 5 ? 1 : 0), ch - 1);
    }
    if (state.pos !== null) { ctx.strokeStyle = cssVar("--mask"); ctx.lineWidth = 2; ctx.strokeRect(state.pos * cw - 1, 13, cw + 1, 20 * ch + 1); }
    const cell = (e) => { const r = c.getBoundingClientRect(); return [Math.floor((e.clientX - r.left) / cw), Math.floor((e.clientY - r.top - 14) / ch)]; };
    c.addEventListener("mousemove", (e) => {
      const [i, r] = cell(e);
      if (i < 0 || i >= L || r < 0 || r >= 20) return hideTip();
      const aa = WHEEL_ORDER[r], wt = d.sequence[i];
      if (aa === wt) return showTip(`<b>${wt}${i + 1}</b> original residue`, e.clientX, e.clientY);
      const v = llr[i][AA.indexOf(aa)], lab = d.lab?.[mname(i, wt, aa)];
      showTip(`<b>${mname(i, wt, aa)}</b><br>model ${fmt(v, 2)} · more unusual than ${pct(100 - pctOf(v))} of mutations${lab ? `<br>lab: did better than ${pct(lab[1])} of tested mutations` : hasLab ? "<br>lab: not tested" : ""}`, e.clientX, e.clientY);
    });
    c.addEventListener("mouseleave", hideTip);
    c.addEventListener("click", async (e) => {
      const [i, r] = cell(e);
      if (i < 0 || i >= L || r < 0 || r >= 20 || WHEEL_ORDER[r] === d.sequence[i]) return;
      hideTip();
      await selectPosition(i, { stay: true });
      await chooseMutant(WHEEL_ORDER[r]);
    });
    const labels = el("div", { class: "map-labels", style: "padding-top:14px" }, ...[...WHEEL_ORDER].map((a) => el("span", { style: `height:${ch}px;line-height:${ch}px` }, a)));
    return el("div", {}, el("div", { class: "map-title" }, title), el("div", { class: "map-row" }, labels, c));
  };
  const maps = [makeMap(`Model · ESM-2 ${state.model}`, (i, aa) => ramp(1 - pctOf(llr[i][AA.indexOf(aa)]) / 100))];
  let agree = null;
  if (hasLab) {
    maps.push(makeMap(`Lab · ${d.assay} (${d.assay_type.toLowerCase()})`, (i, aa) => { const lab = d.lab[mname(i, d.sequence[i], aa)]; return lab ? ramp(1 - lab[1] / 100) : cssVar("--surface-3"); }));
    const m = [], l = [];
    for (const [mut, [score]] of Object.entries(d.lab)) { const i = parseInt(mut.slice(1, -1), 10) - 1, a = AA.indexOf(mut.slice(-1)); if (i >= 0 && i < L && a >= 0) { m.push(llr[i][a]); l.push(score); } }
    agree = el("div", { class: "callout", html: `<b>Agreement for ${d.name}:</b> Spearman ρ = <b>${spearman(m, l).toFixed(2)}</b> over ${m.length.toLocaleString()} tested mutations (0 = no relation, 1 = perfect). Across all 217 experiments, ESM-2 35M averages 0.32.` });
  }
  return [
    el("div", { class: "legend" }, "harmless", el("span", { class: "ramp" }), "damaging", el("span", {}, "· dot = original residue · grey = not tested · rows grouped by chemistry")),
    el("div", { class: "map-scroll" }, maps), agree,
  ];
}

/* ================================================================ selection */
async function selectPosition(pos, { stay = false } = {}) {
  const d = state.protein, req = ++state.req;
  Object.assign(state, { pos, mut: null, predict: null, score: null, explain: null, inside: null });
  updateStrand(); scrollStrandTo(pos); styleStructure(); focusStructure(); renderSiteInfo();
  if (!stay) go(state.stage <= 1 ? 1 : 2); else renderRail();
  try {
    const p = await api("/api/predict", { sequence: d.sequence, position: pos, model: state.model });
    if (req !== state.req) return;
    state.predict = p;
    if (STAGES[state.stage].key === "guess") renderStage();
  } catch (e) { if (req === state.req) $("#stage-body").replaceChildren(el("p", { class: "error-text" }, e.message)); }
}

async function chooseMutant(aa, { stay = false } = {}) {
  const d = state.protein, pos = state.pos;
  if (pos === null || aa === d.sequence[pos]) return;
  const req = ++state.req;
  Object.assign(state, { mut: aa, score: null, explain: null, inside: null });
  if (!stay) go(3); else renderRail();
  try {
    const s = await api("/api/score", { sequence: d.sequence, position: pos, mutant: aa, model: state.model });
    if (req !== state.req) return;
    state.score = s;
    if (STAGES[state.stage].key === "score") renderStage();
  } catch (e) { if (req === state.req) $("#stage-body").replaceChildren(el("p", { class: "error-text" }, e.message)); }
}

/* ================================================================ strand */
function renderStrand() {
  const d = state.protein, L = d.sequence.length;
  const strand = $("#strand");
  strand.style.width = `${L * STRIDE + 40}px`;
  $("#arcs").setAttribute("width", L * STRIDE + 40);
  $("#arcs").setAttribute("height", 76);
  $("#beads").replaceChildren(...[...d.sequence].map((aa, i) => el("button", {
    class: "bead", type: "button", "data-pos": i, style: `--c: var(--cls-${AA_CLASS[aa]})`, "aria-label": `${AA_NAME[aa]} ${i + 1}`,
    onclick: () => selectPosition(i),
    onmouseenter: (e) => showTip(beadTip(i), e.clientX, e.clientY), onmouseleave: hideTip,
  }, el("span", { class: "aa" }, aa), el("span", { class: "tid" }, state.tokens[i + 1].id), (i + 1) % 10 === 0 || i === 0 ? el("span", { class: "num" }, i + 1) : null)));
  $("#strand-scroll").scrollLeft = 0;
  updateStrand();
}
$("#show-ids").addEventListener("change", (e) => { state.showIds = e.target.checked; $("#strand").classList.toggle("show-ids", state.showIds); });

function heatOn() { return STAGES[state.stage].key === "why" && state.explain && !state.explain.error && state.explain.position === state.pos; }

function updateStrand() {
  if (!state.protein) return;
  const imp = heatOn() ? state.explain.attribution : null;
  const vmax = imp ? Math.max(...imp.map(Math.abs)) : 1;
  const topSet = new Set(imp ? state.explain.top.map((r) => r.position) : []);
  document.querySelectorAll("#beads .bead").forEach((b) => {
    const i = +b.dataset.pos, masked = i === state.pos;
    b.classList.toggle("masked", masked);
    b.querySelector(".aa").textContent = masked ? "?" : state.protein.sequence[i];
    b.querySelector(".tid").textContent = masked ? 32 : state.tokens[i + 1].id;
    b.classList.toggle("heat", !!imp && !masked);
    b.classList.toggle("top", topSet.has(i));
    b.style.setProperty("--heat", imp && !masked ? diverging(imp[i], vmax) : "");
  });
  drawArcs(imp, vmax);
  const legend = $("#strand-legend");
  if (imp) {
    legend.replaceChildren(
      el("span", {}, el("i", { style: `background:${cssVar("--down")}` }), "pushed the score down"),
      el("span", {}, el("i", { style: `background:${cssVar("--up")}` }), "pushed it up"),
      state.explain.has_structure ? el("span", {}, el("i", { class: "line", style: `background:${cssVar("--fg")}` }), "touches the site in 3D") : null,
      state.explain.has_structure ? el("span", {}, el("i", { class: "dash" }), "does not touch") : null);
  } else {
    legend.replaceChildren(...Object.entries(CLASS_LABEL).map(([k, v]) => el("span", {}, el("i", { style: `background:var(--cls-${k})` }), v)),
      el("span", { style: "color:var(--faint)" }, state.pos === null ? "· click a bead to hide it" : `· position ${state.pos + 1} is hidden`));
  }
}

function drawArcs(imp, vmax) {
  const g = $("#arcs");
  g.replaceChildren();
  if (!imp) return;
  const cx = (i) => i * STRIDE + 13.5;
  const x1 = cx(state.pos);
  for (const r of [...state.explain.top].reverse()) {
    const x2 = cx(r.position), dx = Math.abs(x2 - x1), h = Math.min(56, 16 + dx * 0.22);
    const path = `M${x1},44 C${x1},${44 - h} ${x2},${64 - h} ${x2},64`;
    const w = 1.5 + 4.5 * (Math.abs(r.value) / vmax);
    const arc = svg("path", { d: path, class: "arc", stroke: diverging(r.value, vmax), "stroke-width": w, "stroke-dasharray": state.explain.has_structure && !r.contact ? "5 5" : null });
    arc.addEventListener("mousemove", (e) => showTip(`<b>${r.aa}${r.position + 1}</b> effect ${fmt(r.value, 3)}<br>${r.seq_separation} positions away${r.distance !== undefined ? `, ${r.distance.toFixed(1)} Å in 3D${r.contact ? " (touching)" : ""}` : ""}`, e.clientX, e.clientY));
    arc.addEventListener("mouseleave", hideTip);
    g.append(arc);
  }
}

function scrollStrandTo(pos) {
  const sc = $("#strand-scroll");
  sc.scrollTo({ left: Math.max(0, pos * STRIDE - sc.clientWidth / 2), behavior: "smooth" });
}

function beadTip(i) {
  const d = state.protein, aa = d.sequence[i], r = d.residues?.[i];
  let html = `<b>${aa}${i + 1}</b> · ${AA_NAME[aa]}<br>${CLASS_LABEL[AA_CLASS[aa]]}`;
  if (r) html += ` · ${SS_NAME[r.ss] || "loop"} · ${r.rsa < 0.25 ? "buried" : "surface"}`;
  if (heatOn() && i !== state.pos) html += `<br>effect on the score: ${fmt(state.explain.attribution[i], 3)}`;
  if (state.showIds) html += `<br>token ${state.tokens[i + 1].id}`;
  return html;
}

function hoverResidue(pos) {
  document.querySelectorAll("#beads .bead").forEach((b) => b.classList.toggle("hl", +b.dataset.pos === pos));
  highlighted = pos;
  styleStructure();
  if (pos !== null) scrollStrandTo(pos);
}

/* ================================================================ 3D structure */
let viewer = null, highlighted = null, mouse = [0, 0];
function initViewer() {
  if (!window.$3Dmol) { $("#viewer").replaceChildren(el("div", { class: "viewer-empty" }, "The 3D viewer could not load.")); return; }
  viewer = $3Dmol.createViewer($("#viewer"), { backgroundColor: cssVar("--surface"), antialias: true });
  $("#viewer").addEventListener("mousemove", (e) => { mouse = [e.clientX, e.clientY]; });
  addEventListener("resize", resizeViewer);
}
function resizeViewer() { if (viewer) setTimeout(() => { viewer.resize(); viewer.render(); }, 40); }

function showStructure() {
  if (!viewer) return;
  const d = state.protein;
  viewer.clear();
  $("#viewer-note").textContent = d.pdb ? `AlphaFold2 model of ${d.name}` : "no structure for this sequence";
  if (!d.pdb) { viewer.render(); styleStructure(); return; }
  viewer.addModel(d.pdb, "pdb");
  viewer.setHoverable({}, true, (atom) => { if (atom) showTip(beadTip(atom.resi - 1), mouse[0], mouse[1]); }, hideTip);
  viewer.setClickable({}, true, (atom) => { hideTip(); selectPosition(atom.resi - 1); });
  const core = (d.residues || []).map((r, i) => (r.plddt >= 70 ? i + 1 : null)).filter((x) => x !== null);
  viewer.zoomTo(core.length > 10 ? { resi: core } : {});
  styleStructure();
}

function styleStructure() {
  if (!viewer || !state.protein) return;
  viewer.setBackgroundColor(cssVar("--surface"));
  const legend = $("#viewer-legend");
  if (!state.protein.pdb) { legend.replaceChildren(el("span", {}, "Paste-in sequences have no 3D structure.")); viewer.render(); return; }
  const base = mix(cssVar("--surface"), cssVar("--muted"), 0.5);
  const imp = heatOn() ? state.explain.attribution : null;
  const vmax = imp ? Math.max(...imp.map(Math.abs)) : 1;
  viewer.setStyle({}, { cartoon: imp ? { colorfunc: (a) => (a.resi - 1 === state.pos ? cssVar("--mask") : diverging(imp[a.resi - 1] || 0, vmax)) } : { color: base } });
  if (imp) viewer.addStyle({ resi: state.explain.top.map((r) => r.position + 1) }, { stick: { colorfunc: (a) => diverging(imp[a.resi - 1], vmax), radius: 0.24 } });
  if (state.pos !== null) viewer.addStyle({ resi: state.pos + 1 }, { stick: { color: cssVar("--mask"), radius: 0.32 }, sphere: { color: cssVar("--mask"), radius: 1.5, opacity: 0.6 } });
  if (highlighted !== null) viewer.addStyle({ resi: highlighted + 1 }, { stick: { color: cssVar("--fg"), radius: 0.34 } });
  viewer.render();
  legend.replaceChildren(...(imp
    ? [el("span", {}, el("i", { style: `background:${cssVar("--mask")}` }), "hidden site"), el("span", {}, el("i", { style: `background:${cssVar("--down")}` }), "pushed the score down"), el("span", {}, el("i", { style: `background:${cssVar("--up")}` }), "pushed it up")]
    : [el("span", {}, el("i", { style: `background:${cssVar("--mask")}` }), state.pos !== null ? "hidden site" : "click the protein or a bead"), el("span", { style: "color:var(--faint)" }, "drag to rotate · scroll to zoom")]));
}

function focusStructure() {
  if (viewer && state.protein?.pdb && state.pos !== null) viewer.zoomTo({ within: { distance: 22, sel: { resi: state.pos + 1 } } }, 700);
}

function renderSiteInfo() {
  const d = state.protein, box = $("#site-info");
  if (!d || state.pos === null) { box.replaceChildren(); return; }
  const i = state.pos, aa = d.sequence[i], r = d.residues?.[i];
  const chips = [el("span", { class: "chip" }, el("b", {}, `${aa}${i + 1}`), ` · ${AA_NAME[aa]}`)];
  if (r) {
    chips.push(el("span", { class: "chip" }, SS_NAME[r.ss] || "loop"));
    chips.push(el("span", { class: "chip" }, r.rsa < 0.25 ? "buried inside" : "on the surface"));
    if (r.conservation !== null) chips.push(el("span", { class: "chip" }, r.conservation > 0.6 ? "highly conserved" : r.conservation > 0.35 ? "moderately conserved" : "variable across species"));
    chips.push(el("span", { class: "chip" }, `AlphaFold confidence ${r.plddt.toFixed(0)}`));
  }
  box.replaceChildren(...chips);
}

/* ================================================================ findings */
const FINDINGS = [
  { fig: "fig1_baseline", tag: "solid", title: "Our scores reproduce the benchmark exactly", stat: "<b>217 / 217</b> experiments match ProteinGym", text: "Before explaining the model we checked that we compute its scores correctly. Mean Spearman with lab data: 0.32 for ESM-2 35M and 0.21 for 8M. The model is best on stability experiments, where structure matters most. On 62,727 clinical variants, 35M separates disease-causing from harmless with AUROC 0.74." },
  { fig: "fig2_probing", tag: "solid", title: "The model learned protein structure on its own", stat: "Helix / strand / loop: <b>0.80</b> vs 0.45 from the letter alone", text: "Simple probes read secondary structure, burial, solvent exposure, contact count and conservation from the model's internal numbers on proteins it never saw. The signal grows with depth and is absent in an untrained copy, so it came from training on sequences only." },
  { fig: "fig3_attention_contacts", tag: "solid", title: "A few late attention heads detect 3D contacts", stat: "Head L11H13: <b>11×</b> chance (0.39 vs 0.035)", text: "Some heads link residues that touch in the folded protein. They sit in the last layers. No head of an untrained model does this." },
  { fig: "fig4_attribution", tag: "care", title: "Explanations point at the site's neighbours", stat: "AUROC <b>0.74</b> for 3D contacts (untrained model 0.46–0.52)", text: "Integrated gradients finds residues that physically surround a mutation, and hiding the top one moves the score 10× more than hiding a random one. But simply taking the nearest residues in the sequence works just as well: the model leans heavily on local context." },
  { fig: "fig5_ablation", tag: "solid", title: "The score depends on a few heads, not the contact heads", stat: "Top 16 heads: accuracy drops <b>0.28</b>; 16 random heads: 0.05", text: "Across 24 proteins, heads were ranked on half of the positions and tested on the other half. Switching off the most important ones hurts in 23 of 24 proteins. Switching off the 16 best 3D-contact heads costs only 0.04, about the same as random heads: seeing the fold and judging a mutation are done by different parts of the model." },
  { fig: "fig6_logit_lens", tag: "solid", title: "The verdict forms in the last layers", stat: "Near zero until layer 8 of 12, then <b>0.34</b>", text: "Reading the model's guess after every layer shows that the early and middle layers build features, and the last third turns them into a judgement about the mutation." },
  { fig: "fig7_evolution_vs_function", tag: "solid", title: "The model knows what breaks a protein, not what improves it", stat: "Ranking within the worst quarter <b>0.18</b> vs best quarter 0.04", text: "Across 209 experiments the model separates damaging mutations well but barely ranks the beneficial ones (74% of experiments, 88% of stability ones). Evolution shows what is forbidden, not what is better. HRas cancer mutations are the extreme case: they make the protein over-active, yet the model calls them unusual. Part of the gap may be measurement noise among near-neutral mutations." },
  { fig: "fig8_head_roles", tag: "solid", title: "The heads that judge mutations look nearby", stat: "<b>60%</b> of their attention on residues 1–4 positions away (others 23%)", text: "From the hidden site, the heads that matter most look only about 5 residues away on average, against 25 for other heads, and none of them is a 3D-contact head. The verdict is built mainly from local sequence context, which is why sequence neighbours explain the score as well as our attributions do." },
];
function renderFindings() {
  $("#findings").replaceChildren(...FINDINGS.map((f, k) => el("article", { class: "finding" },
    el("div", {}, el("span", { class: "fnum" }, String(k + 1).padStart(2, "0")), el("h2", {}, f.title), el("p", { class: "stat", html: f.stat }), el("p", { class: "body" }, f.text),
      el("span", { class: `tag ${f.tag}` }, f.tag === "solid" ? "Strong result" : "Good, with a catch")),
    el("div", { class: "plate" }, el("img", { src: `/figures/${f.fig}.png`, alt: f.title, loading: "lazy" })))));
}

init();
