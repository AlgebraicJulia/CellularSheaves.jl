// The certified corridor planner demo: renders the data written by
// docs/scripts/tropical_viz_export.jl. Built with esbuild (docs/viz/build.sh).

type State = [number, number, number, number];          // qx, qy, vx, vy (1-based positions)
type RoutePoint = [number, number, number, number];     // qx, qy, heading angle, speed band (-1: no mode)

interface Run {
  planner: "blind" | "modes";
  width: number; mode_width: number; used_width: number; used_mode_width: number;
  first_try: boolean; cost: number; bound: number | null; lower: number | null; certified: boolean;
  route: RoutePoint[]; corridor: string; traj: State[];
  V_U: number[]; V_opt: number[]; lower_along: (number | null)[];
}
interface Query {
  id: number; src: State; cost: number; optimal: State[]; V_opt: number[];
  lower_blind: number[]; lower_modes: number[];
  explored: { dijkstra: string; blind: string; modes: string };
  expanded: { dijkstra: number; blind: number; modes: number };
  runs: Run[];
}
interface Row { [key: string]: string | number }
interface Data {
  n: number; vmax: number; headings: number; goal: [number, number][]; walls: string;
  V_rest: string; lower_rest: { blind: string; modes: string };
  states: number; edges: number; reachable: number; queries: Query[];
  phase: { nq: number; goal: number; value: number[]; accel: number[]; trajectories: [number, number][][] };
  layered: Row[]; certified: Row[];
}

const data: Data = JSON.parse((document.getElementById("demo-data") as HTMLScriptElement).textContent || "{}");
const N = data.n;

// ---- decoding --------------------------------------------------------------------
function bytes(b64: string): Uint8Array {
  const s = atob(b64);
  const out = new Uint8Array(s.length);
  for (let i = 0; i < s.length; i++) out[i] = s.charCodeAt(i);
  return out;
}
function bits(b64: string, len: number): Uint8Array {
  const raw = bytes(b64);
  const out = new Uint8Array(len);
  for (let i = 0; i < len; i++) out[i] = (raw[i >> 3] >> (i & 7)) & 1;
  return out;
}
const walls = bits(data.walls, N * N);
const vRest = bytes(data.V_rest);                       // 255: cannot reach the goal
const lowerRest = { blind: bytes(data.lower_rest.blind), modes: bytes(data.lower_rest.modes) };

// ---- theme ----------------------------------------------------------------------------
function css(name: string): string {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}
function isDark(): boolean {
  const t = document.documentElement.getAttribute("data-theme");
  if (t) return t === "dark";
  return window.matchMedia("(prefers-color-scheme: dark)").matches;
}
// Sequential ramps (light → dark), from the reference palette.
const BLUE = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"];
const ORANGE = ["#fde3d6", "#fbcdb6", "#f8b391", "#f4986d", "#ef7f4c", "#eb6834", "#d95926", "#bf4c1f", "#9f3e19", "#7f3114"];
function hexRgb(h: string): [number, number, number] {
  const v = parseInt(h.slice(1), 16);
  return [(v >> 16) & 255, (v >> 8) & 255, v & 255];
}
function ramp(colors: string[], t: number, reverse: boolean): [number, number, number] {
  const c = reverse ? [...colors].reverse() : colors;
  const x = Math.min(1, Math.max(0, t)) * (c.length - 1);
  const i = Math.min(c.length - 2, Math.floor(x));
  const a = hexRgb(c[i]), b = hexRgb(c[i + 1]), f = x - i;
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}

// ---- marching squares for isochrones ----------------------------------------------
type Seg = [number, number, number, number];
function contours(val: (x: number, y: number) => number, w: number, h: number, level: number): Seg[] {
  const segs: Seg[] = [];
  for (let y = 0; y < h - 1; y++) for (let x = 0; x < w - 1; x++) {
    const v = [val(x, y), val(x + 1, y), val(x + 1, y + 1), val(x, y + 1)];
    if (v.some(z => !isFinite(z))) continue;
    let idx = 0;
    v.forEach((z, k) => { if (z >= level) idx |= 1 << k; });
    if (idx === 0 || idx === 15) continue;
    const p = (a: number, b: number, ax: number, ay: number, bx: number, by: number): [number, number] => {
      const t = (level - a) / (b - a);
      return [ax + (bx - ax) * t, ay + (by - ay) * t];
    };
    const e: [number, number][] = [
      p(v[0], v[1], x, y, x + 1, y), p(v[1], v[2], x + 1, y, x + 1, y + 1),
      p(v[3], v[2], x, y + 1, x + 1, y + 1), p(v[0], v[3], x, y, x, y + 1)];
    const edges: Record<number, number[][]> = {
      1: [[0, 3]], 2: [[0, 1]], 3: [[1, 3]], 4: [[1, 2]], 5: [[0, 1], [2, 3]], 6: [[0, 2]], 7: [[2, 3]],
      8: [[2, 3]], 9: [[0, 2]], 10: [[0, 3], [1, 2]], 11: [[1, 2]], 12: [[1, 3]], 13: [[0, 1]], 14: [[0, 3]] };
    for (const [a, b] of edges[idx]) segs.push([e[a][0], e[a][1], e[b][0], e[b][1]]);
  }
  return segs;
}

// ---- state -----------------------------------------------------------------------------
const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const qSel = $<HTMLSelectElement>("q"), pSel = $<HTMLSelectElement>("planner"), cSel = $<HTMLSelectElement>("corr");
const lSel = $<HTMLSelectElement>("layer");
const toggles = ["t-corr", "t-route", "t-ctl", "t-opt", "t-iso"].map(id => $<HTMLInputElement>(id));
data.queries.forEach((q, i) => {
  const o = document.createElement("option");
  const [x, y, vx, vy] = q.src;
  o.value = String(i);
  o.textContent = `#${q.id}  q = (${x}, ${y}), v = (${vx}, ${vy}) · optimum ${q.cost} steps`;
  qSel.appendChild(o);
});
function query(): Query { return data.queries[+qSel.value]; }
function run(): Run {
  const q = query();
  return q.runs.find(r => r.planner === pSel.value && r.width === +cSel.value) as Run;
}

// ---- workspace map ------------------------------------------------------------------
const map = $<HTMLCanvasElement>("map");
const mctx = map.getContext("2d") as CanvasRenderingContext2D;
const S = map.width / N;                                // pixels per cell
const cx = (qx: number) => (qx - 0.5) * S;
const cy = (qy: number) => (N - qy + 0.5) * S;          // qy upward

function layerValues(): { v: Float64Array; max: number; kind: "time" | "count" } {
  const layer = lSel.value;
  const v = new Float64Array(N * N);
  if (layer === "value" || layer === "lower") {
    const src = layer === "value" ? vRest : lowerRest[pSel.value as "blind" | "modes"];
    let max = 0;
    for (let i = 0; i < N * N; i++) { v[i] = src[i] === 255 ? Infinity : src[i]; if (isFinite(v[i])) max = Math.max(max, v[i]); }
    return { v, max, kind: "time" };
  }
  const which = layer.slice(4) as "dijkstra" | "blind" | "modes";
  const src = bytes(query().explored[which]);
  let max = 0;
  for (let i = 0; i < N * N; i++) { v[i] = src[i]; max = Math.max(max, v[i]); }
  return { v, max, kind: "count" };
}

function drawMap(): void {
  const dark = isDark();
  const { v, max, kind } = layerValues();
  const img = mctx.createImageData(N, N);
  const wall = hexRgb(css("--wall")), surf = hexRgb(css("--surface"));
  for (let qy = 1; qy <= N; qy++) for (let qx = 1; qx <= N; qx++) {
    const i = (qy - 1) * N + (qx - 1), p = ((N - qy) * N + (qx - 1)) * 4;
    let c: [number, number, number];
    if (walls[i]) c = wall;
    else if (kind === "time") c = isFinite(v[i]) ? ramp(BLUE, v[i] / max, dark) : surf;
    else c = v[i] > 0 ? ramp(ORANGE, Math.sqrt(v[i] / max), dark) : surf;
    img.data[p] = c[0]; img.data[p + 1] = c[1]; img.data[p + 2] = c[2]; img.data[p + 3] = 255;
  }
  const off = document.createElement("canvas");
  off.width = N; off.height = N;
  (off.getContext("2d") as CanvasRenderingContext2D).putImageData(img, 0, 0);
  mctx.imageSmoothingEnabled = false;
  mctx.clearRect(0, 0, map.width, map.height);
  mctx.drawImage(off, 0, 0, map.width, map.height);

  const r = run(), q = query();
  if (toggles[0].checked) {                               // corridor
    const cor = bits(r.corridor, N * N);
    mctx.fillStyle = css("--corridor");
    for (let i = 0; i < N * N; i++) if (cor[i] && !walls[i]) {
      const qx = (i % N) + 1, qy = Math.floor(i / N) + 1;
      mctx.fillRect((qx - 1) * S, (N - qy) * S, S, S);
    }
    // outline of the corridor
    mctx.strokeStyle = css("--s-low"); mctx.lineWidth = 1.2;
    for (let i = 0; i < N * N; i++) if (cor[i]) {
      const qx = (i % N) + 1, qy = Math.floor(i / N) + 1, x0 = (qx - 1) * S, y0 = (N - qy) * S;
      const inC = (x: number, y: number) => x >= 1 && x <= N && y >= 1 && y <= N && cor[(y - 1) * N + (x - 1)] === 1;
      mctx.beginPath();
      if (!inC(qx - 1, qy)) { mctx.moveTo(x0, y0); mctx.lineTo(x0, y0 + S); }
      if (!inC(qx + 1, qy)) { mctx.moveTo(x0 + S, y0); mctx.lineTo(x0 + S, y0 + S); }
      if (!inC(qx, qy + 1)) { mctx.moveTo(x0, y0); mctx.lineTo(x0 + S, y0); }
      if (!inC(qx, qy - 1)) { mctx.moveTo(x0, y0 + S); mctx.lineTo(x0 + S, y0 + S); }
      mctx.stroke();
    }
  }
  if (toggles[4].checked && kind === "time") {          // isochrones every 5 steps
    mctx.strokeStyle = dark ? "rgba(255,255,255,0.45)" : "rgba(11,11,11,0.38)";
    mctx.lineWidth = 1;
    const val = (x: number, y: number) => { const i = y * N + x; return walls[i] ? Infinity : v[i]; };
    for (let level = 5; level < max; level += 5) {
      mctx.beginPath();
      for (const [x0, y0, x1, y1] of contours(val, N, N, level)) {
        mctx.moveTo((x0 + 0.5) * S, (N - y0 - 0.5) * S);
        mctx.lineTo((x1 + 0.5) * S, (N - y1 - 0.5) * S);
      }
      mctx.stroke();
    }
  }
  // goal
  mctx.fillStyle = css("--good");
  for (const [gx, gy] of data.goal) mctx.fillRect((gx - 1) * S, (N - gy) * S, S, S);
  if (toggles[1].checked) {                               // route with mode headings
    mctx.strokeStyle = css("--s-route"); mctx.lineWidth = 2; mctx.setLineDash([6, 4]);
    mctx.beginPath();
    r.route.forEach(([x, y], k) => k === 0 ? mctx.moveTo(cx(x), cy(y)) : mctx.lineTo(cx(x), cy(y)));
    mctx.stroke(); mctx.setLineDash([]);
    mctx.fillStyle = css("--s-route");
    for (const [x, y, ang, band] of r.route) {
      mctx.beginPath(); mctx.arc(cx(x), cy(y), 2.6, 0, 2 * Math.PI); mctx.fill();
      if (band > 0) {                                       // the route's velocity mode
        const L = 4 + 3.5 * band;
        mctx.beginPath(); mctx.moveTo(cx(x), cy(y));
        mctx.lineTo(cx(x) + L * Math.cos(ang), cy(y) - L * Math.sin(ang));
        mctx.strokeStyle = css("--s-route"); mctx.lineWidth = 1.5; mctx.stroke();
      }
    }
  }
  const path = (states: State[], color: string, width: number) => {
    mctx.strokeStyle = color; mctx.lineWidth = width; mctx.lineJoin = "round";
    mctx.beginPath();
    states.forEach(([x, y], k) => k === 0 ? mctx.moveTo(cx(x), cy(y)) : mctx.lineTo(cx(x), cy(y)));
    mctx.stroke();
    mctx.fillStyle = color;
    for (const [x, y] of states) { mctx.beginPath(); mctx.arc(cx(x), cy(y), 3.2, 0, 2 * Math.PI); mctx.fill(); }
  };
  if (toggles[3].checked) path(q.optimal, css("--s-opt"), 2);
  if (toggles[2].checked) path(r.traj, css("--s-ctl"), 2.5);
  // start
  const [sx, sy, svx, svy] = q.src;
  mctx.strokeStyle = css("--ink"); mctx.lineWidth = 2;
  mctx.beginPath(); mctx.arc(cx(sx), cy(sy), 7, 0, 2 * Math.PI); mctx.stroke();
  if (svx !== 0 || svy !== 0) {
    mctx.beginPath(); mctx.moveTo(cx(sx), cy(sy)); mctx.lineTo(cx(sx + svx), cy(sy + svy)); mctx.stroke();
  }

  // legend and note
  const leg = $("map-legend");
  const rampCss = (colors: string[]) => `linear-gradient(90deg, ${(dark ? [...colors].reverse() : colors).join(",")})`;
  const scale = kind === "time"
    ? `<span class="key"><span class="ramp" style="background:${rampCss(BLUE)}"></span><span class="num">0–${max} steps</span></span>`
    : `<span class="key"><span class="ramp" style="background:${rampCss(ORANGE)}"></span><span class="num">0–${max} states per cell</span></span>`;
  leg.innerHTML = scale +
    `<span class="key"><span class="swatch" style="background:var(--s-ctl);height:3px"></span>certified controller</span>` +
    `<span class="key"><span class="swatch" style="background:var(--s-opt)"></span>optimal</span>` +
    `<span class="key"><span class="swatch" style="background:var(--s-route)"></span>route (ticks: mode heading)</span>` +
    `<span class="key"><span class="swatch box" style="background:var(--corridor);outline:1px solid var(--s-low)"></span>corridor</span>` +
    `<span class="key"><span class="swatch box" style="background:var(--good)"></span>goal, stopped</span>` +
    `<span class="key"><span class="swatch box" style="background:var(--wall)"></span>wall</span>`;
  const notes: Record<string, string> = {
    value: "Minimum time to the goal starting at rest; contours every 5 steps.",
    lower: "Coarse value pulled back along the oplax refinement: never above V*.",
    "exp-dijkstra": `Dijkstra settled ${query().expanded.dijkstra.toLocaleString()} states (all velocities, summed per cell).`,
    "exp-blind": `A* with the velocity-blind bound expanded ${query().expanded.blind.toLocaleString()} states.`,
    "exp-modes": `A* with the velocity-mode bound expanded ${query().expanded.modes.toLocaleString()} states.` };
  $("layer-note").textContent = notes[lSel.value];
}

// hover on the map
const mapTip = $("map-tip");
map.addEventListener("mousemove", ev => {
  const rect = map.getBoundingClientRect();
  const qx = Math.floor((ev.clientX - rect.left) / rect.width * N) + 1;
  const qy = N - Math.floor((ev.clientY - rect.top) / rect.height * N);
  if (qx < 1 || qx > N || qy < 1 || qy > N) { mapTip.hidden = true; return; }
  const i = (qy - 1) * N + (qx - 1);
  const vr = vRest[i], lb = lowerRest[pSel.value as "blind" | "modes"][i];
  const fmt = (z: number) => z === 255 ? "unreachable" : `${z} steps`;
  const ex = bytes(query().explored[(lSel.value.startsWith("exp-") ? lSel.value.slice(4) : "modes") as "dijkstra" | "blind" | "modes"])[i];
  mapTip.innerHTML = walls[i] ? `cell (${qx}, ${qy})<br>wall`
    : `cell (${qx}, ${qy})<br>V* at rest: ${fmt(vr)}<br>lower bound: ${fmt(lb)}` +
      (lSel.value.startsWith("exp-") ? `<br>expanded here: ${ex}` : "");
  mapTip.hidden = false;
  const box = $("map-box").getBoundingClientRect();
  const left = ev.clientX - box.left + 14, top = ev.clientY - box.top + 14;
  mapTip.style.left = `${Math.min(left, box.width - 170)}px`;
  mapTip.style.top = `${Math.min(top, box.height - 80)}px`;
});
map.addEventListener("mouseleave", () => { mapTip.hidden = true; });

// ---- stats ------------------------------------------------------------------------------
function drawStats(): void {
  const r = run(), q = query();
  const ratio = r.cost / q.cost;
  const stats: [string, string][] = [
    ["optimal cost V*", `${q.cost}`],
    ["controller cost V_U", `${r.cost}`],
    ["true cost / optimum", ratio.toFixed(3)],
    ["guaranteed bound", r.bound === null ? "–" : r.bound.toFixed(3)],
    ["lower bound at start", r.lower === null ? "–" : `${r.lower}`],
    ["corridor used", `w ${r.used_width} · m ${r.used_mode_width}`]];
  $("stats").innerHTML = stats.map(([k, v]) => `<div class="stat"><span class="v">${v}</span><span class="k">${k}</span></div>`).join("");
  const pill = (ok: boolean, yes: string, no: string) =>
    `<span class="pill"><span class="dot" style="background:${ok ? "var(--good)" : "var(--serious)"}"></span>${ok ? "✓ " + yes : "△ " + no}</span>`;
  $("pills").innerHTML = pill(r.certified, "Lyapunov certified: positive definite, strict decrease, invariant", "not certified") +
    pill(r.first_try, "requested corridor sufficed", `widened to w ${r.used_width}, m ${r.used_mode_width}`);
}

// ---- sandwich chart ---------------------------------------------------------------------
const SVGNS = "http://www.w3.org/2000/svg";
function el(tag: string, attrs: Record<string, string | number>, parent?: Element): SVGElement {
  const e = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
  if (parent) parent.appendChild(e);
  return e;
}
function niceTicks(max: number, count: number): number[] {
  const step = Math.max(1, Math.ceil(max / count / 5) * 5);
  const t: number[] = [];
  for (let v = 0; v <= max + 1e-9; v += step) t.push(v);
  return t;
}
function drawSandwich(): void {
  const r = run();
  const host = $("sandwich");
  host.innerHTML = "";
  const W = 560, H = 270, m = { l: 40, r: 16, t: 12, b: 34 };
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Cost-to-go along the closed loop" }, host);
  const K = r.V_U.length - 1, ymax = Math.max(...r.V_U, 1);
  const x = (k: number) => m.l + (k / Math.max(K, 1)) * (W - m.l - m.r);
  const y = (v: number) => H - m.b - (v / ymax) * (H - m.t - m.b);
  for (const t of niceTicks(ymax, 5)) {
    el("line", { x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), class: "gridline" }, svg);
    const lab = el("text", { x: m.l - 6, y: y(t) + 4, "text-anchor": "end" }, svg); lab.textContent = String(t);
  }
  for (const t of niceTicks(K, 6)) {
    const lab = el("text", { x: x(t), y: H - m.b + 16, "text-anchor": "middle" }, svg); lab.textContent = String(t);
  }
  const xl = el("text", { x: (m.l + W - m.r) / 2, y: H - 2, "text-anchor": "middle" }, svg); xl.textContent = "step k of the closed loop";
  // band between the lower bound and V_U: the sandwich
  const lower = r.lower_along.map(v => (v === null ? 0 : v));
  let band = "";
  r.V_U.forEach((v, k) => { band += `${k === 0 ? "M" : "L"}${x(k)},${y(v)}`; });
  for (let k = K; k >= 0; k--) band += `L${x(k)},${y(lower[k])}`;
  el("path", { d: band + "Z", fill: css("--s-ctl"), "fill-opacity": 0.08, stroke: "none" }, svg);
  const line = (vals: number[], color: string, dash = "") => {
    let d = "";
    vals.forEach((v, k) => { d += `${k === 0 ? "M" : "L"}${x(k)},${y(v)}`; });
    el("path", { d, fill: "none", stroke: color, "stroke-width": 2, "stroke-dasharray": dash, "stroke-linejoin": "round" }, svg);
    vals.forEach((v, k) => el("circle", { cx: x(k), cy: y(v), r: 3, fill: color, stroke: css("--surface"), "stroke-width": 1.5 }, svg));
  };
  // reference: the slowest decrease the certificate allows
  el("path", { d: `M${x(0)},${y(r.V_U[0])}L${x(Math.min(K, r.V_U[0]))},${y(Math.max(0, r.V_U[0] - Math.min(K, r.V_U[0])))}`,
    stroke: css("--ink-3"), "stroke-dasharray": "3 4", "stroke-width": 1, fill: "none" }, svg);
  line(lower, css("--s-low"));
  line(r.V_opt, css("--s-opt"));
  line(r.V_U, css("--s-ctl"));
  // direct labels at the start
  const lab = (txt: string, v: number, color: string, dy: number) => {
    const t = el("text", { x: x(0) + 8, y: y(v) + dy, fill: color }, svg); t.textContent = txt;
  };
  lab("V_U", r.V_U[0], css("--ink"), -6);
  // hover crosshair
  const tip = $("sand-tip"), cross = el("line", { y1: m.t, y2: H - m.b, stroke: css("--ink-3"), "stroke-width": 1, visibility: "hidden" }, svg);
  const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent" }, svg);
  hit.addEventListener("mousemove", (ev: Event) => {
    const e = ev as MouseEvent, rect = (svg as SVGSVGElement).getBoundingClientRect();
    const px = (e.clientX - rect.left) / rect.width * W;
    const k = Math.max(0, Math.min(K, Math.round((px - m.l) / (W - m.l - m.r) * K)));
    cross.setAttribute("x1", String(x(k))); cross.setAttribute("x2", String(x(k))); cross.setAttribute("visibility", "visible");
    const s = r.traj[k];
    tip.innerHTML = `step ${k}: q = (${s[0]}, ${s[1]}), v = (${s[2]}, ${s[3]})<br>V_U ${r.V_U[k]} · V* ${r.V_opt[k]} · lower ${lower[k]}`;
    tip.hidden = false;
    const host = tip.parentElement as HTMLElement, hb = host.getBoundingClientRect();
    tip.style.left = `${Math.min(e.clientX - hb.left + 12, hb.width - 260)}px`;
    tip.style.top = `${e.clientY - hb.top + 12}px`;
  });
  hit.addEventListener("mouseleave", () => { tip.hidden = true; cross.setAttribute("visibility", "hidden"); });
}

// ---- phase portrait ------------------------------------------------------------------
function drawPhase(): void {
  const c = $<HTMLCanvasElement>("phase"), ctx = c.getContext("2d") as CanvasRenderingContext2D;
  const P = data.phase, nv = 2 * data.vmax + 1, nq = P.nq;
  const pad = { l: 44, r: 12, t: 12, b: 40 };
  const w = c.width - pad.l - pad.r, h = c.height - pad.t - pad.b;
  const X = (q: number) => pad.l + (q - 0.5) / nq * w;
  const Y = (v: number) => pad.t + (data.vmax + 0.5 - v) / nv * h;
  ctx.clearRect(0, 0, c.width, c.height);
  const fills: Record<number, string> = { 1: css("--acc-plus"), [-1]: css("--acc-minus"), 0: css("--acc-zero") };
  for (let iv = 0; iv < nv; iv++) for (let q = 0; q < nq; q++) {
    const a = P.accel[iv * nq + q], v = iv - data.vmax;
    ctx.fillStyle = P.value[iv * nq + q] < 0 ? css("--surface") : fills[a];
    ctx.fillRect(X(q + 1) - w / nq / 2, Y(v) - h / nv / 2, w / nq + 0.5, h / nv + 0.5);
  }
  // isochrones
  const val = (x: number, y: number) => { const z = P.value[y * nq + x]; return z < 0 ? Infinity : z; };
  ctx.strokeStyle = isDark() ? "rgba(255,255,255,0.5)" : "rgba(11,11,11,0.35)"; ctx.lineWidth = 1;
  const vmaxVal = Math.max(...P.value);
  for (let level = 2; level < vmaxVal; level += 2) {
    ctx.beginPath();
    for (const [x0, y0, x1, y1] of contours(val, nq, nv, level)) {
      ctx.moveTo(X(x0 + 1), Y(y0 - data.vmax)); ctx.lineTo(X(x1 + 1), Y(y1 - data.vmax));
    }
    ctx.stroke();
  }
  // axes
  ctx.fillStyle = css("--ink-2"); ctx.font = `22px ${css("--font-mono")}`; ctx.textAlign = "center";
  for (let q = 0; q <= nq; q += 16) ctx.fillText(String(q - P.goal), X(Math.max(q, 1)), c.height - 12);
  ctx.textAlign = "right";
  for (let v = -data.vmax; v <= data.vmax; v += 3) ctx.fillText(String(v), pad.l - 8, Y(v) + 7);
  ctx.textAlign = "left";
  ctx.fillText("position − goal →", pad.l + 6, c.height - 12 - 0);
  // trajectories
  ctx.strokeStyle = css("--ink"); ctx.fillStyle = css("--ink"); ctx.lineWidth = 2.5;
  for (const tr of P.trajectories) {
    ctx.beginPath();
    tr.forEach(([q, v], k) => k === 0 ? ctx.moveTo(X(q), Y(v)) : ctx.lineTo(X(q), Y(v)));
    ctx.stroke();
    for (const [q, v] of tr) { ctx.beginPath(); ctx.arc(X(q), Y(v), 4, 0, 2 * Math.PI); ctx.fill(); }
  }
  ctx.fillStyle = css("--good");
  ctx.beginPath(); ctx.arc(X(P.goal), Y(0), 8, 0, 2 * Math.PI); ctx.fill();
  ctx.save(); ctx.translate(14, pad.t + h / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillStyle = css("--ink-2"); ctx.textAlign = "center"; ctx.fillText("velocity", 0, 0); ctx.restore();

  const tip = $("phase-tip");
  c.onmousemove = ev => {
    const rect = c.getBoundingClientRect();
    const px = (ev.clientX - rect.left) / rect.width * c.width, py = (ev.clientY - rect.top) / rect.height * c.height;
    const q = Math.round((px - pad.l) / w * nq + 0.5), v = Math.round(data.vmax + 0.5 - (py - pad.t) / h * nv - 0.5);
    if (q < 1 || q > nq || v < -data.vmax || v > data.vmax) { tip.hidden = true; return; }
    const i = (v + data.vmax) * nq + (q - 1), a = P.accel[i];
    tip.innerHTML = `q − goal = ${q - P.goal}, v = ${v}<br>minimum time ${P.value[i] < 0 ? "–" : P.value[i]} · a* = ${a > 0 ? "+1" : a < 0 ? "−1" : "0"}`;
    tip.hidden = false;
    const hb = (tip.parentElement as HTMLElement).getBoundingClientRect();
    tip.style.left = `${Math.min(ev.clientX - hb.left + 12, hb.width - 200)}px`;
    tip.style.top = `${ev.clientY - hb.top + 12}px`;
  };
  c.onmouseleave = () => { tip.hidden = true; };
}

// ---- trade-off charts ---------------------------------------------------------------------
const PLANNERS = [
  { key: "velocity-blind", label: "velocity-blind", color: "--s-opt" },
  { key: "8 headings × 3 speeds", label: "8 headings × 3 speeds", color: "--s-ctl" },
  { key: "16 headings × 3 speeds", label: "16 headings × 3 speeds", color: "--s-low" }];
function drawTradeoffs(): void {
  // first-try success by corridor setting, from the layered study
  // One mode per stalk makes the mode width moot: its w=4, 8, 16 runs line up with m=1.
  const settings = ["w=0", "w=1", "w=2", "w=2, m=1", "w=4, m=1", "w=8, m=1", "w=16, m=1"];
  const series = PLANNERS.map(p => {
    const pts: { s: string; first: number; ratio: number; worst: number }[] = [];
    for (const row of data.layered) {
      const name = String(row["method"]);
      if (!name.startsWith(p.key + ":")) continue;
      let s = name.slice(p.key.length + 1).trim().replace(", m=0", "");
      if (p.key === "velocity-blind" && ["w=4", "w=8", "w=16"].includes(s)) s += ", m=1";
      pts.push({ s, first: +row["first_try"], ratio: +row["mean_cost_over_optimal"], worst: +row["max_cost_over_optimal"] });
    }
    return { ...p, pts: pts.filter(q => settings.includes(q.s)).sort((a, b) => settings.indexOf(a.s) - settings.indexOf(b.s)) };
  });
  const used = settings.filter(s => series.some(se => se.pts.some(p => p.s === s)));
  const host = $("firsttry"); host.innerHTML = "";
  const W = 420, H = 250, m = { l: 40, r: 10, t: 10, b: 52 };
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "First-try success by corridor" }, host);
  const x = (s: string) => m.l + (used.indexOf(s) + 0.5) / used.length * (W - m.l - m.r);
  const y = (f: number) => H - m.b - f * (H - m.t - m.b);
  for (const t of [0, 0.25, 0.5, 0.75, 1]) {
    el("line", { x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), class: "gridline" }, svg);
    const l = el("text", { x: m.l - 6, y: y(t) + 4, "text-anchor": "end" }, svg); l.textContent = `${Math.round(t * 100)}%`;
  }
  used.forEach(s => {
    const l = el("text", { x: x(s), y: H - m.b + 14, "text-anchor": "end", transform: `rotate(-35 ${x(s)} ${H - m.b + 14})` }, svg);
    l.textContent = s;
  });
  const tip = $("trade-tip");
  for (const se of series) {
    const col = css(se.color);
    let d = "";
    se.pts.forEach((p, k) => { d += `${k === 0 ? "M" : "L"}${x(p.s)},${y(p.first)}`; });
    el("path", { d, fill: "none", stroke: col, "stroke-width": 2 }, svg);
    for (const p of se.pts) {
      const dot = el("circle", { cx: x(p.s), cy: y(p.first), r: 4.5, fill: col, stroke: css("--surface"), "stroke-width": 2 }, svg);
      dot.addEventListener("mouseenter", (ev: Event) => {
        const e = ev as MouseEvent, hb = (tip.parentElement as HTMLElement).getBoundingClientRect();
        tip.innerHTML = `${se.label}, ${p.s}<br>first try ${(p.first * 100).toFixed(0)}% · cost/opt ${p.ratio.toFixed(3)} (worst ${p.worst.toFixed(3)})`;
        tip.hidden = false; tip.style.left = `${Math.min(e.clientX - hb.left + 12, hb.width - 300)}px`; tip.style.top = `${e.clientY - hb.top + 12}px`;
      });
      dot.addEventListener("mouseleave", () => { tip.hidden = true; });
    }
  }

  // guaranteed bound vs true ratio, from the certified study
  const host2 = $("bounds"); host2.innerHTML = "";
  const svg2 = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Guaranteed bound against true ratio" }, host2);
  const rows = data.certified;
  const groups = PLANNERS.map(p => ({ ...p, rows: rows.filter(r => r["planner"] === p.key) }));
  const ymax = Math.max(...rows.map(r => +r["max_bound"]));
  const yb = (v: number) => H - m.b - (v - 1) / (ymax - 1) * (H - m.t - m.b);
  for (let t = 1; t <= ymax + 1e-9; t += ymax > 4 ? 1 : 0.5) {
    el("line", { x1: m.l, x2: W - m.r, y1: yb(t), y2: yb(t), class: "gridline" }, svg2);
    const l = el("text", { x: m.l - 6, y: yb(t) + 4, "text-anchor": "end" }, svg2); l.textContent = `${t}×`;
  }
  const slots = groups.flatMap(g => g.rows.map(r => ({ g, r })));
  const xs = (i: number) => m.l + (i + 0.5) / slots.length * (W - m.l - m.r);
  slots.forEach(({ g, r }, i) => {
    const col = css(g.color);
    el("line", { x1: xs(i), x2: xs(i), y1: yb(+r["mean_bound"]), y2: yb(+r["max_bound"]), stroke: col, "stroke-width": 2, "stroke-opacity": 0.5 }, svg2);
    el("rect", { x: xs(i) - 6, y: yb(+r["mean_bound"]) - 2, width: 12, height: 4, rx: 1, fill: col }, svg2);
    el("circle", { cx: xs(i), cy: yb(+r["mean_cost_over_optimal"]), r: 4.5, fill: css("--surface"), stroke: col, "stroke-width": 2 }, svg2);
    const l = el("text", { x: xs(i), y: H - m.b + 14, "text-anchor": "end", transform: `rotate(-35 ${xs(i)} ${H - m.b + 14})` }, svg2);
    l.textContent = `w=${r["width"]}`;
    const hit = el("rect", { x: xs(i) - 14, y: m.t, width: 28, height: H - m.t - m.b, fill: "transparent" }, svg2);
    hit.addEventListener("mouseenter", (ev: Event) => {
      const e = ev as MouseEvent, hb = (tip.parentElement as HTMLElement).getBoundingClientRect();
      tip.innerHTML = `${g.label}, w=${r["width"]} m=${r["mode_width"]}<br>true ${(+r["mean_cost_over_optimal"]).toFixed(3)} · bound mean ${(+r["mean_bound"]).toFixed(2)}, worst ${(+r["max_bound"]).toFixed(2)}`;
      tip.hidden = false; tip.style.left = `${Math.min(e.clientX - hb.left + 12, hb.width - 300)}px`; tip.style.top = `${e.clientY - hb.top + 12}px`;
    });
    hit.addEventListener("mouseleave", () => { tip.hidden = true; });
  });

  $("trade-legend").innerHTML = PLANNERS.map(p => `<span class="key"><span class="swatch" style="background:var(${p.color})"></span>${p.label}</span>`).join("") +
    `<span class="key">right: bar = mean bound, whisker = worst, ring = true ratio</span>`;
  const tr = data.certified.map(r => `<tr><td>${r["planner"]}</td><td>${r["width"]}</td><td>${(+r["first_try"] * 100).toFixed(0)}%</td>` +
    `<td>${(+r["mean_cost_over_optimal"]).toFixed(3)}</td><td>${(+r["mean_bound"]).toFixed(3)}</td><td>${(+r["max_bound"]).toFixed(3)}</td>` +
    `<td>${(+r["seconds_per_query"]).toFixed(3)} s</td></tr>`).join("");
  $("trade-table").innerHTML = `<table><thead><tr><th>planner</th><th>width</th><th>first try</th><th>cost/opt</th><th>bound mean</th><th>bound worst</th><th>time/query</th></tr></thead><tbody>${tr}</tbody></table>`;
}

// ---- wiring -----------------------------------------------------------------------------
function render(): void { drawMap(); drawStats(); drawSandwich(); }
for (const s of [qSel, pSel, cSel, lSel]) s.addEventListener("change", render);
for (const t of toggles) t.addEventListener("change", drawMap);
const redrawAll = () => { render(); drawPhase(); drawTradeoffs(); };
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redrawAll);
new MutationObserver(redrawAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
redrawAll();
