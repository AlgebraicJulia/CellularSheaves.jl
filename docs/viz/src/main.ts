// The certified corridor planner demo: renders the data written by
// docs/scripts/tropical_viz_export.jl. Built with esbuild (docs/viz/build.sh).
// No solver runs here: every value comes from the Julia export; the page only
// recomputes geometry (corridors from routes) and the stalkwise maps P and E.

type State = [number, number, number, number];                  // qx, qy, vx, vy
type RoutePoint = [number, number, number, number, number];     // qx, qy, heading, speed band (-1: none), mode

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
  stalk: { modes16: number[]; K16: number; neighbors16: number[][]; modes_cert: number[]; Kcert: number;
           certified_fraction: number; Vup_rest: string };
}

const data: Data = JSON.parse((document.getElementById("demo-data") as HTMLScriptElement).textContent || "{}");
const N = data.n, VMAX = data.vmax, NV = 2 * VMAX + 1, K16 = data.stalk.K16;
const INF = 255;

// ---- decoding --------------------------------------------------------------------
function bytes(b64: string): Uint8Array {
  const s = atob(b64.trim());
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
async function blob(id: string): Promise<Uint8Array> {
  const el = document.getElementById(id);
  if (!el) throw new Error(`missing ${id}`);
  const gz = bytes(el.textContent || "");
  const stream = new Blob([gz]).stream().pipeThrough(new DecompressionStream("gzip"));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}
const walls = bits(data.walls, N * N);
const vRest = bytes(data.V_rest);
const lowerRest = { blind: bytes(data.lower_rest.blind), modes: bytes(data.lower_rest.modes) };
const vupRest = bytes(data.stalk.Vup_rest);
const exploredCache = new Map<string, Uint8Array>();
function explored(q: Query, which: "dijkstra" | "blind" | "modes"): Uint8Array {
  const key = `${q.id}:${which}`;
  let v = exploredCache.get(key);
  if (!v) { v = bytes(q.explored[which]); exploredCache.set(key, v); }
  return v;
}
// Full state fields, filled once the gzipped blobs are decompressed.
let vstar: Uint8Array | null = null, vup: Uint8Array | null = null, lower16: Uint8Array | null = null;
const stateIndex = (qx: number, qy: number, vx: number, vy: number) =>
  (qx - 1) + N * (qy - 1) + N * N * (vx + VMAX) + N * N * NV * (vy + VMAX);
const modeOf = (modes: number[], vx: number, vy: number) => modes[(vy + VMAX) * NV + (vx + VMAX)];

// ---- theme ----------------------------------------------------------------------------
function css(name: string): string { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
function isDark(): boolean {
  const t = document.documentElement.getAttribute("data-theme");
  if (t) return t === "dark";
  return window.matchMedia("(prefers-color-scheme: dark)").matches;
}
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
const rgb = (c: [number, number, number]) => `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`;

// ---- marching squares for isochrones ----------------------------------------------
type Seg = [number, number, number, number];
const CASES: Record<number, number[][]> = {
  1: [[0, 3]], 2: [[0, 1]], 3: [[1, 3]], 4: [[1, 2]], 5: [[0, 1], [2, 3]], 6: [[0, 2]], 7: [[2, 3]],
  8: [[2, 3]], 9: [[0, 2]], 10: [[0, 3], [1, 2]], 11: [[1, 2]], 12: [[1, 3]], 13: [[0, 1]], 14: [[0, 3]] };
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
    for (const [a, b] of CASES[idx]) segs.push([e[a][0], e[a][1], e[b][0], e[b][1]]);
  }
  return segs;
}

// ---- corridors, recomputed from routes exactly as TropicalAbstraction.corridor ----------
const roundHalfEven = (x: number) => {                             // Julia's round(Int, x)
  const f = Math.floor(x), d = x - f;
  return d > 0.5 ? f + 1 : d < 0.5 ? f : (f % 2 === 0 ? f : f + 1);
};
function segmentCells(a: [number, number], b: [number, number]): [number, number][] {
  const steps = Math.max(Math.abs(b[0] - a[0]), Math.abs(b[1] - a[1]));
  const pts: [number, number][] = [a];
  for (let s = 1; s <= steps; s++)
    pts.push([roundHalfEven(a[0] + (b[0] - a[0]) * (s / steps)), roundHalfEven(a[1] + (b[1] - a[1]) * (s / steps))]);
  return pts;
}
function dilateModes(modes: number[], width: number): Set<number> {
  let S = new Set(modes);
  for (let k = 0; k < width; k++) {
    const T = new Set(S);
    for (const m of S) for (const nb of data.stalk.neighbors16[m - 1]) T.add(nb);
    S = T;
  }
  return S;
}
interface Corridor { cells: Uint8Array; modes: Map<number, Set<number>> | null }   // modes: null = every velocity
const corridorCache = new Map<string, Corridor>();
function corridorOf(q: Query, r: Run): Corridor {
  const key = `${q.id}:${r.planner}:${r.width}`;
  const hit = corridorCache.get(key);
  if (hit) return hit;
  const cells = new Uint8Array(N * N);
  const modes = r.planner === "modes" ? new Map<number, Set<number>>() : null;
  const pieces: [RoutePoint, RoutePoint][] = r.route.length === 1 ? [[r.route[0], r.route[0]]]
    : r.route.slice(0, -1).map((p, i) => [p, r.route[i + 1]] as [RoutePoint, RoutePoint]);
  const w = r.used_width;
  for (const [a, b] of pieces) {
    const ms = modes ? dilateModes([a[4], b[4]], r.used_mode_width) : null;
    for (const [x, y] of segmentCells([a[0], a[1]], [b[0], b[1]]))
      for (let yy = Math.max(1, y - w); yy <= Math.min(N, y + w); yy++)
        for (let xx = Math.max(1, x - w); xx <= Math.min(N, x + w); xx++) {
          const i = (yy - 1) * N + (xx - 1);
          cells[i] = 1;
          if (modes && ms) {
            const s = modes.get(i);
            if (s) for (const m of ms) s.add(m); else modes.set(i, new Set(ms));
          }
        }
  }
  // Check against the cells the export projected from the planner's own corridor.
  const exported = bits(r.corridor, N * N);
  let mismatch = 0;
  for (let i = 0; i < N * N; i++) if (exported[i] !== cells[i]) mismatch++;
  if (mismatch > 0) console.warn(`corridor of query ${q.id} (${r.planner}, w=${r.width}): ${mismatch} cells differ from the export`);
  const c = { cells: mismatch > 0 ? exported : cells, modes };
  corridorCache.set(key, c);
  return c;
}

// ---- state -----------------------------------------------------------------------------
const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const qSel = $<HTMLSelectElement>("q"), pSel = $<HTMLSelectElement>("planner"), cSel = $<HTMLSelectElement>("corr");
const lSel = $<HTMLSelectElement>("layer"), fSel = $<HTMLSelectElement>("field");
const T_CORR = $<HTMLInputElement>("t-corr"), T_ROUTE = $<HTMLInputElement>("t-route"), T_CTL = $<HTMLInputElement>("t-ctl");
const T_OPT = $<HTMLInputElement>("t-opt"), T_ISO = $<HTMLInputElement>("t-iso"), T_CMP = $<HTMLInputElement>("t-cmp");
const stepIn = $<HTMLInputElement>("step"), playBtn = $<HTMLButtonElement>("play"), followBtn = $<HTMLButtonElement>("follow");
data.queries.forEach((q, i) => {
  const o = document.createElement("option");
  const [x, y, vx, vy] = q.src;
  o.value = String(i);
  o.textContent = `#${q.id}  q = (${x}, ${y}), v = (${vx}, ${vy}) · optimum ${q.cost} steps`;
  qSel.appendChild(o);
});
const query = (): Query => data.queries[+qSel.value];
const runFor = (planner: string): Run => query().runs.find(r => r.planner === planner && r.width === +cSel.value) as Run;
const run = (): Run => runFor(pSel.value);
const other = (): Run => runFor(pSel.value === "modes" ? "blind" : "modes");
let step = 0;
let pinned: [number, number] | null = null;              // a clicked cell; null: follow the run
function currentCell(): [number, number] {
  if (pinned) return pinned;
  const s = run().traj[Math.min(step, run().traj.length - 1)];
  return [s[0], s[1]];
}

// ---- workspace map ------------------------------------------------------------------
const map = $<HTMLCanvasElement>("map");
const mctx = map.getContext("2d") as CanvasRenderingContext2D;
const S = map.width / N;
const cx = (qx: number) => (qx - 0.5) * S;
const cy = (qy: number) => (N - qy + 0.5) * S;

function layerValues(): { v: Float64Array; max: number; kind: "time" | "count" } {
  const layer = lSel.value;
  const v = new Float64Array(N * N);
  let max = 0;
  if (layer === "value" || layer === "lower" || layer === "cert") {
    const src = layer === "value" ? vRest : layer === "cert" ? vupRest : lowerRest[pSel.value as "blind" | "modes"];
    for (let i = 0; i < N * N; i++) { v[i] = src[i] === INF ? Infinity : src[i]; if (isFinite(v[i])) max = Math.max(max, v[i]); }
    return { v, max, kind: "time" };
  }
  const src = explored(query(), layer.slice(4) as "dijkstra" | "blind" | "modes");
  for (let i = 0; i < N * N; i++) { v[i] = src[i]; max = Math.max(max, v[i]); }
  return { v, max, kind: "count" };
}

function outline(cells: Uint8Array, color: string, dash: number[]): void {
  mctx.strokeStyle = color; mctx.lineWidth = 1.2; mctx.setLineDash(dash);
  mctx.beginPath();
  const inC = (x: number, y: number) => x >= 1 && x <= N && y >= 1 && y <= N && cells[(y - 1) * N + (x - 1)] === 1;
  for (let i = 0; i < N * N; i++) if (cells[i]) {
    const qx = (i % N) + 1, qy = Math.floor(i / N) + 1, x0 = (qx - 1) * S, y0 = (N - qy) * S;
    if (!inC(qx - 1, qy)) { mctx.moveTo(x0, y0); mctx.lineTo(x0, y0 + S); }
    if (!inC(qx + 1, qy)) { mctx.moveTo(x0 + S, y0); mctx.lineTo(x0 + S, y0 + S); }
    if (!inC(qx, qy + 1)) { mctx.moveTo(x0, y0); mctx.lineTo(x0 + S, y0); }
    if (!inC(qx, qy - 1)) { mctx.moveTo(x0, y0 + S); mctx.lineTo(x0 + S, y0 + S); }
  }
  mctx.stroke(); mctx.setLineDash([]);
}
function polyline(states: [number, number][], color: string, width: number, dash: number[] = [], dots = true): void {
  mctx.strokeStyle = color; mctx.lineWidth = width; mctx.lineJoin = "round"; mctx.setLineDash(dash);
  mctx.beginPath();
  states.forEach(([x, y], k) => k === 0 ? mctx.moveTo(cx(x), cy(y)) : mctx.lineTo(cx(x), cy(y)));
  mctx.stroke(); mctx.setLineDash([]);
  if (dots) {
    mctx.fillStyle = color;
    for (const [x, y] of states) { mctx.beginPath(); mctx.arc(cx(x), cy(y), 3.2, 0, 2 * Math.PI); mctx.fill(); }
  }
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
  if (lSel.value === "cert") {                              // hatch what the certificate does not cover
    mctx.strokeStyle = dark ? "rgba(255,255,255,0.18)" : "rgba(11,11,11,0.14)"; mctx.lineWidth = 1;
    mctx.beginPath();
    for (let i = 0; i < N * N; i++) if (!walls[i] && !isFinite(v[i])) {
      const qx = (i % N) + 1, qy = Math.floor(i / N) + 1, x0 = (qx - 1) * S, y0 = (N - qy) * S;
      mctx.moveTo(x0, y0 + S); mctx.lineTo(x0 + S, y0);
    }
    mctx.stroke();
  }

  const r = run(), q = query();
  if (T_CORR.checked) {
    const cor = corridorOf(q, r);
    mctx.fillStyle = css("--corridor");
    for (let i = 0; i < N * N; i++) if (cor.cells[i] && !walls[i]) {
      const qx = (i % N) + 1, qy = Math.floor(i / N) + 1;
      mctx.fillRect((qx - 1) * S, (N - qy) * S, S, S);
    }
    outline(cor.cells, css("--s-low"), []);
  }
  if (T_ISO.checked && kind === "time") {
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
  mctx.fillStyle = css("--good");
  for (const [gx, gy] of data.goal) mctx.fillRect((gx - 1) * S, (N - gy) * S, S, S);
  if (T_CMP.checked) {                                    // the other planner, dashed
    const o = other();
    outline(corridorOf(q, o).cells, css("--ink-2"), [4, 3]);
    polyline(o.traj.map(s => [s[0], s[1]] as [number, number]), css("--ink-2"), 2, [5, 4], false);
  }
  if (T_ROUTE.checked) {
    polyline(r.route.map(p => [p[0], p[1]] as [number, number]), css("--s-route"), 2, [6, 4], true);
    for (const [x, y, ang, band] of r.route) if (band > 0) {
      const L = 4 + 3.5 * band;
      mctx.strokeStyle = css("--s-route"); mctx.lineWidth = 1.5;
      mctx.beginPath(); mctx.moveTo(cx(x), cy(y)); mctx.lineTo(cx(x) + L * Math.cos(ang), cy(y) - L * Math.sin(ang)); mctx.stroke();
    }
  }
  if (T_OPT.checked) polyline(q.optimal.map(s => [s[0], s[1]] as [number, number]), css("--s-opt"), 2);
  if (T_CTL.checked) polyline(r.traj.map(s => [s[0], s[1]] as [number, number]), css("--s-ctl"), 2.5);
  // the vehicle at the current step, with its velocity
  const k = Math.min(step, r.traj.length - 1), [x, y, vx, vy] = r.traj[k];
  mctx.fillStyle = css("--surface"); mctx.strokeStyle = css("--ink"); mctx.lineWidth = 2.5;
  mctx.beginPath(); mctx.arc(cx(x), cy(y), 7, 0, 2 * Math.PI); mctx.fill(); mctx.stroke();
  if (vx !== 0 || vy !== 0) {
    const tx = cx(x + vx), ty = cy(y + vy), ang = Math.atan2(ty - cy(y), tx - cx(x));
    mctx.beginPath(); mctx.moveTo(cx(x), cy(y)); mctx.lineTo(tx, ty); mctx.stroke();
    mctx.beginPath(); mctx.moveTo(tx, ty);
    mctx.lineTo(tx - 9 * Math.cos(ang - 0.4), ty - 9 * Math.sin(ang - 0.4));
    mctx.lineTo(tx - 9 * Math.cos(ang + 0.4), ty - 9 * Math.sin(ang + 0.4)); mctx.closePath();
    mctx.fillStyle = css("--ink"); mctx.fill();
  }
  // the inspected cell
  const [ix, iy] = currentCell();
  mctx.strokeStyle = css("--ink"); mctx.lineWidth = 1.5; mctx.setLineDash(pinned ? [] : [2, 2]);
  mctx.strokeRect((ix - 1) * S - 2, (N - iy) * S - 2, S + 4, S + 4); mctx.setLineDash([]);

  const leg = $("map-legend");
  const rampCss = (colors: string[]) => `linear-gradient(90deg, ${(dark ? [...colors].reverse() : colors).join(",")})`;
  const scale = kind === "time"
    ? `<span class="key"><span class="ramp" style="background:${rampCss(BLUE)}"></span><span class="num">0–${max} steps</span></span>`
    : `<span class="key"><span class="ramp" style="background:${rampCss(ORANGE)}"></span><span class="num">0–${max} states per cell</span></span>`;
  leg.innerHTML = scale +
    `<span class="key"><span class="swatch" style="background:var(--s-ctl)"></span>certified controller</span>` +
    `<span class="key"><span class="swatch" style="background:var(--s-opt)"></span>optimal</span>` +
    `<span class="key"><span class="swatch" style="background:var(--s-route)"></span>route (ticks: mode heading)</span>` +
    `<span class="key"><span class="swatch box" style="background:var(--corridor);outline:1px solid var(--s-low)"></span>corridor</span>` +
    (T_CMP.checked ? `<span class="key"><span class="swatch" style="background:var(--ink-2)"></span>other planner (dashed)</span>` : "") +
    `<span class="key"><span class="swatch box" style="background:var(--good)"></span>goal, stopped</span>` +
    `<span class="key"><span class="swatch box" style="background:var(--wall)"></span>wall</span>`;
  const notes: Record<string, string> = {
    value: "Minimum time to the goal starting at rest; contours every 5 steps.",
    lower: "Coarse value pulled back along the oplax refinement: never above V*.",
    cert: `Starts at rest certified by the ∀∃ coarse Lyapunov function (8 headings × 6 speeds); hatched: not certified. Over all velocities it certifies ${(data.stalk.certified_fraction * 100).toFixed(0)}% of the states that can reach the goal.`,
    "exp-dijkstra": `Dijkstra settled ${q.expanded.dijkstra.toLocaleString()} states (all velocities, summed per cell).`,
    "exp-blind": `A* with the velocity-blind bound expanded ${q.expanded.blind.toLocaleString()} states.`,
    "exp-modes": `A* with the velocity-mode bound expanded ${q.expanded.modes.toLocaleString()} states.` };
  $("layer-note").textContent = notes[lSel.value];
}

const mapTip = $("map-tip");
function mapCell(ev: MouseEvent): [number, number] | null {
  const rect = map.getBoundingClientRect();
  const qx = Math.floor((ev.clientX - rect.left) / rect.width * N) + 1;
  const qy = N - Math.floor((ev.clientY - rect.top) / rect.height * N);
  return qx < 1 || qx > N || qy < 1 || qy > N ? null : [qx, qy];
}
map.addEventListener("mousemove", ev => {
  const c = mapCell(ev);
  if (!c) { mapTip.hidden = true; return; }
  const [qx, qy] = c, i = (qy - 1) * N + (qx - 1);
  const fmt = (z: number) => z === INF ? "unreachable" : `${z} steps`;
  const exp = lSel.value.startsWith("exp-") ? `<br>expanded here: ${explored(query(), lSel.value.slice(4) as "dijkstra" | "blind" | "modes")[i]}` : "";
  mapTip.innerHTML = walls[i] ? `cell (${qx}, ${qy})<br>wall`
    : `cell (${qx}, ${qy}) · click to inspect<br>V* at rest: ${fmt(vRest[i])}<br>lower bound: ${fmt(lowerRest[pSel.value as "blind" | "modes"][i])}` +
      `<br>∀∃ certificate: ${vupRest[i] === INF ? "not certified" : vupRest[i] + " steps"}` + exp;
  mapTip.hidden = false;
  const box = $("map-box").getBoundingClientRect();
  mapTip.style.left = `${Math.min(ev.clientX - box.left + 14, box.width - 190)}px`;
  mapTip.style.top = `${Math.min(ev.clientY - box.top + 14, box.height - 96)}px`;
});
map.addEventListener("mouseleave", () => { mapTip.hidden = true; });
map.addEventListener("click", ev => {
  const c = mapCell(ev);
  if (!c || walls[(c[1] - 1) * N + (c[0] - 1)]) return;
  pinned = c;
  followBtn.setAttribute("aria-pressed", "false");
  drawMap(); drawStalk();
});
followBtn.addEventListener("click", () => { pinned = null; followBtn.setAttribute("aria-pressed", "true"); drawMap(); drawStalk(); });

// ---- stats ------------------------------------------------------------------------------
function drawStats(): void {
  const r = run(), q = query();
  const stats: [string, string][] = [
    ["optimal cost V*", `${q.cost}`],
    ["controller cost V_U", `${r.cost}`],
    ["true cost / optimum", (r.cost / q.cost).toFixed(3)],
    ["guaranteed bound", r.bound === null ? "–" : r.bound.toFixed(3)],
    ["lower bound at start", r.lower === null ? "–" : `${r.lower}`],
    ["corridor used", `w ${r.used_width} · m ${r.used_mode_width}`]];
  $("stats").innerHTML = stats.map(([k, v]) => `<div class="stat"><span class="v">${v}</span><span class="k">${k}</span></div>`).join("");
  const pill = (ok: boolean, yes: string, no: string) =>
    `<span class="pill"><span class="dot" style="background:${ok ? "var(--good)" : "var(--serious)"}"></span>${ok ? "✓ " + yes : "△ " + no}</span>`;
  $("pills").innerHTML = pill(r.certified, "Lyapunov certified: positive definite, strict decrease, invariant", "not certified") +
    pill(r.first_try, "requested corridor sufficed", `widened to w ${r.used_width}, m ${r.used_mode_width}`);
}

// ---- SVG helpers ------------------------------------------------------------------------
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
function placeTip(tip: HTMLElement, e: MouseEvent, width: number): void {
  const hb = (tip.parentElement as HTMLElement).getBoundingClientRect();
  tip.hidden = false;
  tip.style.left = `${Math.max(0, Math.min(e.clientX - hb.left + 12, hb.width - width))}px`;
  tip.style.top = `${e.clientY - hb.top + 12}px`;
}

// ---- sandwich chart ---------------------------------------------------------------------
function drawSandwich(): void {
  const r = run();
  const host = $("sandwich");
  host.innerHTML = "";
  const W = 560, H = 270, m = { l: 40, r: 16, t: 12, b: 34 };
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Cost-to-go along the closed loop" }, host);
  const o = T_CMP.checked ? other() : null;
  const K = Math.max(r.V_U.length - 1, o ? o.V_U.length - 1 : 0, 1);
  const ymax = Math.max(...r.V_U, ...(o ? o.V_U : []), 1);
  const x = (k: number) => m.l + (k / K) * (W - m.l - m.r);
  const y = (v: number) => H - m.b - (v / ymax) * (H - m.t - m.b);
  for (const t of niceTicks(ymax, 5)) {
    el("line", { x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), class: "gridline" }, svg);
    el("text", { x: m.l - 6, y: y(t) + 4, "text-anchor": "end" }, svg).textContent = String(t);
  }
  for (const t of niceTicks(K, 6)) el("text", { x: x(t), y: H - m.b + 16, "text-anchor": "middle" }, svg).textContent = String(t);
  el("text", { x: (m.l + W - m.r) / 2, y: H - 2, "text-anchor": "middle" }, svg).textContent = "step k of the closed loop";
  const lower = r.lower_along.map(v => (v === null ? 0 : v));
  const KR = r.V_U.length - 1;
  let band = "";
  r.V_U.forEach((v, k) => { band += `${k === 0 ? "M" : "L"}${x(k)},${y(v)}`; });
  for (let k = KR; k >= 0; k--) band += `L${x(k)},${y(lower[k])}`;
  el("path", { d: band + "Z", fill: css("--s-ctl"), "fill-opacity": 0.08, stroke: "none" }, svg);
  const line = (vals: number[], color: string, dash = "", dots = true) => {
    let d = "";
    vals.forEach((v, k) => { d += `${k === 0 ? "M" : "L"}${x(k)},${y(v)}`; });
    el("path", { d, fill: "none", stroke: color, "stroke-width": 2, "stroke-dasharray": dash, "stroke-linejoin": "round" }, svg);
    if (dots) vals.forEach((v, k) => el("circle", { cx: x(k), cy: y(v), r: 3, fill: color, stroke: css("--surface"), "stroke-width": 1.5 }, svg));
  };
  el("path", { d: `M${x(0)},${y(r.V_U[0])}L${x(Math.min(KR, r.V_U[0]))},${y(Math.max(0, r.V_U[0] - Math.min(KR, r.V_U[0])))}`,
    stroke: css("--ink-3"), "stroke-dasharray": "3 4", "stroke-width": 1, fill: "none" }, svg);
  if (o) line(o.V_U, css("--ink-2"), "5 4", false);
  line(lower, css("--s-low"));
  line(r.V_opt, css("--s-opt"));
  line(r.V_U, css("--s-ctl"));
  el("text", { x: x(0) + 8, y: y(r.V_U[0]) - 6, fill: css("--ink") }, svg).textContent = "V_U";
  // the current step
  const k0 = Math.min(step, KR);
  el("line", { x1: x(k0), x2: x(k0), y1: m.t, y2: H - m.b, stroke: css("--ink"), "stroke-width": 1.2 }, svg);
  const tip = $("sand-tip"), cross = el("line", { y1: m.t, y2: H - m.b, stroke: css("--ink-3"), "stroke-width": 1, visibility: "hidden" }, svg);
  const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent", style: "cursor:pointer" }, svg);
  const kAt = (e: MouseEvent) => {
    const rect = (svg as SVGSVGElement).getBoundingClientRect();
    const px = (e.clientX - rect.left) / rect.width * W;
    return Math.max(0, Math.min(KR, Math.round((px - m.l) / (W - m.l - m.r) * K)));
  };
  hit.addEventListener("mousemove", (ev: Event) => {
    const e = ev as MouseEvent, k = kAt(e);
    cross.setAttribute("x1", String(x(k))); cross.setAttribute("x2", String(x(k))); cross.setAttribute("visibility", "visible");
    const s = r.traj[k];
    tip.innerHTML = `step ${k}: q = (${s[0]}, ${s[1]}), v = (${s[2]}, ${s[3]})<br>V_U ${r.V_U[k]} · V* ${r.V_opt[k]} · lower ${lower[k]}` +
      (o && k < o.V_U.length ? `<br>other planner V_U ${o.V_U[k]}` : "") + `<br>click to go to this step`;
    placeTip(tip, e, 280);
  });
  hit.addEventListener("mouseleave", () => { tip.hidden = true; cross.setAttribute("visibility", "hidden"); });
  hit.addEventListener("click", (ev: Event) => { setStep(kAt(ev as MouseEvent)); });
}

// ---- the velocity stalk -------------------------------------------------------------------
const stalkC = $<HTMLCanvasElement>("stalk");
const sctx = stalkC.getContext("2d") as CanvasRenderingContext2D;
function stalkValue(field: string, qx: number, qy: number, vx: number, vy: number): number {
  if (!vstar || !vup || !lower16) return NaN;
  const i = stateIndex(qx, qy, vx, vy);
  if (field === "vstar") return vstar[i] === INF ? Infinity : vstar[i];
  if (field === "cert") return vup[i] === INF ? Infinity : vup[i];
  if (field === "lower") {
    if (pSel.value === "blind") { const z = lowerRest.blind[(qy - 1) * N + (qx - 1)]; return z === INF ? Infinity : z; }
    const z = lower16[((qy - 1) * N + (qx - 1)) * K16 + modeOf(data.stalk.modes16, vx, vy) - 1];
    return z === INF ? Infinity : z;
  }
  // E·P·V*: the minimum of V* over the velocity's mode in this cell (the planner's partition)
  const modes = pSel.value === "blind" ? null : data.stalk.modes16, mine = modes ? modeOf(modes, vx, vy) : 1;
  let best = Infinity;
  for (let wy = -VMAX; wy <= VMAX; wy++) for (let wx = -VMAX; wx <= VMAX; wx++)
    if (!modes || modeOf(modes, wx, wy) === mine) {
      const z = vstar[stateIndex(qx, qy, wx, wy)];
      if (z !== INF) best = Math.min(best, z);
    }
  return best;
}
function drawStalk(): void {
  const [qx, qy] = currentCell();
  const field = fSel.value, dark = isDark();
  $("stalk-title").textContent = `Velocity stalk at cell (${qx}, ${qy})` + (pinned ? "" : " · following the run");
  const W = stalkC.width, pad = 40, cs = (W - 2 * pad) / NV;
  sctx.clearRect(0, 0, W, W);
  if (!vstar) { sctx.fillStyle = css("--ink-2"); sctx.font = `18px ${css("--font-body")}`; sctx.fillText("Decompressing the stalk data…", pad, W / 2); return; }
  const vals: number[][] = [];
  let vmaxV = 1;
  for (let vy = -VMAX; vy <= VMAX; vy++) { const row: number[] = []; for (let vx = -VMAX; vx <= VMAX; vx++) {
    const z = stalkValue(field, qx, qy, vx, vy); row.push(z); if (isFinite(z)) vmaxV = Math.max(vmaxV, z); } vals.push(row); }
  const scaleMax = Math.max(vmaxV, ...(field === "ep" || field === "lower" ? [Math.max(...vals.flat().filter(isFinite))] : []));
  const gx = (vx: number) => pad + (vx + VMAX) * cs, gy = (vy: number) => pad + (VMAX - vy) * cs;
  for (let vy = -VMAX; vy <= VMAX; vy++) for (let vx = -VMAX; vx <= VMAX; vx++) {
    const z = vals[vy + VMAX][vx + VMAX];
    sctx.fillStyle = isFinite(z) ? rgb(ramp(BLUE, z / scaleMax, dark)) : css("--surface");
    sctx.fillRect(gx(vx), gy(vy), cs, cs);
    if (!isFinite(z)) { sctx.strokeStyle = css("--rule"); sctx.strokeRect(gx(vx) + 0.5, gy(vy) + 0.5, cs - 1, cs - 1); }
    if (cs > 26 && isFinite(z)) {
      sctx.fillStyle = z / scaleMax > 0.55 !== dark ? "#ffffff" : "#0b0b0b";
      sctx.font = `13px ${css("--font-mono")}`; sctx.textAlign = "center"; sctx.textBaseline = "middle";
      sctx.fillText(String(z), gx(vx) + cs / 2, gy(vy) + cs / 2);
    }
  }
  // corridor: hatch velocities whose (cell, mode) is outside it
  const r = run(), cor = corridorOf(query(), r), ci = (qy - 1) * N + (qx - 1);
  sctx.strokeStyle = dark ? "rgba(255,255,255,0.55)" : "rgba(11,11,11,0.45)"; sctx.lineWidth = 1.2;
  for (let vy = -VMAX; vy <= VMAX; vy++) for (let vx = -VMAX; vx <= VMAX; vx++) {
    const inside = cor.cells[ci] === 1 && (!cor.modes || (cor.modes.get(ci)?.has(modeOf(data.stalk.modes16, vx, vy)) ?? false));
    if (inside) continue;
    sctx.save();
    sctx.beginPath(); sctx.rect(gx(vx), gy(vy), cs, cs); sctx.clip();
    sctx.beginPath();
    for (let t = -cs; t < cs; t += 7) { sctx.moveTo(gx(vx) + t, gy(vy) + cs); sctx.lineTo(gx(vx) + t + cs, gy(vy)); }
    sctx.stroke();
    sctx.restore();
  }
  // mode outlines
  const modes = field === "cert" ? data.stalk.modes_cert : pSel.value === "blind" && field !== "cert" ? null : data.stalk.modes16;
  sctx.strokeStyle = css("--ink"); sctx.lineWidth = 2.2;
  sctx.beginPath();
  if (modes) {
    for (let vy = -VMAX; vy <= VMAX; vy++) for (let vx = -VMAX; vx <= VMAX; vx++) {
      const mHere = modeOf(modes, vx, vy);
      if (vx < VMAX && modeOf(modes, vx + 1, vy) !== mHere) { sctx.moveTo(gx(vx + 1), gy(vy)); sctx.lineTo(gx(vx + 1), gy(vy) + cs); }
      if (vy < VMAX && modeOf(modes, vx, vy + 1) !== mHere) { sctx.moveTo(gx(vx), gy(vy)); sctx.lineTo(gx(vx) + cs, gy(vy)); }
    }
  }
  sctx.rect(pad, pad, NV * cs, NV * cs);
  sctx.stroke();
  // the run's velocity at this cell
  r.traj.forEach(([x, y, vx, vy], k) => {
    if (x !== qx || y !== qy) return;
    sctx.strokeStyle = css("--s-ctl"); sctx.lineWidth = k === Math.min(step, r.traj.length - 1) ? 4 : 2.5;
    sctx.beginPath(); sctx.arc(gx(vx) + cs / 2, gy(vy) + cs / 2, cs * 0.38, 0, 2 * Math.PI); sctx.stroke();
  });
  // axes
  sctx.fillStyle = css("--ink-2"); sctx.font = `14px ${css("--font-mono")}`; sctx.textAlign = "center"; sctx.textBaseline = "alphabetic";
  for (let v = -VMAX; v <= VMAX; v += 3) { sctx.fillText(String(v), gx(v) + cs / 2, W - pad + 20); }
  sctx.textAlign = "right";
  for (let v = -VMAX; v <= VMAX; v += 3) sctx.fillText(String(v), pad - 8, gy(v) + cs / 2 + 5);
  sctx.textAlign = "left"; sctx.fillText("vₓ →", W - pad - 40, W - 6); sctx.fillText("v_y ↑", 4, pad - 14);

  const names: Record<string, string> = { vstar: "V*", ep: "E·P·V*", lower: "E·Vc (lower bound)", cert: "∀∃ certificate" };
  $("stalk-legend").innerHTML =
    `<span class="key"><span class="ramp" style="background:linear-gradient(90deg, ${(dark ? [...BLUE].reverse() : BLUE).join(",")})"></span><span class="num">${names[field]}: 0–${scaleMax}</span></span>` +
    `<span class="key"><span class="swatch box" style="background:var(--surface);outline:1px solid var(--rule)"></span>${field === "cert" ? "not certified" : "cannot reach the goal"}</span>` +
    `<span class="key"><span class="swatch box" style="outline:2px solid var(--ink)"></span>mode boundaries${field === "cert" ? " (8 × 6)" : pSel.value === "blind" ? " (one mode)" : " (16 × 3)"}</span>` +
    `<span class="key"><span class="swatch box" style="border:2px solid var(--s-ctl);border-radius:50%"></span>the run's velocity here</span>`;
}
const stalkTip = $("stalk-tip");
stalkC.addEventListener("mousemove", ev => {
  const rect = stalkC.getBoundingClientRect(), W = stalkC.width, pad = 40, cs = (W - 2 * pad) / NV;
  const px = (ev.clientX - rect.left) / rect.width * W, py = (ev.clientY - rect.top) / rect.height * W;
  const vx = Math.floor((px - pad) / cs) - VMAX, vy = VMAX - Math.floor((py - pad) / cs);
  if (Math.abs(vx) > VMAX || Math.abs(vy) > VMAX || !vstar) { stalkTip.hidden = true; return; }
  const [qx, qy] = currentCell(), f = (z: number) => isFinite(z) ? String(z) : "∞";
  const m16 = modeOf(data.stalk.modes16, vx, vy), mc = modeOf(data.stalk.modes_cert, vx, vy);
  stalkTip.innerHTML = `v = (${vx}, ${vy}) · mode ${pSel.value === "blind" ? 1 : m16} (cert ${mc})<br>` +
    `V* ${f(stalkValue("vstar", qx, qy, vx, vy))} · E·P·V* ${f(stalkValue("ep", qx, qy, vx, vy))}<br>` +
    `E·Vc ${f(stalkValue("lower", qx, qy, vx, vy))} · ∀∃ ${f(stalkValue("cert", qx, qy, vx, vy))}`;
  placeTip(stalkTip, ev, 240);
});
stalkC.addEventListener("mouseleave", () => { stalkTip.hidden = true; });

// ---- the run in phase space ---------------------------------------------------------------
function drawRunPhase(): void {
  const r = run(), host = $("runphase");
  host.innerHTML = "";
  const gxc = data.goal.reduce((s, g) => s + g[0], 0) / data.goal.length;
  const gyc = data.goal.reduce((s, g) => s + g[1], 0) / data.goal.length;
  const ptsX = r.traj.map(s => [s[0] - gxc, s[2]] as [number, number]);
  const ptsY = r.traj.map(s => [s[1] - gyc, s[3]] as [number, number]);
  const W = 560, H = 300, m = { l: 40, r: 14, t: 12, b: 34 };
  const all = [...ptsX, ...ptsY].map(p => p[0]);
  const xmin = Math.min(-20, ...all) - 2, xmax = Math.max(20, ...all) + 2;
  const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Trajectory projections in phase space" }, host);
  const X = (q: number) => m.l + (q - xmin) / (xmax - xmin) * (W - m.l - m.r);
  const Y = (v: number) => m.t + (VMAX + 0.5 - v) / (NV) * (H - m.t - m.b);
  for (let v = -VMAX; v <= VMAX; v += 2) {
    el("line", { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v), class: "gridline" }, svg);
    el("text", { x: m.l - 6, y: Y(v) + 4, "text-anchor": "end" }, svg).textContent = String(v);
  }
  const tickStep = (xmax - xmin) > 150 ? 40 : 20;
  for (let q = Math.ceil(xmin / tickStep) * tickStep; q <= xmax; q += tickStep)
    el("text", { x: X(q), y: H - m.b + 16, "text-anchor": "middle" }, svg).textContent = String(q);
  el("text", { x: (m.l + W - m.r) / 2, y: H - 2, "text-anchor": "middle" }, svg).textContent = "position − goal (cells)";
  el("line", { x1: X(0), x2: X(0), y1: m.t, y2: H - m.b, stroke: css("--rule") }, svg);
  // the lattice switching curve: from (q, v), full braking stops after |v|(|v| − 1)/2 cells
  let d = "";
  for (let v = -VMAX; v <= VMAX; v++) d += `${v === -VMAX ? "M" : "L"}${X(-Math.sign(v) * Math.abs(v) * (Math.abs(v) - 1) / 2)},${Y(v)}`;
  el("path", { d, fill: "none", stroke: css("--ink-3"), "stroke-width": 2, "stroke-dasharray": "4 3" }, svg);
  const series = (pts: [number, number][], color: string) => {
    let p = "";
    pts.forEach(([q, v], k) => { p += `${k === 0 ? "M" : "L"}${X(q)},${Y(v)}`; });
    el("path", { d: p, fill: "none", stroke: color, "stroke-width": 2, "stroke-linejoin": "round" }, svg);
    pts.forEach(([q, v], k) => el("circle", { cx: X(q), cy: Y(v), r: k === Math.min(step, pts.length - 1) ? 6 : 3,
      fill: k === Math.min(step, pts.length - 1) ? css("--surface") : color, stroke: color, "stroke-width": 2 }, svg));
  };
  series(ptsY, css("--s-route"));
  series(ptsX, css("--s-ctl"));
  const tip = $("runphase-tip");
  const hit = el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: "transparent" }, svg);
  hit.addEventListener("mousemove", (ev: Event) => {
    const e = ev as MouseEvent, rect = (svg as SVGSVGElement).getBoundingClientRect();
    const px = (e.clientX - rect.left) / rect.width * W, py = (e.clientY - rect.top) / rect.height * H;
    let best = 0, dist = Infinity, which = "x";
    [ptsX, ptsY].forEach((pts, j) => pts.forEach(([q, v], k) => {
      const dd = (X(q) - px) ** 2 + (Y(v) - py) ** 2;
      if (dd < dist) { dist = dd; best = k; which = j === 0 ? "x" : "y"; }
    }));
    if (dist > 400) { tip.hidden = true; return; }
    const s = r.traj[best];
    tip.innerHTML = `step ${best}, ${which} axis<br>q − goal = ${which === "x" ? (s[0] - gxc).toFixed(0) : (s[1] - gyc).toFixed(0)}, v = ${which === "x" ? s[2] : s[3]}`;
    placeTip(tip, e, 200);
  });
  hit.addEventListener("mouseleave", () => { tip.hidden = true; });
}

// ---- the one-axis phase portrait -------------------------------------------------------
function drawPhase(): void {
  const c = $<HTMLCanvasElement>("phase"), ctx = c.getContext("2d") as CanvasRenderingContext2D;
  const P = data.phase, nq = P.nq;
  const pad = { l: 44, r: 12, t: 12, b: 40 };
  const w = c.width - pad.l - pad.r, h = c.height - pad.t - pad.b;
  const X = (q: number) => pad.l + (q - 0.5) / nq * w;
  const Y = (v: number) => pad.t + (VMAX + 0.5 - v) / NV * h;
  ctx.clearRect(0, 0, c.width, c.height);
  const fills: Record<number, string> = { 1: css("--acc-plus"), [-1]: css("--acc-minus"), 0: css("--acc-zero") };
  for (let iv = 0; iv < NV; iv++) for (let q = 0; q < nq; q++) {
    const a = P.accel[iv * nq + q], v = iv - VMAX;
    ctx.fillStyle = P.value[iv * nq + q] < 0 ? css("--surface") : fills[a];
    ctx.fillRect(X(q + 1) - w / nq / 2, Y(v) - h / NV / 2, w / nq + 0.5, h / NV + 0.5);
  }
  const val = (x: number, y: number) => { const z = P.value[y * nq + x]; return z < 0 ? Infinity : z; };
  ctx.strokeStyle = isDark() ? "rgba(255,255,255,0.5)" : "rgba(11,11,11,0.35)"; ctx.lineWidth = 1;
  const vmaxVal = Math.max(...P.value);
  for (let level = 2; level < vmaxVal; level += 2) {
    ctx.beginPath();
    for (const [x0, y0, x1, y1] of contours(val, nq, NV, level)) {
      ctx.moveTo(X(x0 + 1), Y(y0 - VMAX)); ctx.lineTo(X(x1 + 1), Y(y1 - VMAX));
    }
    ctx.stroke();
  }
  ctx.fillStyle = css("--ink-2"); ctx.font = `22px ${css("--font-mono")}`; ctx.textAlign = "center";
  for (let q = 0; q <= nq; q += 16) ctx.fillText(String(q - P.goal), X(Math.max(q, 1)), c.height - 12);
  ctx.textAlign = "right";
  for (let v = -VMAX; v <= VMAX; v += 3) ctx.fillText(String(v), pad.l - 8, Y(v) + 7);
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
    const q = Math.round((px - pad.l) / w * nq + 0.5), v = Math.round(VMAX + 0.5 - (py - pad.t) / h * NV - 0.5);
    if (q < 1 || q > nq || v < -VMAX || v > VMAX) { tip.hidden = true; return; }
    const i = (v + VMAX) * nq + (q - 1), a = P.accel[i];
    tip.innerHTML = `q − goal = ${q - P.goal}, v = ${v}<br>minimum time ${P.value[i] < 0 ? "–" : P.value[i]} · a* = ${a > 0 ? "+1" : a < 0 ? "−1" : "0"}`;
    placeTip(tip, ev, 210);
  };
  c.onmouseleave = () => { tip.hidden = true; };
}

// ---- trade-off charts ---------------------------------------------------------------------
const PLANNERS = [
  { key: "velocity-blind", label: "velocity-blind", color: "--s-opt" },
  { key: "8 headings × 3 speeds", label: "8 headings × 3 speeds", color: "--s-ctl" },
  { key: "16 headings × 3 speeds", label: "16 headings × 3 speeds", color: "--s-low" }];
function drawTradeoffs(): void {
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
    el("text", { x: m.l - 6, y: y(t) + 4, "text-anchor": "end" }, svg).textContent = `${Math.round(t * 100)}%`;
  }
  used.forEach(s => {
    el("text", { x: x(s), y: H - m.b + 14, "text-anchor": "end", transform: `rotate(-35 ${x(s)} ${H - m.b + 14})` }, svg).textContent = s;
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
        tip.innerHTML = `${se.label}, ${p.s}<br>first try ${(p.first * 100).toFixed(0)}% · cost/opt ${p.ratio.toFixed(3)} (worst ${p.worst.toFixed(3)})`;
        placeTip(tip, ev as MouseEvent, 300);
      });
      dot.addEventListener("mouseleave", () => { tip.hidden = true; });
    }
  }
  const host2 = $("bounds"); host2.innerHTML = "";
  const svg2 = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Guaranteed bound against true ratio" }, host2);
  const rows = data.certified;
  const groups = PLANNERS.map(p => ({ ...p, rows: rows.filter(r => r["planner"] === p.key) }));
  const ymax = Math.max(...rows.map(r => +r["max_bound"]));
  const yb = (v: number) => H - m.b - (v - 1) / (ymax - 1) * (H - m.t - m.b);
  for (let t = 1; t <= ymax + 1e-9; t += ymax > 4 ? 1 : 0.5) {
    el("line", { x1: m.l, x2: W - m.r, y1: yb(t), y2: yb(t), class: "gridline" }, svg2);
    el("text", { x: m.l - 6, y: yb(t) + 4, "text-anchor": "end" }, svg2).textContent = `${t}×`;
  }
  const slots = groups.flatMap(g => g.rows.map(r => ({ g, r })));
  const xs = (i: number) => m.l + (i + 0.5) / slots.length * (W - m.l - m.r);
  slots.forEach(({ g, r }, i) => {
    const col = css(g.color);
    el("line", { x1: xs(i), x2: xs(i), y1: yb(+r["mean_bound"]), y2: yb(+r["max_bound"]), stroke: col, "stroke-width": 2, "stroke-opacity": 0.5 }, svg2);
    el("rect", { x: xs(i) - 6, y: yb(+r["mean_bound"]) - 2, width: 12, height: 4, rx: 1, fill: col }, svg2);
    el("circle", { cx: xs(i), cy: yb(+r["mean_cost_over_optimal"]), r: 4.5, fill: css("--surface"), stroke: col, "stroke-width": 2 }, svg2);
    el("text", { x: xs(i), y: H - m.b + 14, "text-anchor": "end", transform: `rotate(-35 ${xs(i)} ${H - m.b + 14})` }, svg2).textContent = `w=${r["width"]}`;
    const hit = el("rect", { x: xs(i) - 14, y: m.t, width: 28, height: H - m.t - m.b, fill: "transparent" }, svg2);
    hit.addEventListener("mouseenter", (ev: Event) => {
      tip.innerHTML = `${g.label}, w=${r["width"]} m=${r["mode_width"]}<br>true ${(+r["mean_cost_over_optimal"]).toFixed(3)} · bound mean ${(+r["mean_bound"]).toFixed(2)}, worst ${(+r["max_bound"]).toFixed(2)}`;
      placeTip(tip, ev as MouseEvent, 300);
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

// ---- animation and wiring -------------------------------------------------------------
let timer: number | null = null;
function setStep(k: number): void {
  step = Math.max(0, Math.min(k, run().traj.length - 1));
  stepIn.value = String(step);
  const s = run().traj[step];
  $("step-label").textContent = `step ${step} / ${run().traj.length - 1} · v = (${s[2]}, ${s[3]})`;
  drawMap(); drawSandwich(); drawRunPhase(); drawStalk();
}
function stop(): void { if (timer !== null) { clearInterval(timer); timer = null; } playBtn.textContent = "Play"; playBtn.setAttribute("aria-pressed", "false"); }
playBtn.addEventListener("click", () => {
  if (timer !== null) { stop(); return; }
  if (step >= run().traj.length - 1) setStep(0);
  playBtn.textContent = "Pause"; playBtn.setAttribute("aria-pressed", "true");
  const reduce = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  timer = window.setInterval(() => {
    if (step >= run().traj.length - 1) { stop(); return; }
    setStep(step + 1);
  }, reduce ? 900 : 380);
});
stepIn.addEventListener("input", () => { stop(); setStep(+stepIn.value); });
function render(): void {
  stop();
  stepIn.max = String(run().traj.length - 1);
  drawStats();
  setStep(Math.min(step, run().traj.length - 1));
}
qSel.addEventListener("change", () => { step = 0; pinned = null; render(); });
for (const s of [pSel, cSel]) s.addEventListener("change", render);
lSel.addEventListener("change", drawMap);
fSel.addEventListener("change", drawStalk);
for (const t of [T_CORR, T_ROUTE, T_CTL, T_OPT, T_ISO]) t.addEventListener("change", drawMap);
T_CMP.addEventListener("change", () => { drawMap(); drawSandwich(); });
followBtn.setAttribute("aria-pressed", "true");
const redrawAll = () => { render(); drawPhase(); drawTradeoffs(); };
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redrawAll);
new MutationObserver(redrawAll).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
redrawAll();
Promise.all([blob("blob-vstar"), blob("blob-vup"), blob("blob-lower16")]).then(([a, b, c]) => {
  vstar = a; vup = b; lower16 = c;
  drawStalk();
}).catch(err => {
  console.error(err);
  $("stalk-title").textContent = "Velocity stalk (this browser cannot decompress the stalk data)";
});
