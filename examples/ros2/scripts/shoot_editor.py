#!/usr/bin/env python3
"""Screenshot the scenario editor, driving the real page in headless Chrome.

The editor is one self-contained HTML file, so a shot is taken by appending a
short bootstrap script to a copy of it: load a scenario, place the camera, run
the preview solver forward to an interesting moment, then let Chrome capture it.
Nothing is mocked, so the images cannot drift from what the tool actually does.

    python3 scripts/shoot_editor.py            # writes docs/img/*.png
"""
from __future__ import annotations
import json, os, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
REPO = os.path.dirname(os.path.dirname(ROOT))
EDITOR = os.path.join(ROOT, "tools", "swarm_editor.html")
OUT = os.path.join(REPO, "docs", "img")

BOOT = """
<script>
(async () => {
  const S = %s;
  try {
    if (S.cfg) { io.value = JSON.stringify(S.cfg); loadFromText(); }
    if (S.tab) document.querySelectorAll("#tabs .tab").forEach(b => {
      if (b.dataset.tab === S.tab) b.click(); });
    if (S.mode) document.querySelectorAll("button.mode").forEach(b => {
      if (b.dataset.mode === S.mode) b.click(); });
    if (S.select) selected = S.select;
    if (S.cam) Object.assign(cam, S.cam);
    if (S.t) {           // run the preview solver to a chosen moment
      anim.t = 0; simReset();
      const dt = 0.02;
      for (let u = 0; u < S.t; u += dt) { simStep(dt); anim.t += dt; }
      anim.playing = true;              // so the readout and slots show
    }
    if (S.prompt) document.getElementById("nlPrompt").value = S.prompt;
    if (S.status) document.getElementById("nlStatus").innerHTML = S.status;
    render();
  } catch (e) { document.title = "SHOOT ERROR: " + e.message; }
  document.title = "ready";
})();
</script>
"""

def shoot(name, state, w=1600, h=1000):
    html = open(EDITOR).read()
    page = html.replace("</body>", BOOT % json.dumps(state) + "</body>")
    tmp = tempfile.NamedTemporaryFile("w", suffix=".html", dir=os.path.join(ROOT, "tools"),
                                      delete=False)
    tmp.write(page); tmp.close()
    dest = os.path.join(OUT, name + ".png")
    cmd = ["google-chrome", "--headless=new", "--disable-gpu", "--hide-scrollbars",
           f"--window-size={w},{h}", f"--screenshot={dest}",
           "--virtual-time-budget=6000", "--no-sandbox", "file://" + tmp.name]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    os.unlink(tmp.name)
    ok = os.path.exists(dest) and os.path.getsize(dest) > 20000
    print(f"  {'ok ' if ok else 'FAIL'} {name}.png"
          + (f"  ({os.path.getsize(dest)//1024} kB)" if ok else "  " + r.stderr[-200:]))
    return ok

def main():
    os.makedirs(OUT, exist_ok=True)
    cfgs = json.load(open(os.path.join(HERE, "shots.json")))
    good = 0
    for shot in cfgs:
        good += bool(shoot(shot.pop("name"), shot))
    print(f"{good}/{len(cfgs)} screenshots written to docs/img/")
    return 0 if good == len(cfgs) else 1

if __name__ == "__main__":
    sys.exit(main())
