#!/usr/bin/env python3
"""Local app: opens a window-like page in your web browser where you choose the analysis, point to a folder of
data files (any names), tick the analyses you want, and press Run. Everything stays on your computer."""
import argparse
import contextlib
import json
import mimetypes
import os
import sys
import threading
import traceback
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
sys.path.insert(0, str(ROOT))

from biomarker import autodetect, catalog, runner  # noqa: E402

OUT_ROOT = ROOT / "outputs"
JOBS = {}
_lock = threading.Lock()


class _Log:
    def __init__(self, job):
        self.job = job

    def write(self, s):
        if s:
            self.job["buf"] = self.job.get("buf", "") + s
            *done, self.job["buf"] = self.job["buf"].split("\n")
            self.job["lines"].extend(done)
        return len(s)

    def flush(self):
        pass


def _run(job):
    log = _Log(job)
    try:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            rep = runner.run_job(job["settings"], str(OUT_ROOT))
        job["report"] = "/" + str(Path(rep).resolve().relative_to(ROOT)).replace(os.sep, "/")
        job["folder"] = str(Path(rep).resolve().parent)
        job["state"] = "done"
    except Exception as e:                                  # shown to the user in plain words
        job["lines"].append("")
        job["lines"].append("".join(traceback.format_exception_only(type(e), e)).strip())
        job["error"] = str(e) or type(e).__name__
        job["trace"] = traceback.format_exc()
        job["state"] = "error"


class H(BaseHTTPRequestHandler):
    server_version = "BiomarkerStudio/1.0"

    def log_message(self, *a):
        pass

    def _json(self, obj, code=200):
        b = json.dumps(obj, default=str).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(b)))
        self.end_headers()
        self.wfile.write(b)

    def _body(self):
        n = int(self.headers.get("Content-Length") or 0)
        return self.rfile.read(n) if n else b""

    def _file(self, path, ctype=None):
        data = Path(path).read_bytes()
        self.send_response(200)
        self.send_header("Content-Type", ctype or mimetypes.guess_type(str(path))[0] or "application/octet-stream")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            if u.path == "/":
                return self._file(ROOT / "biomarker" / "ui" / "index.html", "text/html; charset=utf-8")
            if u.path == "/api/catalog":
                return self._json(dict(analyses=catalog.ANALYSES, modules=catalog.MODULES, roles=catalog.ROLES,
                                       default_folder=str(ROOT / "input"), outputs=str(OUT_ROOT)))
            if u.path == "/api/scan":
                return self._json(dict(files=autodetect.scan_folder(q["folder"])))
            if u.path == "/api/fs":
                p = Path(q.get("path") or ROOT).expanduser().resolve()
                if not p.is_dir():
                    p = p.parent
                dirs = sorted(d.name for d in p.iterdir() if d.is_dir() and not d.name.startswith("."))
                return self._json(dict(path=str(p), parent=str(p.parent), dirs=dirs,
                                       files=sum(1 for f in p.iterdir() if f.is_file())))
            if u.path == "/api/models":
                ms = sorted(OUT_ROOT.glob("*/model.joblib"), key=lambda m: m.stat().st_mtime, reverse=True) if OUT_ROOT.exists() else []
                return self._json(dict(models=[str(m) for m in ms]))
            if u.path == "/api/job":
                j = JOBS.get(q["id"])
                if not j:
                    return self._json(dict(error="unknown job"), 404)
                since = int(q.get("since", 0))
                return self._json(dict(state=j["state"], lines=j["lines"][since:], total=len(j["lines"]),
                                       report=j.get("report"), folder=j.get("folder"), error=j.get("error")))
            if u.path.startswith("/outputs/"):
                p = (ROOT / unquote(u.path.lstrip("/"))).resolve()
                if OUT_ROOT.resolve() in p.parents and p.is_file():
                    return self._file(p)
            self.send_error(404)
        except Exception as e:
            self._json(dict(error=f"{type(e).__name__}: {e}"), 400)

    def do_POST(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            if u.path == "/api/upload":
                folder = Path(q["folder"])
                folder.mkdir(parents=True, exist_ok=True)
                dest = folder / Path(unquote(q["name"])).name
                n = int(self.headers.get("Content-Length") or 0)
                with open(dest, "wb") as f:
                    left = n
                    while left > 0:
                        chunk = self.rfile.read(min(1 << 20, left))
                        if not chunk:
                            break
                        f.write(chunk)
                        left -= len(chunk)
                return self._json(dict(ok=True, path=str(dest)))
            body = json.loads(self._body() or b"{}")
            if u.path == "/api/groups":
                return self._json(autodetect.group_info(body["expression"], body.get("sheet"), body.get("group_column")))
            if u.path == "/api/run":
                with _lock:
                    if any(j["state"] == "running" for j in JOBS.values()):
                        return self._json(dict(error="An analysis is already running."), 409)
                    jid = str(len(JOBS) + 1)
                    JOBS[jid] = job = dict(state="running", lines=[], settings=body)
                threading.Thread(target=_run, args=(job,), daemon=True).start()
                return self._json(dict(id=jid))
            if u.path == "/api/open":
                p = Path(body["path"]).resolve()
                if OUT_ROOT.resolve() in [p, *p.parents]:
                    if sys.platform.startswith("win"):
                        os.startfile(str(p))                      # noqa: S606
                    elif sys.platform == "darwin":
                        os.system(f'open "{p}"')                  # noqa: S605
                    else:
                        os.system(f'xdg-open "{p}" >/dev/null 2>&1 &')  # noqa: S605
                return self._json(dict(ok=True))
            if u.path == "/api/quit":
                self._json(dict(ok=True))
                threading.Thread(target=self.server.shutdown, daemon=True).start()
                return
            self.send_error(404)
        except Exception as e:
            self._json(dict(error=f"{type(e).__name__}: {e}"), 400)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--no-browser", action="store_true")
    a = ap.parse_args()
    (ROOT / "input").mkdir(exist_ok=True)
    srv = None
    for port in range(a.port, a.port + 20):
        try:
            srv = ThreadingHTTPServer(("127.0.0.1", port), H)
            break
        except OSError:
            continue
    if srv is None:
        sys.exit("Could not open a local port. Close other copies of the app and try again.")
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    print(f"\n  RNA-seq Biomarker Studio is running at {url}\n  (It opens in your web browser. Close this window or press Ctrl+C to stop.)\n")
    if not a.no_browser:
        threading.Timer(0.8, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    print("Stopped.")


if __name__ == "__main__":
    main()
