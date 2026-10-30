import requests, os
url = "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000005640_9606_HUMAN_v6.tar"
out = "data/raw/UP000005640_9606_HUMAN_v6.tar"
os.makedirs("data/raw", exist_ok=True)
resume = os.path.exists(out)
headers = {}
mode = "wb"
if resume:
    size = os.path.getsize(out)
    headers["Range"] = f"bytes={size}-"
    mode = "ab"
r = requests.get(url, stream=True, headers=headers, timeout=120)
print("status", r.status_code, "mode", mode, flush=True)
if resume and r.status_code == 200:
    raise SystemExit("server ignored Range; refusing to append (would corrupt tar)")
if not resume and r.status_code != 200:
    raise SystemExit(f"unexpected status {r.status_code}")
if resume and r.status_code not in (206, 416):
    raise SystemExit(f"unexpected status {r.status_code}")
if r.status_code == 416:
    print("already complete", flush=True)
    raise SystemExit(0)
with open(out, mode) as f:
    for chunk in r.iter_content(chunk_size=8*1024*1024):
        if chunk:
            f.write(chunk)
            f.flush()
print("done", flush=True)

