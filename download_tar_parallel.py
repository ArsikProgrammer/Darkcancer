import os, requests, threading
url = "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000005640_9606_HUMAN_v6.tar"
out = "data/raw/UP000005640_9606_HUMAN_v6.tar"
os.makedirs("data/raw", exist_ok=True)
total = 5177506304
print("total", total, flush=True)
parts = 4
chunk = total // parts
ranges = []
for i in range(parts):
    start = i * chunk
    end = total - 1 if i == parts - 1 else (i + 1) * chunk - 1
    ranges.append((i, start, end))

def worker(i, start, end):
    part = f"{out}.part{i}"
    if os.path.exists(part) and os.path.getsize(part) == (end - start + 1):
        print(f"part {i} already complete", flush=True)
        return
    headers = {"Range": f"bytes={start}-{end}"}
    r = requests.get(url, headers=headers, stream=True, timeout=120)
    if r.status_code != 206:
        raise RuntimeError(f"part {i}: unexpected status {r.status_code}")
    with open(part, "wb") as f:
        for c in r.iter_content(chunk_size=8*1024*1024):
            if c:
                f.write(c)
    if os.path.getsize(part) != (end - start + 1):
        raise RuntimeError(f"part {i}: size mismatch")
    print(f"part {i} done", flush=True)

threads = []
for i, s, e in ranges:
    t = threading.Thread(target=worker, args=(i, s, e))
    threads.append(t)
    t.start()
for t in threads:
    t.join()

for i, (start, end) in enumerate(ranges):
    part = f"{out}.part{i}"
    expected = end - start + 1
    if not os.path.exists(part) or os.path.getsize(part) != expected:
        raise SystemExit(f"part {i} incomplete; not touching {out}")

with open(out, "wb") as out_f:
    for i in range(parts):
        with open(f"{out}.part{i}", "rb") as p:
            while True:
                b = p.read(8*1024*1024)
                if not b:
                    break
                out_f.write(b)
print("combined", os.path.getsize(out), flush=True)
for i in range(parts):
    os.remove(f"{out}.part{i}")
print("done", flush=True)