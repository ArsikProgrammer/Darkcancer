from pathlib import Path
import gzip, shutil

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "raw" / "reference" / "hg19.fa.gz"
DST = ROOT / "data" / "raw" / "reference" / "hg19.fa"

with gzip.open(SRC, "rb") as fin, DST.open("wb") as fout:
    shutil.copyfileobj(fin, fout, length=1024 * 1024)
print("Decompressed:", DST)
