"""Download TCGA LUAD (PanCanAtlas) mutations + clinical data.

The datahub MAF is a broken Git-LFS pointer (object missing on the server),
so mutations are fetched from the cBioPortal REST API
(POST /molecular-profiles/{id}/mutations/fetch, same route as the JHU cohort)
and written as a MAF-like TSV. Clinical tables come from the datahub raw URLs.
Output: data/raw/luad_data_mutations.txt, luad_data_clinical_sample.txt,
luad_data_clinical_patient.txt.
"""
from pathlib import Path
import requests

STUDY = "luad_tcga_pan_can_atlas_2018"
PROFILE = f"{STUDY}_mutations"
API = "https://www.cbioportal.org/api"
BASE_RAW = f"https://raw.githubusercontent.com/cBioPortal/datahub/master/public/{STUDY}"
RAW = Path(__file__).resolve().parents[0] / "data" / "raw"
HEAD = {"Accept": "application/json", "Content-Type": "application/json"}

COLUMNS = ["Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Start_Position",
           "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
           "Variant_Classification", "Variant_Type", "Protein_position",
           "HGVSp_Short", "Transcript_ID", "NCBI_Build"]


def fetch_clinical(name):
    out = RAW / f"luad_{name}"
    if out.exists() and out.stat().st_size > 1000:
        print(f"exists: {out.name} ({out.stat().st_size} bytes)")
        return
    for url in (f"https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/{STUDY}/{name}",
                f"{BASE_RAW}/{name}"):
        r = requests.get(url, timeout=300)
        if r.status_code == 200 and not r.text.startswith("version https://git-lfs"):
            break
    else:
        print(f"FAILED clinical {name}")
        return
    out.write_bytes(r.content)
    print(f"downloaded: {out.name} ({len(r.content)} bytes)")


def all_sample_ids():
    ids = []
    page = 0
    while True:
        r = requests.get(f"{API}/studies/{STUDY}/samples",
                         params={"projection": "ID", "pageSize": 1000, "pageNumber": page},
                         headers=HEAD, timeout=120)
        r.raise_for_status()
        batch = [s["sampleId"] for s in r.json()]
        if not batch:
            return ids
        ids.extend(batch)
        page += 1


def fetch_mutations(sample_ids):
    rows = []
    for i in range(0, len(sample_ids), 50):
        chunk = sample_ids[i:i + 50]
        r = requests.post(f"{API}/molecular-profiles/{PROFILE}/mutations/fetch?projection=DETAILED",
                          json={"sampleIds": chunk}, headers=HEAD, timeout=600)
        r.raise_for_status()
        rows.extend(r.json())
        print(f"  batch {i // 50 + 1}: +{len(r.json())} rows (total {len(rows)})", flush=True)
    return rows


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    out = RAW / "luad_data_mutations.txt"
    if not (out.exists() and out.stat().st_size > 1000):
        ids = all_sample_ids()
        print("samples:", len(ids))
        rows = fetch_mutations(ids)
        with open(out, "w", encoding="utf-8", newline="") as f:
            f.write("\t".join(COLUMNS) + "\n")
            for m in rows:
                f.write("\t".join(str(x) for x in [
                    m.get("gene", {}).get("hugoGeneSymbol", ""),
                    m.get("sampleId", ""),
                    m.get("chr", ""),
                    m.get("startPosition", ""),
                    m.get("endPosition", ""),
                    m.get("referenceAllele", ""),
                    m.get("variantAllele", ""),
                    m.get("mutationType", ""),
                    m.get("variantType", ""),
                    m.get("proteinPosStart", ""),
                    m.get("proteinChange", ""),
                    m.get("refseqMrnaId", ""),
                    m.get("ncbiBuild", ""),
                ]) + "\n")
        print(f"downloaded: {out.name} ({out.stat().st_size} bytes, {len(rows)} rows)")
    else:
        print(f"exists: {out.name} ({out.stat().st_size} bytes)")
    fetch_clinical("data_clinical_sample.txt")
    fetch_clinical("data_clinical_patient.txt")
    # smoking history (pack-years) is not in the PanCanAtlas clinical; it lives
    # in the Firehose Legacy study for the same TCGA patients
    if not (RAW / "luad_smoking_data_clinical_patient.txt").exists():
        r = requests.get("https://media.githubusercontent.com/media/cBioPortal/datahub"
                         "/master/public/luad_tcga/data_clinical_patient.txt", timeout=300)
        if r.status_code == 200 and not r.text.startswith("version https://git-lfs"):
            (RAW / "luad_smoking_data_clinical_patient.txt").write_bytes(r.content)
            print(f"downloaded: luad_smoking_data_clinical_patient.txt ({len(r.content)} bytes)")
        else:
            print("FAILED smoking clinical")
    print("done")


if __name__ == "__main__":
    main()
