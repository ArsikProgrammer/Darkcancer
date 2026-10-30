"""Resume download of AlphaFold human proteome tar via FTP."""
import ftplib
import os
from pathlib import Path

LOCAL = Path("data/raw/UP000005640_9606_HUMAN_v6.tar")
REMOTE = "UP000005640_9606_HUMAN_v6.tar"
DIR = "/pub/databases/alphafold/v6"

def main():
    ftp = ftplib.FTP("ftp.ebi.ac.uk")
    ftp.login()
    ftp.cwd(DIR)
    size = ftp.size(REMOTE)
    local_size = LOCAL.stat().st_size if LOCAL.exists() else 0
    print(f"server={size}, local={local_size}, remaining={size-local_size}", flush=True)
    if local_size >= size:
        print("Already complete.")
        ftp.quit()
        return
    # resume
    ftp.sendcmd(f"REST {local_size}")
    with open(LOCAL, "ab") as f:
        ftp.retrbinary(f"RETR {REMOTE}", f.write, blocksize=1024*1024)
    ftp.quit()
    print("Download finished. Size:", LOCAL.stat().st_size)

if __name__ == "__main__":
    main()
