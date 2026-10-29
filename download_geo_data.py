
from pathlib import Path
import requests
import GEOparse

PROJECT_DIR = Path(__file__).resolve().parent
RAW_DIR = PROJECT_DIR / "data" / "raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)

GEO_ACCESSION = "GSE10072"
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/"
    "GSE10nnn/GSE10072/matrix/GSE10072_series_matrix.txt.gz"
)
MATRIX_FILE = RAW_DIR / "GSE10072_series_matrix.txt.gz"


def download_matrix():
    if MATRIX_FILE.exists() and MATRIX_FILE.stat().st_size > 0:
        print(f"Series matrix already exists: {MATRIX_FILE}")
        return

    temporary_file = MATRIX_FILE.with_suffix(".gz.part")

    print(f"Downloading series matrix from GEO:\n{MATRIX_URL}")

    try:
        with requests.get(MATRIX_URL, stream=True, timeout=(30, 180)) as response:
            response.raise_for_status()

            with temporary_file.open("wb") as output:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        output.write(chunk)

        if temporary_file.stat().st_size == 0:
            raise RuntimeError("Downloaded series matrix is empty.")

        temporary_file.replace(MATRIX_FILE)
        print(f"Saved series matrix: {MATRIX_FILE}")

    finally:
        if temporary_file.exists():
            temporary_file.unlink()


def main():
    print(f"Preparing GEO dataset {GEO_ACCESSION}")

    # Download the brief GEO series record used by the metadata scripts.
    GEOparse.get_GEO(
        geo=GEO_ACCESSION,
        destdir=str(RAW_DIR),
        how="brief",
    )

    download_matrix()

    print("\nDownload preparation complete.")
    print(f"Raw data directory: {RAW_DIR}")
    print("The GPL96 platform annotation is prepared separately.")


if __name__ == "__main__":
    main()
