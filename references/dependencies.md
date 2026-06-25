# Dependencies

## Python packages

```
biopython>=1.81
matplotlib>=3.7
numpy>=1.24
pandas>=2.0
```

Install with:
```bash
pip install biopython matplotlib numpy pandas
```

## Optional CLI tools

For higher-quality multiple sequence alignment:

- **MAFFT** — https://mafft.cbrc.jp/alignment/software/
  - On Debian/Ubuntu: `sudo apt install mafft`
  - On macOS: `brew install mafft`
  - On Windows: download the binary from the MAFFT site or use WSL.
- **MUSCLE** — https://www.drive5.com/muscle/
  - `conda install -c bioconda muscle` (recommended)

The pipeline automatically detects `mafft` on `PATH`. If neither `mafft` nor `muscle` is available, the pairwise identity for panel B is still computed (alignment of every sequence vs reference via Biopython's `pairwise2`), but the multi-sequence alignment for panel D will be skipped.

## NCBI access

- An email address is required by NCBI's E-utilities policy (Entrez.email). Pass it via `--email`.
- An optional API key raises the rate limit from 3 to 10 requests/second. Get one at https://www.ncbi.nlm.nih.gov/account/settings/.

## Internet access

- The script must reach `eutils.ncbi.nlm.nih.gov` on port 443.
- If you are behind a corporate proxy, set the standard `HTTP_PROXY` / `HTTPS_PROXY` environment variables.

## Tested environment

- Python 3.10 - 3.13
- Windows 10/11, macOS 13+, Ubuntu 20.04+
- Biopython 1.81+, matplotlib 3.7+, numpy 1.24+, pandas 2.0+
