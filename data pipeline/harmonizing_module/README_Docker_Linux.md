# Harmonizing Module - Linux Docker Container

This module has been configured to run as a Linux Docker container.

## Key Changes Made for Linux Compatibility

### 1. **Dockerfile Changes**
- Changed from `python:3.9-windowsservercore-ltsc2022` to `python:3.9-slim` (Linux base)
- Updated all paths from Windows format (`C:\path`) to Linux format (`/path`)
- Removed `.exe` extensions from PLINK binary paths
- Added system dependencies and security improvements

### 2. **Configuration Changes**
- Updated `config.py` default paths to use Linux filesystem structure
- Changed PLINK executable names from `plink.exe`/`plink2.exe` to `plink`/`plink2`
- Updated temp directory from `C:\Users\frezz\pipeline_testing` to `/tmp`

### 3. **File Structure**
```
/app/                          # Application code
├── dependencies/              # PLINK binaries and test files
│   ├── plink                 # Linux PLINK 1.9 binary (no .exe)
│   ├── plink2                # Linux PLINK 2.0 binary (no .exe)
│   └── genome_*.txt          # Test genome file
└── main.py                   # Main script

/data/                        # Volume mount for reference data
├── harmonized_gwas_catalog_risk_snps.txt
├── plink.GRCh37.map/
├── hs37d5.fa.zst
├── beagle_references/
└── all_phase3.pvar.zst

/tmp/                         # Working directory for processing
```

## Building the Docker Image

```bash
docker build -t harmonizer-module .
```

## Running the Container

```bash
docker run -v /host/path/to/reference/data:/data \
           -v /host/path/to/output:/tmp \
           harmonizer-module
```

## Environment Variables

You can override the default paths using environment variables:

```bash
docker run -e TEMP_DIR=/custom/temp \
           -e PLINK_MAP_DIR=/data/custom/plink/maps \
           -v /host/data:/data \
           harmonizer-module
```

## Required Dependencies

Make sure your `/app/dependencies/` folder contains:
- `plink` (Linux binary, executable)
- `plink2` (Linux binary, executable)
- Test genome file

Make sure your `/data/` volume contains:
- Reference genome files
- PLINK map files
- Beagle reference files
- GWAS catalog data

## Notes

1. **PLINK Binaries**: You'll need Linux versions of PLINK 1.9 and 2.0. Download them from:
   - PLINK 1.9: https://www.cog-genomics.org/plink/
   - PLINK 2.0: https://www.cog-genomics.org/plink/2.0/

2. **File Permissions**: The container runs as a non-root user for security. Make sure volume-mounted directories have appropriate permissions.

3. **Memory**: Genomic data processing can be memory-intensive. Consider using `--memory` flag when running the container.


quick test: docker run -it -v "C:\Users\frezz\Desktop\dependencies:/data" -v "C:\Users\frezz\Desktop\genome_scraping\pipeline\SNPster\data pipeline\harmonizing_module\docker_test:/workspace" -w /workspace harmonizer-module:latest python /app/main.py