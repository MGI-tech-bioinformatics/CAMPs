# CAMPs

CAMPs designs multiplex PCR primers for bisulfite methylation targets.

Primer terminology follows the manuscript: GSP (formerly USP), snGSP (formerly DSP), and tagGSP (formerly TSP). Detailed Excel columns, primer IDs and alignment filenames use these new names. Existing workbooks must have their USP/DSP/TSP column prefixes updated before reuse. The distance parameter is now `--snGSP2tar`.

## Requirements

- Linux or Windows with WSL (the pipeline invokes Unix tools including awk, grep, cp, cat and zip)
- Python 3.9+
- Python packages in `requirements.txt`
- `bwa` and `samtools` on `PATH`
- The external `primersDB` executable and its database/kernel files
- The external `catPos` executable placed in `src/` with execute permission
- Reference, SNP and CpG index files configured in `src/conf.ini`

The repository does not include genome references or the external primersDB binary. Set those paths for your environment before running a design.

## Run

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
cp src/conf.ini conf.local.ini
# Edit conf.local.ini before running a design.
python src/CAMPs.py --help
python src/CAMPs.py --config conf.local.ini -i examples/targets.tsv -o camp_demo --ref_name hg19 --kernel CPU
```

Input rows are `chromosome<TAB>start<TAB>end<TAB>name`. Use `python CAMPs.py --help` for all options. For a custom reference, pass `--ref_name contol_ref --ref_path /path/to/reference-directory`.

Results are written to the output directory, including `Output.xlsx`, BED files and an HTML summary report.

## Configuration

Copy `src/conf.ini` to `conf.local.ini`, update `[ref]`, `[snp]`, `[cg]` and `[software]`, and pass `--config conf.local.ini`. Paths are relative to the current working directory; absolute paths are recommended. Use paths and output names without spaces or shell metacharacters because the pipeline invokes shell commands. `--ref_path` and `--snp_dir` override the configuration defaults. BWA and samtools must also be on PATH.

For a reference `hg19.fa`, `python src/genome_convert.py /path/to/hg19` generates `/path/to/hg19.tran.fa` with CT and GA contigs. This helper loads the entire genome into memory. Build its BWA index with `bwa index /path/to/hg19.tran.fa`, then configure `[ref] hg19` to this file. SNP data require chromosome-specific `<chromosome>.snp` files matching `PrepareData.load_snp`; CpG-number input additionally requires the CpG index. These datasets are not bundled. GPU mode requires `primersDB_CL`, `kernel.cl` and a compatible OpenCL runtime. No public download source or installer for primersDB/catPos is provided.

The example targets illustrate input format only. Interval ends are treated as inclusive in the pipeline; BED-named outputs should not be assumed to use standard BED half-open coordinates. Results include `<order_num>.html`, `CAMPs.log` and a ZIP archive. The legacy `src/ATOPlexMethDesigner.py` entry point remains supported.

## License

The original contact-based license notice is retained. No open-source license has been selected. The rights holder should choose redistribution terms before public release. External tools and datasets have their own terms. See [PUBLISHING.md](PUBLISHING.md). Full biological validation requires a complete run with external dependencies and references.
