import subprocess
from pathlib import Path
from config import QMD_PATH, OUT_DIR, DATA_PATH

def render_report(qmd_path: str, out_dir: str, data_path: str, query: str):
    qmd = Path(qmd_path)
    out = Path(out_dir)

    cmd = [
        "quarto", "render", str(qmd),
        "--output-dir", str(out),
        "-P", f"data_path:{data_path}",
        "-P", f"query:{query}",
    ]

    # check=True raises if quarto fails; capture_output gives you logs
    result = subprocess.run(cmd, check=True, text=True, capture_output=True)
    return result.stdout, result.stderr


stdout, stderr = render_report(
        qmd_path=QMD_PATH,
        out_dir=OUT_DIR,
        data_path=DATA_PATH,
        query="region='EU' AND date >= '2026-01-01'",
    )

print(stdout)