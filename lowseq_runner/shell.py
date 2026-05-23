import shlex
import subprocess
from pathlib import Path
from typing import Sequence


def quote_cmd(cmd: Sequence[str]) -> str:
    return " ".join(shlex.quote(str(x)) for x in cmd)


def run_to_log(cmd: Sequence[str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)

    with log_path.open("w") as log_file:
        proc = subprocess.Popen(
            [str(x) for x in cmd],
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return proc.wait()
