"""有理接触链实验复用已有输入、资源限制和只追加产物接口。"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "chain-radius"))
from common import ROOT, TARGET, AXIS, read, write, sha, bounded, binary_environment
sys.path.insert(0, str(HERE))

FACTOR = HERE.parent / "three-parameter-elimination/c0-critical/factor-branch-verified.json"
