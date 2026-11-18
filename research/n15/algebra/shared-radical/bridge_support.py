"""共享判别式实验的固定父输入与只追加输出接口。"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent / "rational-chain"
sys.path.insert(0, str(PARENT))
from support import ROOT, TARGET, AXIS, FACTOR, read, write, sha, bounded, binary_environment
sys.path.insert(0, str(HERE))

PARENT_SYSTEM = PARENT / "system-runs/20261015-critical10-final/circuit.json.gz"
