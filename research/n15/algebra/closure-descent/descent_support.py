"""相对闭环消元的固定父输入与只追加输出接口。"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
SHARED = HERE.parent / "shared-radical"
sys.path.insert(0, str(SHARED))
from bridge_support import ROOT, TARGET, AXIS, FACTOR, read, write, sha, bounded, binary_environment
sys.path.insert(0, str(HERE))

SOURCE = SHARED / "system-runs/20261015-selected/circuit.json.gz"
