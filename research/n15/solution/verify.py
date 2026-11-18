"""从归档输入重新验证 N=15 的全部数学交付，并核对公开格式的整数系数。"""

import argparse
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
ALG=ROOT/"n15/algebra/four-directions"
GLOBAL=ROOT/"n15/proofs/global-final"
sys.path.insert(0,str(ALG))
from chart import read
from export_solution import export
from verify import require


def check(folder,threads):
    folder=folder.resolve()
    folder.mkdir(parents=True,exist_ok=False)
    def run(script,*args):
        subprocess.run([sys.executable,"-X","utf8","-B",str(script),*map(str,args)],cwd=ROOT,check=True)
    candidate=ALG/"reconstruction-p256.json.gz"
    coordinates=ALG/"scalars-p256.json.gz"
    require(read(HERE/"coefficients.json")==read(candidate)["candidate_primitive_coefficients_ascending"],
            "交付系数与实际认证多项式不相同")
    run(GLOBAL/"replay_archived_inputs.py","--output",folder/"global.json")
    run(GLOBAL/"replay_archived_inputs.py","--axis-chain","--output",folder/"chain.json")
    run(ALG/"replay.py")
    run(ALG/"check_reconstructed.py","--source",candidate,"--output",folder/"factors.json")
    run(ALG/"check_scalar_identities.py","--output",folder/"scalar-identities.json")
    run(ALG/"scalar_target.py","--output",folder/"target.json")
    run(ALG/"prove_integer_identity.py","--candidate",candidate,"--irreducibility",folder/"factors.json",
        "--coordinates",coordinates,"--threads",threads,"--output",folder/"identity.json.gz")
    run(ALG/"bind_radius.py","--candidate",candidate,"--check",folder/"factors.json",
        "--coordinates",coordinates,"--target",folder/"target.json","--identity",folder/"identity.json.gz",
        "--output",folder/"binding.json")
    require(read(folder/"global.json")["global_optimality_certificate_chain_closed"] and
            read(folder/"chain.json")["verified"] and read(folder/"binding.json")["root_equals_global_optimal_radius_proved"],
            "全量重放没有证明最优半径等于指定实根")
    answer=export(folder/"delivery")
    for name in ("coefficients.json","ordered-centers.json","solution.json"):
        require((folder/"delivery"/name).read_bytes()==(HERE/name).read_bytes(),"交付文件逐字节复核失败："+name)
    print(answer)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-directory",type=Path,required=True,help="保存新生成报告的独立目录")
    p.add_argument("--threads",type=int,choices=range(1,7),default=5)
    args=p.parse_args()
    check(args.output_directory,args.threads)
