"""在数学证书齐备后导出题面所需的整系数数组、有序圆心及官方协议指纹。"""

import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import gcd
from pathlib import Path
import subprocess
import sys

from chart import HERE,read
from context import ROOT
from sample_plane import digest
from verify import require


def exact_file(path,content):
    raw=content.encode("utf-8")
    if path.exists():
        require(path.read_bytes()==raw,"现存交付件与重算结果不一致："+str(path))
    else:
        with path.open("xb") as stream:
            stream.write(raw)


def export(folder):
    source=ROOT/"n15/proofs/global-final"
    candidate=HERE/"reconstruction-p256.json.gz"
    certificate=HERE/"integer-candidate-check.json"
    identity=HERE/"integer-identity.json.gz"
    binding=HERE/"radius-binding.json"
    coords=HERE/"scalars-p256.json.gz"
    s=read(candidate)
    p=s["candidate_primitive_coefficients_ascending"]
    irreducible=read(certificate)
    proven=read(identity)
    bound=read(binding)
    global_report=read(source/"master-20261118.json")
    chain=read(source/"axis-chain-20261118.json")
    require(global_report["global_optimality_certificate_chain_closed"] and global_report["unresolved_leaves"]==0
            and global_report["enumeration_total"]==global_report["metric_excluded"]+global_report["noncandidate_cases"]+1,
            "全局最优性证书未闭合")
    require(irreducible["source_sha256"]==proven["candidate_sha256"]==bound["candidate_sha256"]==digest(candidate)
            and irreducible["irreducibility"]["irreducibility_proved"] and irreducible["root_count"]["unique_root_proved"],
            "单变量多项式的不可约性或根数不符")
    require(proven["irreducibility_sha256"]==digest(certificate) and
            proven["coordinates_sha256"]==bound["coordinate_sha256"]==digest(coords)
            and sum(proven["integer_point_bits"])==proven["sum_lower_bits"]>proven["upper_bits"]
            and proven["all_six_rational_function_identities_proved"],"整数恒等式证书不符")
    require(bound["identity_sha256"]==digest(identity) and bound["root_box_sha256"]==digest(HERE/"scalar-target.json")
            and bound["root_equals_global_optimal_radius_proved"] and bound["real_embedding_in_unique_source_box"],
            "整系数实根没有连接到几何根")
    require(len(p)==2870 and p[-1]>0 and gcd(*p)==1,"导出的系数未规范化")
    folder.mkdir(parents=True,exist_ok=True)
    file=folder/"coefficients.json"
    exact_file(file,"[\n"+",\n".join(str(i) for i in p)+"\n]\n")
    decimals=global_report["decimal_intervals"]
    intervals=[decimals["r_lower"],decimals["r_upper"]]
    require(list(map(Q,intervals))==list(map(Q,irreducible["root_count"]["endpoints"])),"全局半径与单变量隔离端点不同")
    core=source/"replay/core"
    root=read(core/"root15.json")
    require(len(root["xnum"])==148 and chain["verified"] and chain["source_replayed"]
            and chain["variables"]==25 and chain["global_root_sha256"]==digest(core/"root15.json"),
            "圆心来源有理盒或 148 元至 25 元根链不符")
    centers=[{"index":i+1,"coordinates":[20+2*i,21+2*i],"box_center_numerators":
              [str(root["xnum"][20+2*i]),str(root["xnum"][21+2*i])]} for i in range(15)]
    centerfile=folder/"ordered-centers.json"
    centerdoc={"coordinate_order":"从零编号的 148 元唯一根坐标","root_data":"../proofs/global-final/replay/core/root15.json",
               "system":"../proofs/global-final/replay/core/system15.py","denominator":str(root["Qx"]),
               "box_half_width":"1/"+str(10**80),"radius_squared_coordinate":76,"centers":centers}
    exact_file(centerfile,json.dumps(centerdoc,ensure_ascii=False,indent=2)+"\n")
    upstream="disk-covering-benchmark-122ba78c748303f4d590d20c7df16630dba51fd2"
    official=ROOT/"references/upstream"/upstream/"tools/diskcover_answer.py"
    digests=[]
    for endpoint in intervals:
        out=subprocess.check_output([sys.executable,"-X","utf8","-B",str(official),"--n","15","--coefficients",
                                     str(file),"--radius",endpoint])
        report=json.loads(out)
        require(report["answer"]["coefficients"]==p and report["answer"]["n"]==15,"官方脚本改变了规范系数")
        encoded=json.dumps(report["answer"],sort_keys=True,separators=(",",":"),ensure_ascii=True)
        require(encoded==report["canonical_json"] and hashlib.sha256(encoded.encode()).hexdigest()==report["sha256"],
                "独立重算的协议指纹不符")
        digests.append(report["sha256"])
    require(digests[0]==digests[1],"隔离区间跨过协议舍入边界")
    local_record={"n":15,"degree":len(p)-1,"coefficients_order":"ascending","coefficients_sha256":digest(file),
                  "interval":intervals,"endpoint_signs":irreducible["root_count"]["endpoint_signs"],
                  "ordered_centers_sha256":digest(centerfile),"root_sha256":digest(core/"root15.json"),
                  "system_sha256":digest(core/"system15.py"),"global_replay_sha256":digest(source/"master-20261118.json"),
                  "axis_chain_sha256":digest(source/"axis-chain-20261118.json"),"irreducibility_sha256":digest(certificate),
                  "integer_identity_sha256":digest(identity),"binding_sha256":digest(binding),
                  "protocol":"diskcover-minpoly-v1","protocol_script_sha256":digest(official),
                  "rounded_radius":report["answer"]["radius"],"answer_sha256":digests[0],
                  "reference_status":read(ROOT/"n15/problem-metadata.zh.json")["reference_status"],"matches_reference":None}
    exact_file(folder/"solution.json",json.dumps(local_record,ensure_ascii=False,indent=2)+"\n")
    return digests[0]


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--folder",type=Path,default=ROOT/"n15/solution")
    args=p.parse_args()
    from windows_compute_limits import apply_limits
    _,close=apply_limits()
    try:
        print(export(args.folder))
    finally:
        close()
