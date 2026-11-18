"""重算根盒、消元恒等式、实际纤维以及四个完整的模截面，逐字段核对归档。"""

from fractions import Fraction as Q
from chart import HERE, read
from verify_chart import prove
from verify_elimination import verify_projection, characteristic_zero_slice
from binary_projection import verify
from target_fiber import analyze
from radius_slice import run
from windows_compute_limits import apply_limits


def compare(name, computed):
    if read(HERE / name) != computed:
        raise ValueError("重算结果不同：" + name)
    print(name, flush=True)


if __name__ == "__main__":
    _, close = apply_limits()
    try:
        compare("target.json", prove())
        projection = verify_projection()
        projection["slice"] = characteristic_zero_slice()
        compare("projection.json", projection)
        compare("binary.json", verify())
        compare("target-fiber.json", analyze())
        for name, fixed, vary_radius in (("slice-9-16", Q(9, 16), False),
                                         ("slice-2-5", Q(2, 5), False),
                                         ("slice-a3", Q(3), True),
                                         ("slice-a-third", Q(1,3), True)):
            compare(name+".json.gz", run(fixed, vary_radius))
    finally:
        close()
