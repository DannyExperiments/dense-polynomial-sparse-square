"""Verify the parameters and evaluate the paper's support bounds with integers."""
import csv
from itertools import islice
import json
from math import isqrt
from pathlib import Path
import sys
from proth import certificates, check


def certificate():
    with Path(__file__).with_name("proth-periods.csv").open() as source:
        rows = [tuple(int(row[key]) for key in ("odd", "prime", "witness"))
                for row in csv.DictReader(source)]
    t = len(rows)
    assert rows == list(islice(certificates(), 128))
    period = check(rows)
    base = 1 + isqrt((14*t + 3)*period - 1)
    repetitions = period // 400
    filler_degree = period + 4*sum(period//p for _, p, _ in rows) + base*(base+1)
    degree = repetitions*period + filler_degree
    local_quotient_bound = (18*t + 1)*period*7**t // 10**t
    repaired_quotient_bound = (5*local_quotient_bound + 18*base*10**(6*t)
                               + 3*(8*base + 2)*9**t)
    support_bound = 5*(2*repetitions - 1) + 3*repaired_quotient_bound
    assert 400*repaired_quotient_bound < period
    assert 9*t*period + 2 <= filler_degree < repetitions*period
    assert support_bound**2 < degree
    assert (1000*support_bound)**2 < 501**2*degree
    assert len(str(degree)) == 5798
    return {
        "prime_count": t, "period": str(period), "auxiliary_base": str(base),
        "repetitions": str(repetitions), "degree": str(degree),
        "square_weight_upper_bound": str(support_bound),
        "degree_decimal_digits": len(str(degree)),
        "square_weight_over_sqrt_degree_strict_upper_bound": "501/1000",
    }


if __name__ == "__main__":
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(20000)
    data = certificate()
    recorded = Path(__file__).with_name("certificate.json")
    if sys.argv[1:] == ["--write"]:
        recorded.write_text(json.dumps(data, indent=2) + "\n")
    else:
        assert json.loads(recorded.read_text()) == data
    print("Verified all 128 periods, the 5798-digit degree, and the strict 501/1000 bound.")
