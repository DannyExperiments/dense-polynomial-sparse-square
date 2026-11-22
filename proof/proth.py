"""Generate or verify 128 explicit prime periods with a one-congruence certificate."""
import csv
from itertools import count
from math import isqrt, prod
from pathlib import Path


def certificates():
    for odd in count(7, 2):
        if odd >= 2**64:
            raise ValueError("Proth interval exhausted")
        candidate = odd * 2**64 + 1
        for witness in (3, 5, 7, 11, 13, 17, 19, 23):
            if pow(witness, (candidate - 1) // 2, candidate) == candidate - 1:
                yield odd, candidate, witness
                break


def check(rows):
    assert len(rows) == 128
    primes = set()
    for odd, prime, witness in rows:
        assert odd > 0 and odd % 2 == 1 and odd < 2**64
        assert prime == odd*2**64 + 1 and prime > 10**20
        assert prime not in primes
        primes.add(prime)
        assert pow(witness, (prime-1)//2, prime) == prime-1
        assert isqrt(prime) < 2**64 + 1
    return prod(primes)


if __name__ == "__main__":
    import sys
    from itertools import islice
    path = Path(__file__).with_name("proth-periods.csv")
    if "--make" in sys.argv:
        rows = list(islice(certificates(), 128))
        check(rows)
        with path.open("w", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(["odd", "prime", "witness"])
            writer.writerows(rows)
    with path.open() as handle:
        rows = [tuple(int(row[k]) for k in ("odd", "prime", "witness"))
                for row in csv.DictReader(handle)]
    N = check(rows)
    print("Distinct primes certified:", len(rows))
    print("First and last:", rows[0][1], rows[-1][1])
    print("Digits in period N:", len(str(N)))
    target = (14*len(rows) + 3)*N
    b = isqrt(target - 1) + 1
    degree = (N//400)*N + N + 4*sum(N//prime for _,prime,_ in rows) + b*(b+1)
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(20000)
    degree_string = str(degree)
    print("Exact number of degree digits:", len(degree_string))
    print("Degree prefix:", degree_string[:24])
    print("Degree suffix:", degree_string[-24:])
