"""Export the actual lacunary numerators used by the compact witness."""
import csv
import json
from functools import lru_cache
from math import comb
from pathlib import Path
from local_sections import context, section

here=Path(__file__).parent
with (here/'proth-periods.csv').open() as stream:
    periods=[int(r['prime']) for r in csv.DictReader(stream)]
records=[]
for p in periods:
    @lru_cache(None)
    def rows(r):
        return section(p,r)

    def coefficient(k):
        if k<0 or k>=9*p:return 0
        layer,r=divmod(k,p)
        return rows(r)[layer]

    e=context(p)[0][0]
    positions={layer*p+cut+i for layer in range(9)
               for cut in e.piece_boundaries for i in range(4)}
    mask=[]
    for k in sorted(positions):
        value=sum((-1)**i*comb(4,i)*coefficient(k-i) for i in range(5))
        if value:
            mask.append([k,value])
    records.append({'prime':p,'terms':mask})
with (here/'compact-local-masks.json').open('w') as out:
    json.dump(records,out,separators=(',',':'))
    out.write('\n')
print('Exported',len(records),'masks;',sum(len(r['terms']) for r in records),'total terms')
