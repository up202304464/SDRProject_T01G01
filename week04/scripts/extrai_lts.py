# extrai_lts.py
import re, numpy as np
t = open("/mnt/c/Users/luisa/gr-ieee802-11/lib/sync_long.cc").read()
m = re.search(r"LONG\s*=\s*\{(.*?)\};", t, re.S)
p = re.findall(r"gr_complex\(\s*([-+0-9.eE]+)f?\s*,\s*([-+0-9.eE]+)f?\s*\)", m.group(1))
taps = np.array([complex(float(a), float(b)) for a, b in p], dtype=np.complex64)
print(len(taps), taps[0], taps[63])
np.save("/mnt/c/Users/luisa/SDRProject_T01G01/week04/data/lts_taps.npy", taps)