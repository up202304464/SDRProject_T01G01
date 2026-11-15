# verifica_amostras.py   uso: python3 verifica_amostras.py caminho/*.bin
import numpy as np, sys
for f in sys.argv[1:]:
    x = np.fromfile(f, dtype=np.complex64)
    print(f"{f}: {len(x)} amostras = {len(x)/20e6*1e3:.1f} ms | potência média {np.mean(np.abs(x)**2):.3g} | "
          f"|x|max {np.abs(x).max():.3g} | NaN: {np.isnan(x).any()}")