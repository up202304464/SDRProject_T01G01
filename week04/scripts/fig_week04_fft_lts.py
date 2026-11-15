# fig_week04_fft_lts.py   uso: python3 fig_week04_fft_lts.py [fft_lts.bin] [--out fig.png]
# Lê a saída do FFT (vetores de 64, shift sim) e mostra a FFT da 1.ª cópia da LTS: 52 subportadoras ativas e 12 nulas.
import argparse, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/mnt/c/Users/luisa/SDRProject_T01G01/week04"
ap = argparse.ArgumentParser()
ap.add_argument("ficheiro", nargs="?", default=f"{BASE}/data/fft_lts.bin")
ap.add_argument("--out", default=f"{BASE}/figures/fft_lts.png")
ap.add_argument("--vetor", type=int, default=None, help="índice do vetor de 64 a mostrar (por defeito: 1.ª cópia da LTS encontrada)")
a = ap.parse_args()

X = np.fromfile(a.ficheiro, dtype=np.complex64)
n = len(X) // 64
print(f"{len(X)} amostras = {n} vetores de 64 (resto: {len(X) % 64})")
if n < 2: raise SystemExit("poucos vetores: confirma que o File Sink tem vlen 64 e que a trama passou no Sync Long")
X = X[:n * 64].reshape(n, 64)
k = np.arange(-32, 32)                      # índice da subportadora (com shift: DC no meio)

def limiar(v, rel=0.4):
    m = np.abs(v)
    return rel * np.median(np.sort(m)[-48:])        # 40 % da mediana das 48 maiores (as subportadoras ativas)

def ativas(v):
    return np.abs(v) > limiar(v)

def ativas_n(v):
    return int(ativas(v).sum())

# Como distinguir a LTS de símbolos de dados? Os símbolos BPSK/QPSK também têm módulo constante nas 52
# subportadoras, por isso o módulo NÃO serve. O que distingue a LTS é que a 2.ª cópia é IGUAL à 1.ª (em complexo,
# a menos de uma rotação de fase comum). Dados consecutivos são diferentes (rho ~ 1/sqrt(52) ~ 0,14).
num = np.abs(np.sum(X[:-1] * np.conj(X[1:]), axis=1))
rho = num / (np.linalg.norm(X[:-1], axis=1) * np.linalg.norm(X[1:], axis=1) + 1e-30)   # 1 = vetor i igual ao i+1
pares = np.flatnonzero(rho > 0.9)
print(f"pares de vetores consecutivos iguais (rho > 0,9, = 1.ª e 2.ª cópia da LTS): {len(pares)} de {n - 1}")
print("primeiros 8 vetores -> nº ativas | rho com o vetor seguinte:")
for i in range(min(8, n - 1)):
    print(f"  vetor {i}: {ativas_n(X[i])} ativas | rho = {rho[i]:.3f}")
if len(pares):
    print("primeiros pares (índice do vetor):", pares[:10].tolist())
    d = np.diff(pares); vals, cont = np.unique(d, return_counts=True)
    print("distância entre pares consecutivos (em vetores) -> nº de vezes:", dict(zip(vals.tolist()[:5], cont.tolist()[:5])))

if a.vetor is None:
    a.vetor = int(pares[0]) if len(pares) else 0
    print("a usar o vetor", a.vetor, "(1.ª cópia da LTS)" if len(pares) else "(nenhum par igual encontrado: a usar 0)")
v = X[a.vetor]; v2 = X[a.vetor + 1]
at = ativas(v)
print(f"vetor {a.vetor}: {at.sum()} subportadoras ativas e {64 - at.sum()} nulas (limiar = 40 % da mediana das 48 maiores = {limiar(v):.1f})")
m = np.abs(v); print(f"nulas: máx {m[~at].max():.1f} ({20*np.log10(m[~at].max()/np.median(m[at])):.1f} dB rel. à mediana das ativas) | ativas: mín {m[at].min():.1f}, mediana {np.median(m[at]):.1f}, máx {m[at].max():.1f}")
print("nulas nas posições:", k[~at].tolist())
print(f"2.ª cópia (vetor {a.vetor + 1}) vs 1.ª: rho = {rho[a.vetor]:.3f} (1 = iguais)")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
ax1.stem(k, np.abs(v), basefmt=" ")
ax1.axhline(limiar(v), c="r", ls="--", lw=1, label="limiar de 'ativa'"); ax1.legend(loc="lower center")
ax1.set_ylabel("|X[k]|"); ax1.set_title(f"FFT da cópia {a.vetor + 1} da LTS: {at.sum()} ativas, {64 - at.sum()} nulas")
ax2.stem(k, np.angle(v) * (at), basefmt=" ")
ax2.set_ylabel("fase (rad)"); ax2.set_xlabel("índice da subportadora k (DC em k = 0)")
for ax in (ax1, ax2): ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(a.out, dpi=150); print("guardado em", a.out)