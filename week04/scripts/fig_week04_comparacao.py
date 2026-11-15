import argparse, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/mnt/c/Users/luisa/SDRProject_T01G01/week04"
ap = argparse.ArgumentParser()
ap.add_argument("ficheiro")
ap.add_argument("--taps", default=f"{BASE}/data/lts_taps.npy")
ap.add_argument("--out", default=f"{BASE}/figures/comparacao_c_vs_mf.png")
ap.add_argument("--thr", type=float, default=0.56)
ap.add_argument("--plateau", type=int, default=2, help="min_plateau do Sync Short (exige min_plateau + 1 amostras seguidas acima do limiar)")
ap.add_argument("--minlen", type=int, default=100, help="duração mínima do plateau (a STS dá ~130 amostras)")
ap.add_argument("--frame", type=int, default=0, help="qual das tramas com LTS confirmada usar (0 = primeira)")
ap.add_argument("--pmr", type=float, default=5.0, help="os 2 picos a 64 têm de ser >= pmr x mediana da saída do filtro")
ap.add_argument("--todas", action="store_true", help="não exigir LTS confirmada (usa qualquer plateau longo)")
ap.add_argument("--semcfo", action="store_true", help="NÃO corrigir o desvio de frequência grosseiro antes do filtro adaptado (por defeito corrige, como o Sync Short)")
ap.add_argument("--win", type=int, default=48)
ap.add_argument("--start", type=int, default=0, help="primeira amostra a ler")
ap.add_argument("--max", type=int, default=2_000_000, help="nº máx. de amostras a ler (2e6 = 100 ms)")
a = ap.parse_args()

x = np.fromfile(a.ficheiro, dtype=np.complex64, count=a.max, offset=a.start * 8)
print(f"lidas {len(x)} amostras (a partir da amostra {a.start})")
if len(x) < 1000: raise SystemExit("poucas amostras lidas: verifica --start e o ficheiro")
taps = np.load(a.taps)                      # LONG: LTS conjugada e invertida (como no sync_long.cc)

# --- c[n], eqs. (1)-(3): janela 48 para a[n], janela 48+16 para p[n] (como no fluxograma) ---
def msum(v, n):                              # soma deslizante: r[i] = sum(v[i:i+n])
    c = np.concatenate(([0], np.cumsum(v, dtype=np.complex128 if np.iscomplexobj(v) else np.float64)))
    return c[n:] - c[:-n]
prod = np.conj(x[:-16]) * x[16:]
A = msum(prod, a.win)
P = msum(np.abs(x) ** 2, a.win + 16)
L = min(len(A), len(P))
c = np.abs(A[:L]) / np.maximum(P[:L], 1e-30)
c = np.concatenate((np.zeros(a.win + 15), c))      # janela CAUSAL, como as médias móveis do GNU Radio
del prod, A, P

# --- plateaus: sequências de amostras acima do limiar ---
acima = np.concatenate(([False], c > a.thr, [False]))
d = np.diff(acima.astype(np.int8))
ini, fim = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
dur = fim - ini
n_det = int(np.sum(dur >= a.plateau + 1))          # o que o Sync Short (min_plateau) consideraria "trama"
validos = np.flatnonzero(dur >= a.minlen)          # plateaus com duração de STS
print(f"plateaus com >= {a.plateau + 1} amostras seguidas acima do limiar (min_plateau = {a.plateau}): {n_det}")
print(f"plateaus com >= {a.minlen} amostras (duração típica da STS): {len(validos)}")

def cfo_rad(n0):                                    # como o Sync Short: arg(a[n]) / 16, medido dentro da STS
    return np.angle(np.sum(np.conj(x[n0:n0 + 100]) * x[n0 + 16:n0 + 116])) / 16

def alinha(n0):
    seg = x[n0:n0 + 320]
    if not a.semcfo:                                # correção grosseira (Sync Short) ANTES do filtro adaptado
        seg = seg * np.exp(-1j * cfo_rad(n0) * np.arange(len(seg)))
    mf = np.abs(np.convolve(seg, taps)[:len(seg)])      # FIR causal: pico no FIM da janela (n+63 da eq. 6)
    picos = []
    for i in np.argsort(mf)[::-1]:
        if all(abs(i - j) >= 10 for j in picos): picos.append(int(i))
        if len(picos) == 3: break
    return mf, sorted(picos)

def tem_par64(p, mf):                                   # par de picos a 64 (±1) E proeminentes (não ruído)
    lim = a.pmr * np.median(mf)
    return any(abs((q - r) - 64) <= 1 and mf[q] >= lim and mf[r] >= lim for r in p for q in p if q > r)

# quantos plateaus longos têm LTS confirmada (par de picos a 64, ambos >= pmr x mediana)?
boas = [i for i in validos[:300] if tem_par64(*reversed(alinha(int(ini[i]))))]
print(f"{len(boas)} de {min(len(validos), 300)} plateaus longos têm LTS confirmada (par a 64, picos >= {a.pmr}x a mediana; "
      f"{'SEM' if a.semcfo else 'COM'} correção de frequência)")
if boas: print("posições (amostra do ficheiro):", [a.start + int(ini[i]) for i in boas[:10]])
print(f"RESUMO,{a.ficheiro.split('/')[-1]},{len(x)},{n_det},{len(validos)},{min(len(validos), 300)},{len(boas)},{'semcfo' if a.semcfo else 'cfo'}")
lista = list(validos) if a.todas else boas
if len(lista) <= a.frame:
    raise SystemExit("nenhuma trama com LTS confirmada neste troço; tenta --start com outro valor (ou --todas para ver os plateaus sem LTS)")
iq = lista[a.frame]
n0 = int(ini[iq])
mf, picos = alinha(n0)
par = [(r, q) for r in picos for q in picos if abs(q - r - 64) <= 1]
n_dados_eq7 = picos[-1] + 1                              # eq. 7 à letra: último dos 3 maiores picos (+1: o filtro dá o pico no fim da janela)
n_dados = (par[-1][1] + 1) if par else n_dados_eq7       # 2.º pico do ÚLTIMO par a 64 (os 2 picos das cópias da LTS; o 1.º pico, a meia altura, vem do GI2)
print(f"trama {a.frame}: plateau de {int(dur[iq])} amostras, deteção em n0 = {n0} (amostra {a.start + n0} do ficheiro)")
print(f"desvio de frequência estimado na STS: {cfo_rad(n0) * 20e6 / (2 * np.pi) / 1e3:.1f} kHz "
      f"({'NÃO corrigido' if a.semcfo else 'corrigido antes do filtro'})")
print("3 maiores picos (amostras desde a deteção):", picos, "| início dos dados (par a 64):", n_dados,
      "| eq. 7 à letra (último pico):", n_dados_eq7)
print("espaçamentos:", np.diff(picos), "| LTS confirmada:", tem_par64(picos, mf))

# --- figura ---
ant = 60
i0 = max(n0 - ant, 0)
cc = c[i0: n0 + 380]
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
ax1.plot(np.arange(i0 - n0, i0 - n0 + len(cc)), cc, lw=1.2)
ax1.axhline(a.thr, ls="--", c="r", lw=1, label=f"limiar = {a.thr}")
ax1.axvline(0, c="k", lw=0.8, ls=":", label="deteção (início do plateau)")
ax1.set_ylabel("c[n]"); ax1.set_title("Deteção: autocorrelação normalizada (lag 16)"); ax1.legend(loc="upper right")
ax2.plot(np.arange(len(mf)), mf, lw=1.2)
ax2.plot(picos, mf[picos], "ro", label="3 maiores picos")
ax2.axvline(n_dados, c="g", ls="--", label=f"início dos dados (n = {n_dados})")
if n_dados_eq7 != n_dados:
    ax2.axvline(n_dados_eq7, c="orange", ls=":", label=f"eq. 7 à letra (n = {n_dados_eq7})")
ax2.set_ylabel("|saída do filtro adaptado|")
ax2.set_xlabel("amostras desde a deteção (20 Msps: 20 amostras = 1 µs)")
ax2.set_title("Alinhamento: correlação com a LTS (64 coeficientes)"); ax2.legend(loc="upper left")
fig.tight_layout(); fig.savefig(a.out, dpi=150)
print("guardado em", a.out)