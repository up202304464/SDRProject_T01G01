#!/usr/bin/env python3
"""
Week 05 - Parte C: fase dos pilotos, regressao linear por simbolo e
comparacao da constelacao sem correcao / fase constante / regressao linear.

Uso (a partir de week05/scripts):
    python3 comparacao_constelacao.py                       # vetores_sim.bin, 14 simbolos, BPSK
    python3 comparacao_constelacao.py ../data/vetores_s6.bin --nsym 20 --mod qpsk --tag s6
    python3 comparacao_constelacao.py --log ../data/run_sim.log   # valida contra o 'beta' do Bloessl

Estrutura assumida de cada trama no ficheiro (vetores de 64 pos-FFT, shift=Yes):
    LTS1, LTS2, SIGNAL, dado_1 ... dado_N      (+ possivel lixo ate ao fim das amostras do Sync Short)
"""
import argparse
import os
import re
import shutil
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

def _copiar(tmp, destino):
    """Grava primeiro em /tmp (Linux) e copia para o destino: contorna o Errno 22 do WSL em /mnt/c."""
    try:
        shutil.copyfile(tmp, destino)
    except OSError as e:
        alt = Path.cwd() / Path(destino).name
        shutil.copyfile(tmp, alt)
        print(f"AVISO: não consegui gravar em {destino} ({e}); ficou em {alt}")
    finally:
        Path(tmp).unlink(missing_ok=True)

def guardar_fig(destino):
    tmp = Path(tempfile.gettempdir()) / f"_tmp_{os.getpid()}_{Path(destino).name}"
    plt.savefig(tmp, dpi=150)
    plt.close()
    _copiar(tmp, destino)

# ----------------------------------------------------------------- argumentos
ap = argparse.ArgumentParser()
ap.add_argument("bin", nargs="?", default=None)
ap.add_argument("--nsym", type=int, default=14, help="simbolos de DADOS por trama (sem contar o SIGNAL)")
ap.add_argument("--mod", choices=["bpsk", "qpsk", "qam16", "qam64"], default="bpsk", help="modulacao dos dados")
ap.add_argument("--log", default=None, help="run_*.log do GNU Radio (opcional, valida o beta)")
ap.add_argument("--tag", default="sim")
ap.add_argument("--det", choices=["energia", "corr"], default="energia",
                help="deteção de tramas: rajadas de energia (default) ou correlação LTS1-LTS2")
ap.add_argument("--out", default="../figures")
ap.add_argument("--evm-max", type=float, default=0.5, help="EVM maximo do SIGNAL para aceitar a trama")
args = ap.parse_args()

# ----------------------------------------------------------------- dados
if args.bin:
    data_path = Path(args.bin)
else:
    data_path = Path("../data/vetores_sim.bin")
    if not data_path.exists():
        data_path = Path("../../week05/data/vetores_sim.bin")
v = np.fromfile(data_path, np.complex64).reshape(-1, 64).astype(np.complex128)

# ----------------------------------------------------------------- IEEE 802.11a
def _scrambler(estado, n):
    out = np.zeros(n, dtype=np.int8)
    for i in range(n):
        b = ((estado >> 6) ^ (estado >> 3)) & 1
        estado = ((estado << 1) | b) & 0x7F
        out[i] = b
    return out

POL = 1 - 2 * _scrambler(0x7F, 127)              # p_0 (SIGNAL), p_1 (1o dado), ...
PIL_K = np.array([-21, -7, 7, 21])
PIL_V = np.array([1, 1, 1, -1])
DADOS_K = np.array([k for k in range(-26, 27) if k != 0 and k not in PIL_K])
IDX_PIL = PIL_K + 32
IDX_DADOS = DADOS_K + 32

L = np.array([1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,
              0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1])
LFULL = np.zeros(64)
LFULL[np.arange(-26, 27) + 32] = L

# ----------------------------------------------------------------- 1. deteccao
usadas = [i for i in range(6, 59) if i != 32]
a_, b_ = v[:-1][:, usadas], v[1:][:, usadas]
c = np.abs(np.sum(np.conj(a_) * b_, axis=1)) / (
    np.linalg.norm(a_, axis=1) * np.linalg.norm(b_, axis=1) + 1e-12)
NS = args.nsym
ini_corr = np.flatnonzero((c > 0.9) & np.r_[True, c[:-1] <= 0.9])
if args.det == "corr":
    ini = ini_corr
else:
    # rajadas de energia: o ficheiro tem ~270 vetores por trama, a maioria ruído a seguir à trama
    e = np.sum(np.abs(v[:, usadas]) ** 2, axis=1)
    thr = 0.2 * np.percentile(e, 99)
    m = e > thr
    ini = np.flatnonzero(m & ~np.r_[False, m[:-1]])
    fim = np.flatnonzero(m & ~np.r_[m[1:], False])
    comp = fim - ini + 1
    vv, cc = np.unique(comp, return_counts=True)
    top = np.argsort(cc)[::-1][:8]
    print(f"Rajadas de energia: {len(ini)} | comprimentos (valor: nº):",
          {int(vv[i]): int(cc[i]) for i in top})
    print(f"  c[início da rajada] > 0.9 (LTS1≈LTS2): {np.mean(c[np.minimum(ini, len(c)-1)] > 0.9)*100:.0f}% das rajadas")
    longos = comp[comp >= 10]
    if len(longos):
        modo = int(np.bincount(longos).argmax())
        print(f"  comprimento mais comum das rajadas longas: {modo} vetores -> --nsym {modo - 3} ?")
    ini = ini[comp >= 3 + NS]
ini = ini[ini + 3 + NS <= len(v)]                 # LTS1+LTS2+SIGNAL+NS dados cabem no ficheiro

print(f"Vetores FFT: {len(v)} | tramas usadas: {len(ini)} (correlação sozinha daria {len(ini_corr)})")
d = np.diff(ini)
val, cnt = np.unique(d, return_counts=True)
top = np.argsort(cnt)[::-1][:8]
print("Distâncias entre tramas (valor: ocorrências):",
      {int(val[i]): int(cnt[i]) for i in top})
if len(ini) == 0:
    raise SystemExit("Nenhuma trama detetada.")

# ----------------------------------------------------------------- 2. equalizacao + pilotos (vetorizado)
idx = ini[:, None] + np.arange(3 + NS)[None, :]
V = v[idx]                                        # (F, 3+NS, 64)
H = (V[:, 0] + V[:, 1]) / 2 * LFULL               # H_k  (L=+-1  ->  1/L = L)
H[:, LFULL == 0] = 1.0
H[np.abs(H) < 1e-9] = 1.0
X = V[:, 2:] / H[:, None, :]                      # (F, NS+1, 64): n=0 SIGNAL, n>=1 dados

pol = POL[np.arange(NS + 1) % 127]
z = X[:, :, IDX_PIL] * PIL_V * pol[None, :, None]  # pilotos "desmodulados": todos com fase = erro de fase
soma = z.sum(-1)
ref = np.angle(soma)                               # fase constante (o que o Bloessl faz)
res = np.angle(z * np.exp(-1j * ref[..., None]))   # residuo por piloto (pequeno, sem saltos)
slope = res @ PIL_K / (PIL_K @ PIL_K)              # = polyfit (k simetrico => media de k = 0)
icpt = res.mean(-1)
ord_tot = ref + icpt                               # ordenada total da reta de fase (rad)

Xd = X[:, :, IDX_DADOS]
fase_reg = ref[..., None] + icpt[..., None] + slope[..., None] * DADOS_K
sem = Xd
const = Xd * np.exp(-1j * ref[..., None])
reg = Xd * np.exp(-1j * fase_reg)

# ----------------------------------------------------------------- 3. EVM e filtro de tramas
def fatia(x, mod):
    if mod == "bpsk":
        return np.sign(x.real) + 0j
    if mod == "qpsk":
        return (np.sign(x.real) + 1j * np.sign(x.imag)) / np.sqrt(2)
    lv = (np.array([-3, -1, 1, 3]) / np.sqrt(10) if mod == "qam16"
          else np.arange(-7, 8, 2) / np.sqrt(42))
    q = lambda u: lv[np.abs(u[..., None] - lv).argmin(-1)]
    return q(x.real) + 1j * q(x.imag)

def evm(x, mod):
    return np.mean(np.abs(x - fatia(x, mod)) ** 2, axis=-1)

evm_sig = evm(reg[:, 0], "bpsk")                   # SIGNAL e sempre BPSK
ok = evm_sig < args.evm_max
print(f"Tramas com SIGNAL aceitável (EVM < {args.evm_max}): {ok.sum()} / {len(ok)}")
if ok.sum() == 0:
    print("AVISO: nenhuma trama passa -> o problema é a deteção/alinhamento, não a correção de fase.")
    ok[:] = True

for nome, x in [("sem correção", sem), ("fase constante", const), ("regressão", reg)]:
    print(f"EVM dados [{nome:>14}]: {np.mean(evm(x[ok][:, 1:], args.mod)):.4f}")
h, _ = np.histogram(np.abs(reg[ok][:, 1:].real).ravel(), bins=np.linspace(0, 1.6, 17))
print("histograma de |I| dos dados (0 a 1.6, 16 classes):", h,
      "\n   BPSK: 1 pico em 1.0 | QPSK: 1 pico em 0.71 | 16-QAM: picos em 0.32 e 0.95")
print(f"|soma pilotos|/4 (dados): média {np.mean(np.abs(soma[ok][:, 1:]))/4:.3f}, "
      f"mínimo {np.min(np.abs(soma[ok][:, 1:]))/4:.3f}")
print(f"declive: rms {np.sqrt(np.mean(slope[ok]**2)):.4f} rad/subportadora | "
      f"ordenada: rms {np.sqrt(np.mean(icpt[ok]**2)):.4f} rad | ref: rms {np.sqrt(np.mean(ref[ok]**2)):.4f} rad")

# ----------------------------------------------------------------- 3b. diagnostico dos dados
print("\n--- Diagnóstico (tramas aceites) ---")
f0 = np.flatnonzero(ok)[0]
print(f"\nPrimeira trama aceite (vetor {ini[f0]}): subportadoras de dados 0..5 depois da correção (I, Q)")
for n in range(min(4, NS + 1)):
    print(f"  n={n}:", " ".join(f"({x.real:+.2f},{x.imag:+.2f})" for x in reg[f0, n, :6]))
print(f"|x| médio: SIGNAL {np.mean(np.abs(reg[ok][:, 0])):.2f} | dados {np.mean(np.abs(reg[ok][:, 1:])):.2f}"
      f" ; |Q| médio: SIGNAL {np.mean(np.abs(reg[ok][:, 0].imag)):.2f} | dados {np.mean(np.abs(reg[ok][:, 1:].imag)):.2f}")

e_n = [np.mean(evm(reg[ok][:, n], "bpsk" if n == 0 else args.mod)) for n in range(NS + 1)]
print("EVM por índice de símbolo (0 = SIGNAL):", np.round(e_n, 2))
print(f"|Σ pilotos|/4: SIGNAL {np.mean(np.abs(soma[ok][:, 0]))/4:.3f} | dados {np.mean(np.abs(soma[ok][:, 1:]))/4:.3f}")
for sh in range(-2, 3):                            # |Σ| nao depende da polaridade; cos(fase) sim
    pp = POL[(np.arange(NS + 1) + sh) % 127]
    zz = (X[ok][:, 1:, :][:, :, IDX_PIL] * PIL_V * pp[None, 1:, None]).sum(-1)
    print(f"  polaridade deslocada {sh:+d}: média cos(fase pilotos) dados = {np.mean(np.cos(np.angle(zz))):+.3f}")
print(f"vetores por trama (ficheiro/deteções): {len(v)/len(ini):.1f}  (esperado ~{3+NS} se só houver tramas)")

# ----------------------------------------------------------------- 4. validacao opcional com o beta do log
if args.log:
    txt = Path(args.log).read_text(errors="ignore")
    beta = np.array([float(x) for x in re.findall(r"Key: beta\s+Value: (\S+)", txt)])
    print(f"beta no log: {len(beta)} valores | tramas aqui: {len(ref)}")
    m = min(len(beta), len(ref))
    cs = [np.corrcoef(beta[:m], ref[:m, n])[0, 1] for n in range(NS + 1)]
    print("  corr(beta_log, ref) por índice de símbolo (0 = SIGNAL):", np.round(cs, 2))
    nb = int(np.argmax(np.abs(cs)))
    ganho = np.polyfit(ref[:m, nb], beta[:m], 1)[0]
    print(f"  melhor: n={nb} (corr {cs[nb]:+.3f}); beta ≈ {ganho:+.2f} × ref")
    if len(beta) != len(ref):
        print("  (nº diferente: se a correlação for ~0 pode ser só desalinhamento de tramas)")

# ----------------------------------------------------------------- 5. figuras
out = Path(args.out)
out.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(0)
Fok = np.flatnonzero(ok)

# 5a. fase dos pilotos vs k
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
amostras = [(rng.choice(Fok), rng.integers(1, NS + 1)) for _ in range(12)]
for f, n in amostras:
    ax[0].plot(PIL_K, np.angle(z[f, n]), "o-", alpha=0.6)
ax[0].set(title="Fase dos 4 pilotos (12 símbolos aleatórios)", xlabel="k", ylabel="fase [rad]")
for f, n in amostras[:4]:
    p = ax[1].plot(PIL_K, res[f, n], "o", label=f"trama {f}, sím. {n}")
    kk = np.linspace(-22, 22, 2)
    ax[1].plot(kk, icpt[f, n] + slope[f, n] * kk, "-", color=p[0].get_color(), alpha=0.7)
ax[1].set(title="Resíduo (fase − fase média) + reta de regressão", xlabel="k", ylabel="rad")
ax[1].legend(fontsize=8)
for a in ax:
    a.grid(True, ls="--", alpha=0.5)
plt.tight_layout()
guardar_fig(out / f"week05_fase_pilotos_vs_k_{args.tag}.png")

# 5b. declive, ordenada e |soma| ao longo da trama
fig, ax = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
n_ax = np.arange(NS + 1)
for a, y, t in zip(ax, [ord_tot, slope, np.abs(soma) / 4],
                   ["Ordenada total b = ref + resíduo [rad]", "Declive a [rad/subportadora]", "|Σ pilotos| / 4"]):
    for f in Fok[:150]:
        a.plot(n_ax, y[f], color="gray", alpha=0.12)
    a.errorbar(n_ax, y[ok].mean(0), y[ok].std(0), color="C3", marker="o", capsize=3, label="média ± desvio")
    a.set_ylabel(t)
    a.grid(True, ls="--", alpha=0.5)
ax[0].set_title("Parâmetros da regressão por símbolo (0 = SIGNAL)")
ax[0].legend()
ax[-1].set_xlabel("índice do símbolo na trama")
plt.tight_layout()
guardar_fig(out / f"week05_regressao_por_simbolo_{args.tag}.png")

# 5c. constelacoes
def hist(a, x, titulo):
    x = x.ravel()
    h = a.hist2d(x.real, x.imag, bins=160, range=[[-2, 2], [-2, 2]], norm=LogNorm(), cmin=1, cmap="viridis")
    a.set(title=titulo, xlabel="I", ylabel="Q", aspect="equal")
    a.axhline(0, color="k", lw=0.6)
    a.axvline(0, color="k", lw=0.6)

fig, axes = plt.subplots(2, 3, figsize=(14, 9))
nomes = ["Sem correção (só LTS)", "Fase constante (pilotos)", "Regressão linear (pilotos)"]
for j, (x, nm) in enumerate(zip([sem, const, reg], nomes)):
    hist(axes[0, j], x[ok][:, 1:], f"DADOS – {nm}\nEVM = {np.mean(evm(x[ok][:, 1:], args.mod)):.3f}")
    hist(axes[1, j], x[ok][:, 0], f"SIGNAL – {nm}\nEVM = {np.mean(evm(x[ok][:, 0], 'bpsk')):.3f}")
plt.suptitle(f"Week 05 – Parte C ({args.tag}, {ok.sum()} tramas, {NS} símbolos de dados)", fontweight="bold")
plt.tight_layout()
guardar_fig(out / f"week05_comparacao_correcao_pilotos_{args.tag}.png")

# 5d. CSV por trama (para ligar a Semana 3 / D = 336)
tab = np.column_stack([np.arange(len(ini)), ini, evm_sig,
                       np.abs(soma[:, 0]) / 4, np.abs(soma[:, 1:]).mean(1) / 4,
                       ref[:, 0], np.angle(np.exp(1j * ord_tot[:, 1:]).mean(1)), slope[:, 1:].mean(1)])
tmp_csv = Path(tempfile.gettempdir()) / f"_tmp_{os.getpid()}.csv"
np.savetxt(tmp_csv, tab, delimiter=",",
           header="trama,inicio,evm_signal,abs_soma_signal,abs_soma_dados,ref_signal,b_med,a_med", comments="")
_copiar(tmp_csv, out / f"pilotos_por_trama_{args.tag}.csv")
print(f"Figuras e CSV gravados em: {out.resolve()}")