import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = "/mnt/c/Users/luisa/SDRProject_T01G01/week01/data/week01_signal.bin"
out = "/mnt/c/Users/luisa/SDRProject_T01G01/week02/figures"
s0 = np.fromfile(path, dtype=np.complex64)[:1560192]
W, LAG, THR = 48, 16, 0.56
print("amostras lidas:", len(s0))

def msum(x, L):
    c = np.cumsum(x, dtype=np.complex128)
    o = c.copy()
    o[L:] -= c[:-L]
    return o

def corr(noise, seed=0):
    rng = np.random.default_rng(seed)
    N = len(s0)
    n = noise/np.sqrt(2)*(rng.standard_normal(N) + 1j*rng.standard_normal(N))
    s = s0.astype(np.complex128) + n
    prod = np.zeros(N, dtype=np.complex128)
    prod[LAG:] = s[LAG:] * np.conj(s[:-LAG])
    a = np.abs(msum(prod, W))
    p = msum(np.abs(s)**2, W + LAG).real
    return a / np.maximum(p, 1e-12)

c0 = corr(0)
i0 = int(np.argmax(c0 > THR))
print("primeira frame no indice", i0, "| pico de c[n] sem ruido:", round(float(c0.max()), 3))

for noise in (0, 3.5, 4):
    c = corr(noise)[i0-100:i0+300]
    plt.figure(figsize=(6, 3.5))
    plt.plot(np.arange(-100, 300), c)
    plt.axhline(THR, ls="--", c="r", label=f"threshold = {THR}")
    plt.ylim(0, 1.2)
    plt.xlabel("Índice da amostra")
    plt.ylabel("c[n]")
    plt.title(f"noise = {noise}")
    plt.legend()
    plt.tight_layout()
    f = f"{out}/c_noise{noise}.png"
    plt.savefig(f, dpi=200)
    plt.close()
    print("guardado:", f)