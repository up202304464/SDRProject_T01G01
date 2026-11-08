#!/usr/bin/env python3
"""
Week 03 - simulação em Python (numpy) do Sync Short + Sync Long sobre a gravação.

NÃO é o GNU Radio: reimplementa a lógica de sync_short.cc e sync_long.cc
(gr-ieee802-11, maint-3.10) para medir o efeito do desvio de frequência,
do atraso D e do sync_length sem precisar do GNU Radio.

Métrica de "descodificação": o campo SIGNAL (1.º símbolo depois da LTS) é
equalizado com a estimativa de canal da LTS, desmodulado (BPSK), desentrelaçado
e descodificado com Viterbi (K=7, r=1/2). Conta como bom se a paridade estiver
certa, a taxa for válida e a cauda for zero. Não descodifica o payload nem
verifica o CRC (isso é o que o Decode MAC faz no GNU Radio).

Uso:
  python3 sim_week03.py --input-file ../../week01/data/week01_signal.bin --varrimento freq
  (varrimentos: base, freq, delay, synclength, threshold)
"""
import argparse
import csv
import zlib
from pathlib import Path

import numpy as np

FS = 20e6

# ---------------------------------------------------------------------------
# Sync Long: coeficientes do filtro adaptado (linhas 257-290 do sync_long.cc)
LONG = np.array([
    -0.0455-1.0679j, 0.3528-0.9865j, 0.8594+0.7348j, 0.1874+0.2475j,
    0.5309-0.7784j, -1.0218-0.4897j, -0.3401-0.9423j, 0.8657-0.2298j,
    0.4734+0.0362j, 0.0088-1.0207j, -1.2142-0.4205j, 0.2172-0.5195j,
    0.5207-0.1326j, -0.1995+1.4259j, 1.0583-0.0363j, 0.5547-0.5547j,
    0.3277+0.8728j, -0.5077+0.3488j, -1.1650+0.5789j, 0.7297+0.8197j,
    0.6173+0.1253j, -0.5353+0.7214j, -0.5011-0.1935j, -0.3110-1.3392j,
    -1.0818-0.1470j, -1.1300-0.1820j, 0.6663-0.6571j, -0.0249+0.4773j,
    -0.8155+1.0218j, 0.8140+0.9396j, 0.1090+0.8662j, -1.3868-0.0000j,
    0.1090-0.8662j, 0.8140-0.9396j, -0.8155-1.0218j, -0.0249-0.4773j,
    0.6663+0.6571j, -1.1300+0.1820j, -1.0818+0.1470j, -0.3110+1.3392j,
    -0.5011+0.1935j, -0.5353-0.7214j, 0.6173-0.1253j, 0.7297-0.8197j,
    -1.1650-0.5789j, -0.5077-0.3488j, 0.3277-0.8728j, 0.5547+0.5547j,
    1.0583+0.0363j, -0.1995-1.4259j, 0.5207+0.1326j, 0.2172+0.5195j,
    -1.2142+0.4205j, 0.0088+1.0207j, 0.4734-0.0362j, 0.8657+0.2298j,
    -0.3401+0.9423j, -1.0218+0.4897j, 0.5309+0.7784j, 0.1874-0.2475j,
    0.8594-0.7348j, 0.3528+0.9865j, -0.0455+1.0679j, 1.3868-0.0000j,
], dtype=np.complex128)
LTS_T = np.conj(LONG[::-1])          # LTS no tempo: LONG[k] = conj(L[63-k])
LTS_F = np.fft.fft(LTS_T)            # LTS na frequência (índice k mod 64)

# subportadoras
PILOTOS = {-21, -7, 7, 21}
DADOS = [k for k in range(-26, 27) if k != 0 and k not in PILOTOS]   # 48, ordem crescente
USADAS = [k for k in range(-26, 27) if k != 0]

# SIGNAL (BPSK, N_CBPS = 48): entrelaçador
_k = np.arange(48)
INTERL = (48 // 16) * (_k % 16) + _k // 16          # bit codificado k -> posição INTERL[k]
TAXAS = {(1, 1, 0, 1): 6, (1, 1, 1, 1): 9, (0, 1, 0, 1): 12, (0, 1, 1, 1): 18,
         (1, 0, 0, 1): 24, (1, 0, 1, 1): 36, (0, 0, 0, 1): 48, (0, 0, 1, 1): 54}


def _paridade(x):
    return bin(x).count("1") & 1


# Viterbi K=7, g0=133, g1=171 (octal); bit 6 do registo = entrada actual (D^0)
_G0, _G1 = 0o133, 0o171
_SAIDA = np.zeros((64, 2, 2), dtype=np.int8)
for _s in range(64):
    for _b in range(2):
        _r = (_b << 6) | _s
        _SAIDA[_s, _b] = (_paridade(_r & _G0), _paridade(_r & _G1))


# Para cada estado seguinte ns: bit de entrada b = ns >> 5 e dois antecessores s = ((ns & 31) << 1) | x
_NS = np.arange(64)
_B = _NS >> 5
_PRED = np.stack([((_NS & 31) << 1) | 0, ((_NS & 31) << 1) | 1])          # (2, 64)
_SGN = 2.0 * np.stack([_SAIDA[_PRED[x], _B] for x in range(2)]) - 1.0    # (2, 64, 2): ±1


def viterbi(soft):
    """Viterbi de decisão suave. soft[i] > 0 significa bit codificado 1. Devolve os bits de informação."""
    n = len(soft) // 2
    y = np.asarray(soft[:2 * n], float).reshape(n, 2)
    met = np.full(64, 1e18)
    met[0] = 0.0
    escolha = np.zeros((n, 64), dtype=np.int8)
    for t in range(n):
        c = met[_PRED] - (_SGN[..., 0] * y[t, 0] + _SGN[..., 1] * y[t, 1])   # (2, 64)
        x = np.argmin(c, axis=0)
        escolha[t] = x
        met = c[x, _NS]
    s = int(np.argmin(met))
    out = np.zeros(n, dtype=np.int8)
    for t in range(n - 1, -1, -1):
        out[t] = s >> 5
        s = _PRED[escolha[t, s], s]
    return out


def _scrambler(estado, n):
    """Sequência do scrambler 802.11 (x^7 + x^4 + 1) a partir de um estado de 7 bits."""
    out = np.zeros(n, dtype=np.int8)
    for i in range(n):
        b = ((estado >> 6) ^ (estado >> 3)) & 1
        estado = ((estado << 1) | b) & 0x7F
        out[i] = b
    return out


POLARIDADE = 1 - 2 * _scrambler(0x7F, 127)      # p_0..p_126 dos pilotos
PIL_K = [-21, -7, 7, 21]
PIL_V = np.array([1, 1, 1, -1])
_SEQS = {e: _scrambler(e, 7) for e in range(1, 128)}


def descodifica_signal(sym, H):
    """sym: 64 amostras (sem prefixo). Devolve (ok, taxa, comprimento)."""
    X = np.fft.fft(sym)
    eq = np.array([X[k % 64] / H[k % 64] for k in DADOS])
    soft_rx = eq.real                                 # BPSK: +1 -> bit 1
    coded = soft_rx[INTERL]                           # desentrelaçar
    bits = viterbi(coded)
    if len(bits) != 24:
        return False, None, None
    paridade_ok = (bits[:18].sum() % 2) == 0
    cauda_ok = not bits[18:].any()
    taxa = TAXAS.get(tuple(int(b) for b in bits[:4]))
    comp = int(sum(int(b) << i for i, b in enumerate(bits[5:17])))
    return bool(paridade_ok and cauda_ok and taxa is not None), taxa, comp


def equaliza(sym, H, n):
    """FFT + equalização + correção de fase pelos pilotos (símbolo n; SIGNAL = 0). Devolve as 48 de dados."""
    X = np.fft.fft(sym)
    eq = np.array([X[k % 64] / H[k % 64] for k in DADOS])
    pil = np.array([X[k % 64] / H[k % 64] for k in PIL_K])
    theta = np.angle(np.sum(pil * PIL_V * POLARIDADE[n % 127]))
    return eq * np.exp(-1j * theta)


def descodifica_payload(simbolos, H, comp):
    """Só 6 Mbit/s (BPSK r=1/2), que é o caso de todas as tramas desta gravação. Devolve True se o CRC bater."""
    soft = np.concatenate([equaliza(s, H, n + 1).real[INTERL] for n, s in enumerate(simbolos)])
    bits = viterbi(soft)
    estado = next((e for e, q in _SEQS.items() if not np.any(bits[:7] ^ q)), None)
    if estado is None:
        return False
    desc = bits ^ _scrambler(estado, len(bits))
    d = desc[16:16 + 8 * comp]
    if len(d) < 8 * comp or comp < 4:
        return False
    dados = np.packbits(d.reshape(-1, 8)[:, ::-1], axis=1).ravel().tobytes()
    return zlib.crc32(dados[:-4]) == int.from_bytes(dados[-4:], "little")


# ---------------------------------------------------------------------------
def media_movel(x, n):
    c = np.concatenate(([0], np.cumsum(x)))
    out = np.empty(len(x), dtype=x.dtype)
    out[n - 1:] = c[n:] - c[:-n]
    out[:n - 1] = c[1:n]                              # arranque como no GNU Radio (soma parcial)
    return out


def sync_short(s, thr, min_plateau, W=48, min_gap=480, max_samples=540 * 80):
    """Devolve lista de (índice na entrada atrasada, estimativa grosseira ω1) - lógica do sync_short.cc."""
    s16 = np.concatenate((np.zeros(16, s.dtype), s[:-16]))            # Delay 16
    a = s * np.conj(s16)
    A = media_movel(a, W)
    P = media_movel(np.abs(s) ** 2, W + 16)
    cor = np.abs(A) / np.maximum(P, 1e-12)
    acima = cor > thr
    det = []
    estado, plateau, copied = "SEARCH", 0, 0
    idx = np.flatnonzero(acima)
    # percorre só as zonas acima do limiar + controla d_copied por diferença de índices
    i_ant = -1
    for i in idx:
        if i != i_ant + 1:
            plateau = 0
        if estado == "COPY":
            copied_i = i - inicio
            if copied_i >= max_samples:
                estado = "SEARCH"
        if plateau < min_plateau:
            plateau += 1
        else:
            if estado == "SEARCH" or (i - inicio) > min_gap:
                det.append((i, np.angle(A[i]) / 16))
                estado, inicio, plateau = "COPY", i, 0
        i_ant = i
    return det, s16


def sync_long(seg, sync_length):
    """seg: amostras desde o início da cópia (já com correção grosseira). Devolve (frame_start, d_freq_offset)."""
    n = sync_length
    if len(seg) < n + 63:
        return None
    # saída do FIR na posição p = sum_m conj(L[m]) seg[p+m]
    cor = np.correlate(seg[:n + 63], LTS_T, mode="valid")[:n]
    ordem = np.argsort(-np.abs(cor), kind="stable")[:4]
    start, off = n, 0.0
    for i in range(3):
        for k in range(i + 1, 4):
            pi, pk = ordem[i], ordem[k]
            first, second = (cor[pk], cor[pi]) if pi > pk else (cor[pi], cor[pk])
            d = abs(pi - pk)
            if d == 64:
                return min(pi, pk), np.angle(first * np.conj(second)) / 64
            if d in (63, 65):
                start, off = min(pi, pk), np.angle(first * np.conj(second)) / d
    return start, off


def corre(x, freq_offset=0.0, thr=0.56, min_plateau=2, sync_length=320, delay=None, noise=0.0, seed=0):
    if delay is None:
        delay = sync_length
    n = np.arange(len(x))
    s = x.astype(np.complex128) * np.exp(2j * np.pi * freq_offset / FS * n)
    if noise > 0:
        rng = np.random.default_rng(seed)
        s = s + noise / np.sqrt(2) * (rng.standard_normal(len(s)) + 1j * rng.standard_normal(len(s)))
    det, s16 = sync_short(s, thr, min_plateau)
    ok, crc_ok, wtf, estimativas = 0, 0, 0, []
    delta = delay - sync_length
    for j, (i0, w1) in enumerate(det):
        if j + 1 < len(det) and det[j + 1][0] - i0 < sync_length:
            wtf += 1                           # tag nova durante SYNC -> exceção "wtf" no bloco real
            continue
        k = np.arange(-400, sync_length + 3000)
        idx = i0 + k
        valido = (idx >= 0) & (idx < len(s16))
        seg_all = np.zeros(len(k), complex)
        seg_all[valido] = s16[idx[valido]] * np.exp(-1j * w1 * k[valido])
        seg = seg_all[400:]                    # seg[0] = primeira amostra copiada
        r = sync_long(seg, sync_length)
        if r is None:
            continue
        fs_, w2 = r
        estimativas.append((w1 - w2) * FS / (2 * np.pi))          # desvio total estimado (Hz)
        # COPY: amostra de saída com rel lê in_delayed -> seg[fs_ + rel - delta]
        def amostra(rel):
            p = fs_ + rel - delta + 400
            if p < 0 or p >= len(seg_all):     # fora do segmento (LENGTH errado numa trama mal sincronizada)
                return 0j
            return seg_all[p] * np.exp(1j * (fs_ + rel) * w2)
        lts1 = np.array([amostra(r) for r in range(0, 64)])
        lts2 = np.array([amostra(r) for r in range(64, 128)])
        sig = np.array([amostra(r) for r in range(144, 208)])
        H = (np.fft.fft(lts1) + np.fft.fft(lts2)) / 2 / np.where(LTS_F == 0, 1, LTS_F)
        boa, taxa, comp = descodifica_signal(sig, H)
        ok += boa
        if boa and taxa == 6:
            n_sym = -(-(16 + 8 * comp + 6) // 24)
            simbolos = [np.array([amostra(r) for r in range(128 + 80 * m + 16, 128 + 80 * m + 80)])
                        for m in range(1, n_sym + 1)]
            crc_ok += descodifica_payload(simbolos, H, comp)
    est = np.array(estimativas)
    return dict(detetadas=len(det), signal_ok=ok, crc_ok=crc_ok, wtf=wtf,
                est_mediana_khz=float(np.median(est) / 1e3) if len(est) else float("nan"),
                est_desvio_khz=float(np.std(est) / 1e3) if len(est) else float("nan"))


VARRIMENTOS = {
    "base": ("base", [0]),
    "freq": ("freq_offset", sorted({s * v for v in [0, 50e3, 100e3, 150e3, 200e3, 300e3, 400e3, 500e3,
                                                     550e3, 600e3, 620e3, 630e3, 650e3, 700e3, 800e3]
                                    for s in (1, -1)})),
    "delay": ("delay", [290, 300, 310, 314, 316, 318, 319, 320, 321, 322, 326, 330, 334, 336, 338, 340, 344, 350, 360]),
    "synclength": ("sync_length", [96, 128, 160, 192, 224, 256, 288, 320, 400, 480, 600, 1000, 1500, 2000, 2100]),
}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-file", required=True)
    p.add_argument("--varrimento", choices=VARRIMENTOS, default="base")
    p.add_argument("--noise", type=float, default=0.0)
    p.add_argument("--saida")
    a = p.parse_args()
    x = np.fromfile(a.input_file, np.complex64)
    nome, valores = VARRIMENTOS[a.varrimento]
    saida = Path(a.saida) if a.saida else Path(__file__).resolve().parent.parent / "data" / f"sim_{a.varrimento}.csv"
    saida.parent.mkdir(parents=True, exist_ok=True)
    with open(saida, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parametro", "valor", "detetadas", "signal_ok", "crc_ok", "wtf", "est_mediana_khz", "est_desvio_khz"])
        for v in valores:
            kw = {} if nome == "base" else {nome: v}
            r = corre(x, noise=a.noise, **kw)
            print(f"{nome}={v:>10g}  " + "  ".join(f"{k}={r[k]:.4g}" if isinstance(r[k], float) else f"{k}={r[k]}" for k in r), flush=True)
            w.writerow([nome, v, r["detetadas"], r["signal_ok"], r["crc_ok"], r["wtf"],
                        f"{r['est_mediana_khz']:.3f}", f"{r['est_desvio_khz']:.3f}"])
    print(f"Gravado em {saida}")


if __name__ == "__main__":
    main()
