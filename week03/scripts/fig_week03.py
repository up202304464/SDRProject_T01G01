#!/usr/bin/env python3
"""
Week 03 - gráficos a partir dos CSV do varrimento.py.

Uso:
  python3 fig_week03.py ../data/varrimento_freq.csv
  python3 fig_week03.py ../data/varrimento_delay.csv
  python3 fig_week03.py ../data/varrimento_synclength.csv
A figura é gravada em ../figures/<nome do csv>.png
"""
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

AZUL, LARANJA = "#2a78d6", "#eb6834"      # séries: descodificadas, detetadas
REF = "#898781"                            # linhas de referência teóricas

CONFIG = {
    "freq-offset": dict(
        xlabel="Desvio de frequência injetado (kHz)", escala=1e-3,
        refs=[(-625, "−625 kHz\n(fs/32)"), (625, "625 kHz\n(fs/32)"),
              (-156.25, "−156 kHz\n(fs/128)"), (156.25, "156 kHz\n(fs/128)")]),
    "sync-delay": dict(
        xlabel="Atraso D da entrada in_delayed (amostras)", escala=1,
        refs=[(320, "D = 320\n(sync_length)"), (336, "D = 336\n(+ prefixo de 16)")]),
    "sync-length": dict(
        xlabel="sync_length (= atraso D) (amostras)", escala=1,
        refs=[(320, "320\n(por omissão)"), (480, "480\n(MIN_GAP)")]),
}


def le(caminho):
    linhas = list(csv.DictReader(open(caminho)))
    param = linhas[0]["parametro"]
    pts = [(float(l["valor"]), l["detetadas"], l["descodificadas"], l["estado"]) for l in linhas]
    return param, pts


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    caminho = Path(sys.argv[1])
    param, pts = le(caminho)
    cfg = CONFIG[param]

    ok = [(v * cfg["escala"], int(d), int(c)) for v, d, c, e in pts if e == "ok"]
    falhas = [v * cfg["escala"] for v, d, c, e in pts if e != "ok"]
    x = [p[0] for p in ok]

    fig, ax = plt.subplots(figsize=(8, 4.2), dpi=150)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(axis="y", color="#e4e3df", linewidth=0.8)
    ax.set_axisbelow(True)

    ax.plot(x, [p[1] for p in ok], color=LARANJA, linewidth=2, marker="o", markersize=5,
            label="Detetadas (Sync Short)")
    ax.plot(x, [p[2] for p in ok], color=AZUL, linewidth=2, marker="o", markersize=5,
            label="Descodificadas (CRC certo)")
    if falhas:
        ax.plot(falhas, [0] * len(falhas), linestyle="none", marker="x", color="#1f1f1d",
                markersize=7, label="Erro / timeout")

    topo = max([p[1] for p in ok] + [1])
    for xr, txt in cfg["refs"]:
        ax.axvline(xr, color=REF, linewidth=1, linestyle="--")
        ax.text(xr, topo * 1.04, txt, ha="center", va="bottom", fontsize=7, color="#52514e",
                bbox=dict(facecolor="white", edgecolor="none", pad=1.5))

    ax.set_xlabel(cfg["xlabel"])
    ax.set_ylabel("N.º de tramas")
    ax.set_ylim(0, topo * 1.2)
    ax.legend(frameon=False, loc="lower center", fontsize=8, ncol=3, bbox_to_anchor=(0.5, -0.32))
    fig.tight_layout()

    saida = Path(__file__).resolve().parent.parent / "figures" / (caminho.stem + ".png")
    saida.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(saida, bbox_inches="tight")
    print(f"Gravado em {saida}")


if __name__ == "__main__":
    main()
