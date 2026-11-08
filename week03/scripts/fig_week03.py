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
        refs=[(235, "≈ 235\n(2.º pico da LTS)"), (2014, "2014\n(menor intervalo\nentre tramas)")]),
}


SIM_PARAM = {"freq_offset": "freq-offset", "delay": "sync-delay", "sync_length": "sync-length"}


def le(caminho):
    """Lê o CSV do varrimento.py (GNU Radio) ou do sim_week03.py (simulação numpy)."""
    linhas = list(csv.DictReader(open(caminho)))
    param = linhas[0]["parametro"]
    if "crc_ok" in linhas[0]:                          # CSV da simulação
        param = SIM_PARAM[param]
        pts = [(float(l["valor"]), l["detetadas"], l["crc_ok"],
                "erro" if int(l["wtf"]) > int(l["detetadas"]) // 2 else "ok") for l in linhas]
        return param, pts, True, linhas
    pts = [(float(l["valor"]), l["detetadas"], l["descodificadas"], l["estado"]) for l in linhas]
    return param, pts, False, linhas


def fig_estimativa(linhas, saida):
    """Desvio total estimado (Sync Short + Sync Long) em função do desvio injetado."""
    x = [float(l["valor"]) / 1e3 for l in linhas]
    y = [float(l["est_mediana_khz"]) for l in linhas]
    fig, ax = plt.subplots(figsize=(6.5, 4.5), dpi=150)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(color="#e4e3df", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.plot([-800, 800], [-800, 800], color=REF, linewidth=1, linestyle="--", label="estimado = injetado")
    ax.plot(x, y, color=AZUL, linewidth=2, marker="o", markersize=5, label="Estimado (mediana das tramas)")
    for xr in (-625, 625):
        ax.axvline(xr, color=REF, linewidth=1, linestyle=":")
    ax.set_xlabel("Desvio injetado (kHz)")
    ax.set_ylabel("Desvio estimado (kHz)")
    ax.set_title("Simulação numpy: desvio estimado vs injetado", fontsize=9, loc="left", color="#52514e")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(saida, bbox_inches="tight")
    print(f"Gravado em {saida}")


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    caminho = Path(sys.argv[1])
    param, pts, sim, linhas = le(caminho)
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
                markersize=7, label="Erro (exceção \"wtf\")" if sim else "Erro / timeout")

    topo = max([p[1] for p in ok] + [1])
    for xr, txt in cfg["refs"]:
        ax.axvline(xr, color=REF, linewidth=1, linestyle="--")
        ax.text(xr, topo * 1.04, txt, ha="center", va="bottom", fontsize=7, color="#52514e",
                bbox=dict(facecolor="white", edgecolor="none", pad=1.5))

    if sim:
        ax.set_title("Simulação numpy (não é o GNU Radio)", fontsize=9, loc="left", color="#52514e")
    ax.set_xlabel(cfg["xlabel"])
    ax.set_ylabel("N.º de tramas")
    ax.set_ylim(0, topo * 1.2)
    ax.legend(frameon=False, loc="lower center", fontsize=8, ncol=3, bbox_to_anchor=(0.5, -0.32))
    fig.tight_layout()

    saida = Path(__file__).resolve().parent.parent / "figures" / (caminho.stem + ".png")
    saida.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(saida, bbox_inches="tight")
    print(f"Gravado em {saida}")
    if sim and param == "freq-offset":
        fig_estimativa(linhas, saida.with_name(caminho.stem + "_estimativa.png"))


if __name__ == "__main__":
    main()
