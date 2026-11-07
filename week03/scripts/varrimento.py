#!/usr/bin/env python3
"""
Week 03 - varrimentos. Corre o week_03_rx.py uma vez por valor e grava um CSV.

Uso:
  python3 varrimento.py freq       --input-file ../../week01/data/week01_signal.bin
  python3 varrimento.py delay      --input-file ...
  python3 varrimento.py synclength --input-file ...

Argumentos a mais são passados ao recetor, por exemplo, para repetir com ruído:
  python3 varrimento.py freq --input-file ... --saida ../data/varrimento_freq_noise3.csv -- --noise 3

Cada corrida tem um limite de tempo; se o recetor bloquear ou der erro
(por exemplo, a exceção "wtf" do Sync Long), o ponto fica com estado "erro"/"timeout".
"""
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
RX = AQUI / "week_03_rx.py"

_f = [0, 50e3, 100e3, 150e3, 200e3, 300e3, 400e3, 500e3, 550e3, 600e3, 625e3, 650e3, 700e3, 800e3]
PRESETS = {
    # nome: (opção do recetor, valores, ficheiro de saída por omissão)
    "freq": ("--freq-offset", sorted({s * v for v in _f for s in (1, -1)}), "varrimento_freq.csv"),
    "delay": ("--sync-delay", [290, 300, 310, 314, 318, 320, 322, 326, 330, 334, 336, 338, 340, 350, 360],
              "varrimento_delay.csv"),
    # sync-length muda também o atraso (o recetor põe sync-delay = sync-length)
    "synclength": ("--sync-length", [128, 160, 192, 224, 256, 288, 320, 400, 480, 520, 600],
                   "varrimento_synclength.csv"),
}

PADRAO = re.compile(r"DETETADAS=(\d+)\s+DESCODIFICADAS=(\d+)")


def corre(opcao, valor, input_file, extra, timeout):
    cmd = [sys.executable, str(RX), "--input-file", input_file, opcao, str(valor)] + extra
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return "", "", "timeout"
    m = PADRAO.search(r.stdout)
    if r.returncode != 0 or not m:
        erro = (r.stderr.strip().splitlines() or ["?"])[-1]
        print(f"    erro: {erro}", file=sys.stderr)
        return "", "", "erro"
    return int(m.group(1)), int(m.group(2)), "ok"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("varrimento", choices=PRESETS)
    p.add_argument("--input-file", required=True)
    p.add_argument("--saida", help="CSV de saída (por omissão ../data/<varrimento>.csv)")
    p.add_argument("--timeout", type=float, default=300, help="segundos por corrida")
    # tudo o que não for reconhecido aqui (por exemplo --noise 3) segue para o week_03_rx.py
    a, extra = p.parse_known_args()
    extra = [x for x in extra if x != "--"]

    opcao, valores, nome = PRESETS[a.varrimento]
    saida = Path(a.saida) if a.saida else AQUI.parent / "data" / nome
    saida.parent.mkdir(parents=True, exist_ok=True)

    with open(saida, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parametro", "valor", "detetadas", "descodificadas", "estado"])
        for v in valores:
            det, dec, estado = corre(opcao, v, a.input_file, extra, a.timeout)
            print(f"{opcao} {v:>10g}  detetadas={det!s:>5}  descodificadas={dec!s:>5}  {estado}")
            w.writerow([opcao.lstrip("-"), v, det, dec, estado])
            f.flush()

    print(f"\nGravado em {saida}")


if __name__ == "__main__":
    main()
