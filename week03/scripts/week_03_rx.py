#!/usr/bin/env python3
"""
Week 03 - recetor 802.11a (Sync Short + Sync Long + descodificação) para varrimentos.

Montagem igual à do week_02_sweep.grc, com:
  * desvio de frequência injetado (Signal Source x Multiply)    --freq-offset
  * ruído opcional (Noise Source + Add)                          --noise
  * cadeia Sync Long -> FFT -> Frame Equalizer -> Decode MAC
  * atraso da entrada in_delayed do Sync Long                    --sync-delay
  * janela de procura do Sync Long                               --sync-length

No fim imprime UMA linha, que o varrimento.py lê:
  DETETADAS=<tramas detetadas pelo Sync Short> DESCODIFICADAS=<tramas com CRC certo>

Exemplo:
  python3 week_03_rx.py --input-file ../../week01/data/week01_signal.bin --freq-offset 200e3
"""
import argparse
import sys
import time

import numpy as np
import pmt
from gnuradio import analog, blocks, fft, gr
from gnuradio.fft import window

import ieee802_11


class ContaTags(gr.sync_block):
    """Passa o sinal sem alterações e conta as tags 'wifi_start' (= tramas detetadas pelo Sync Short)."""

    def __init__(self, chave="wifi_start"):
        gr.sync_block.__init__(self, name="conta_tags",
                               in_sig=[np.complex64], out_sig=[np.complex64])
        self.chave = pmt.intern(chave)
        self.n = 0

    def work(self, input_items, output_items):
        n = len(input_items[0])
        output_items[0][:] = input_items[0]
        for t in self.get_tags_in_window(0, 0, n):
            if pmt.eq(t.key, self.chave):
                self.n += 1
        return n


class ContaMensagens(gr.basic_block):
    """Conta as mensagens publicadas pelo Decode MAC (só publica tramas com CRC certo)."""

    def __init__(self):
        gr.basic_block.__init__(self, name="conta_mensagens", in_sig=None, out_sig=None)
        self.n = 0
        self.message_port_register_in(pmt.intern("in"))
        self.set_msg_handler(pmt.intern("in"), self._recebe)

    def _recebe(self, msg):
        self.n += 1


def equalizador_ls():
    """O nome do enum muda entre versões dos bindings."""
    if hasattr(ieee802_11, "LS"):
        return ieee802_11.LS
    return ieee802_11.Equalizer.LS


class Week03(gr.top_block):

    def __init__(self, a):
        gr.top_block.__init__(self, "Week 03")
        fs = a.samp_rate
        W = a.window_size

        # --- fonte + desvio de frequência (+ ruído) -------------------------
        src = blocks.file_source(gr.sizeof_gr_complex, a.input_file, False)
        lo = analog.sig_source_c(fs, analog.GR_COS_WAVE, a.freq_offset, 1, 0)
        desvio = blocks.multiply_cc(1)
        self.connect(src, (desvio, 0))
        self.connect(lo, (desvio, 1))
        sinal = desvio

        if a.noise > 0:
            ruido = analog.noise_source_c(analog.GR_GAUSSIAN, a.noise, a.seed)
            soma = blocks.add_cc(1)
            self.connect(sinal, (soma, 0))
            self.connect(ruido, (soma, 1))
            sinal = soma

        # --- deteção (igual à semana 2) -------------------------------------
        atraso16 = blocks.delay(gr.sizeof_gr_complex, 16)
        conj = blocks.conjugate_cc()
        mult = blocks.multiply_cc(1)
        media_c = blocks.moving_average_cc(W, 1, 4000, 1)
        modulo = blocks.complex_to_mag(1)
        pot = blocks.complex_to_mag_squared(1)
        media_p = blocks.moving_average_ff(W + 16, 1, 4000, 1)
        divide = blocks.divide_ff(1)
        sync_short = ieee802_11.sync_short(a.threshold, a.min_plateau, False, False)

        self.connect(sinal, pot, media_p, (divide, 1))
        self.connect(sinal, atraso16, conj, (mult, 1))
        self.connect(sinal, (mult, 0))
        self.connect(mult, media_c, modulo, (divide, 0))
        self.connect(atraso16, (sync_short, 0))
        self.connect(media_c, (sync_short, 1))
        self.connect(divide, (sync_short, 2))

        # --- Sync Long e descodificação -------------------------------------
        self.detecoes = ContaTags()
        atraso_long = blocks.delay(gr.sizeof_gr_complex, a.sync_delay)
        sync_long = ieee802_11.sync_long(a.sync_length, False, False)
        s2v = blocks.stream_to_vector(gr.sizeof_gr_complex, 64)
        transf = fft.fft_vcc(64, True, window.rectangular(64), True, 1)
        equal = ieee802_11.frame_equalizer(equalizador_ls(), a.freq, fs, False, False)
        decode = ieee802_11.decode_mac(False, False)
        self.tramas = ContaMensagens()

        self.connect(sync_short, self.detecoes)
        self.connect(self.detecoes, (sync_long, 0))
        self.connect(self.detecoes, atraso_long, (sync_long, 1))
        self.connect(sync_long, s2v, transf, equal, decode)
        self.msg_connect((decode, "out"), (self.tramas, "in"))


def argumentos():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-file", required=True, help="gravação .bin (complex64)")
    p.add_argument("--samp-rate", type=float, default=20e6)
    p.add_argument("--freq", type=float, default=5.89e9,
                   help="frequência da portadora da gravação (Hz); usada pelo Frame Equalizer")
    p.add_argument("--window-size", type=int, default=48)
    p.add_argument("--threshold", type=float, default=0.56)
    p.add_argument("--min-plateau", type=int, default=2)
    p.add_argument("--freq-offset", type=float, default=0.0, help="desvio injetado (Hz)")
    p.add_argument("--noise", type=float, default=0.0, help="amplitude do Noise Source (0 = sem ruído)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--sync-length", type=int, default=320)
    p.add_argument("--sync-delay", type=int, default=None,
                   help="atraso da entrada in_delayed (por omissão = sync-length)")
    a = p.parse_args()
    if a.sync_delay is None:
        a.sync_delay = a.sync_length
    return a


def main():
    a = argumentos()
    tb = Week03(a)
    tb.start()
    tb.wait()
    time.sleep(0.5)  # deixa chegar as últimas mensagens do Decode MAC
    print(f"DETETADAS={tb.detecoes.n} DESCODIFICADAS={tb.tramas.n}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
