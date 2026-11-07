# Week 03 — OFDM Sync Long e correção do desvio de frequência

Código analisado: `gr-ieee802-11/lib/sync_long.cc` (ramo `maint-3.10`). Os números de linha referem-se a esse ficheiro.
Parâmetros do sistema: fs = 20 MHz, sequência curta (STS) com período 16, sequência longa (LTS) com 2 × 64 amostras precedidas de um prefixo de 32.

---

## 1. O que faz o bloco OFDM Sync Long

O bloco tem **duas entradas** (linha 43): `in` (saída do Sync Short) e `in_delayed` (a mesma saída, atrasada por um bloco Delay no fluxograma). Funciona como uma máquina de três estados (linha 233):

| Estado | Linhas | O que faz |
|---|---|---|
| **SYNC** | 102–124 | Passa as primeiras `SYNC_LENGTH` (= 320) amostras da trama por um filtro adaptado à LTS (`d_fir`, linha 103) e guarda cada valor de correlação com a sua posição (linha 108). Ao chegar a 320 amostras chama `search_frame_start()` (linha 114). |
| **COPY** | 126–148 | Copia a trama a partir de `in_delayed`, corrige o desvio de frequência fino (linha 140) e **remove o prefixo cíclico** de cada símbolo (linha 139): mantém as 128 amostras da LTS e, de cada símbolo de 80 amostras seguinte, só as últimas 64. A saída são símbolos de 64 amostras prontos para a FFT. |
| **RESET** | 150–163 | Quando chega a tag de uma nova trama (linhas 89–90), escreve zeros até a saída estar alinhada a um múltiplo de 64 (linha 152) e volta a SYNC. Serve para o Stream to Vector seguinte nunca ficar desalinhado. |

**Filtro adaptado.** Os 64 coeficientes `LONG` (linhas 257–290) são a LTS **conjugada e invertida no tempo**: `LONG[k] = conj(L[63−k])` (por exemplo, `LONG[63] = 1.3868` ↔ `L[0]`, `LONG[62] = conj(L[1])`). A saída do filtro na posição *p* é a correlação Σ conj(L[m])·s[p+m], que tem um pico quando *p* é o início de um símbolo LTS.

**`search_frame_start()` (linhas 189–230).**
1. Ordena as 320 correlações por módulo (linha 194).
2. Entre os 4 maiores picos, procura um par separado de **64 amostras** (as duas cópias da LTS) (linhas 203–215).
3. Início da trama = o pico mais cedo (linha 216).
4. Desvio fino = `arg(first · conj(second)) / 64` (linha 217).
5. Tolera pares a 63 ou 65 amostras (linhas 221–227). Se não encontrar nenhum par, usa `d_frame_start = SYNC_LENGTH` (linha 201), o que dá uma trama mal cortada.

---

## 2. Onde é implementada a correção do desvio de frequência?

Um desvio Δf entre emissor e recetor roda cada amostra de ω = 2π·Δf/fs rad. Duas amostras que deviam ser iguais, separadas por um período *L*, diferem só pela rotação ω·L. O ângulo do produto, dividido por *L*, estima ω.

Nesta versão a correção é feita **em duas etapas, em dois blocos**:

| Etapa | Bloco | Sequência usada | Estimativa | Gama sem ambiguidade |
|---|---|---|---|---|
| Grosseira | Sync Short (linhas 81, 108, 119 do `sync_short.cc`) | STS, período 16 | ω̂₁ = arg(a[n]) / 16 | \|Δf\| < fs/32 = **625 kHz** |
| Fina | Sync Long (linha 217, aplicada na linha 140) | LTS, período 64 | ω̂₂ = arg(pico₁ · pico₂*) / 64 | \|Δf\| < fs/128 = **156,25 kHz** |

- **Compromisso:** um período maior dá mais precisão (o ruído pesa menos no ângulo), mas uma gama menor. A etapa curta remove o desvio grande e a longa estima o **resíduo** que sobra.
- **Sinal:** com o primeiro pico em *p* e o segundo em *p+64*, `first · conj(second)` tem ângulo −64ω. Por isso `d_freq_offset` = −ω, e a linha 140 multiplica por exp(+j·n·d_freq_offset) = exp(−jωn), anulando a rotação.
- **Passagem da estimativa entre blocos:** o Sync Short manda a sua estimativa no valor da tag `wifi_start`, que o Sync Long lê na linha 92 (`d_freq_offset_short`). O Sync Long volta a emitir a tag com o **desvio total** `d_freq_offset_short − d_freq_offset` (linha 135) para os blocos seguintes. *(A confirmar na semana 5/6: o Frame Equalizer deve usar este valor.)*
- Conversão para Hz: Δf = ω · fs / (2π).

---

## 3. Que parâmetros se podem variar?

| Parâmetro | Onde | Efeito esperado |
|---|---|---|
| Desvio de frequência injetado `freq_offset` | Fluxograma: *Signal Source* × *Multiply* à entrada | Experiência principal (ver secção 5). |
| `sync_length` | Parâmetro do bloco (linhas 41, 50) | Janela onde se procuram os picos da LTS. Se for curta demais, o 2.º pico (ou as 63 amostras que o filtro precisa a seguir) fica fora da janela e a trama perde-se. Se for longa demais, aumenta a latência e, **acima de 480**, uma trama nova pode chegar durante SYNC e disparar a exceção `"wtf"` (linhas 86–87). O MIN_GAP = 480 do Sync Short da semana 2 garante que isto não acontece com 320. *(hipótese a testar)* |
| Atraso da entrada `in_delayed` | Fluxograma: bloco Delay | Tem de ser igual a `sync_length` (ver secção 4). |
| Tolerância 63/65 e n.º de picos (4) | Código, linhas 203–227 | Só se recompilarmos. Mais tolerância ajuda com ruído, mas aceita mais pares errados. |
| Correção grosseira ligada/desligada | Código do Sync Short | Desligando-a, a gama total desce de ≈ 625 kHz para ≈ 156 kHz. Mostra para que serve cada etapa. |
| Ruído (σ do Noise Source) | Fluxograma | Repetir o varrimento de frequência com ruído. |

---

## 4. Porque existe uma entrada atrasada? O que acontece se mudarmos o atraso?

**Porquê.** Na fase SYNC o bloco consome as 320 primeiras amostras da trama (linhas 110–111, 169–170) **antes** de saber onde a trama começa. Um bloco de streaming não pode voltar atrás, por isso a solução é ter uma segunda cópia do sinal **atrasada exatamente 320 amostras**:

- No fim de SYNC, `d_offset` volta a 0 (linha 116).
- Nesse momento, `in_delayed[0]` é a amostra que tinha `d_offset = 0` em SYNC, ou seja, a primeira amostra da trama.
- Assim, `d_offset` em COPY indexa **as mesmas amostras** que em SYNC, e `rel = d_offset − d_frame_start` (linha 129) é 0 exatamente no início da LTS.

Conclusão: o atraso tem de ser igual a `SYNC_LENGTH` (320). O bloco "volta atrás no tempo" à custa de 320 amostras de latência e de memória.

**Se mudarmos o atraso D**, todos os cortes ficam deslocados de Δ = D − 320 amostras:

| Caso | Janela da FFT | Previsão |
|---|---|---|
| D = 320 | Correta | Descodifica. |
| 320 < D ≤ ≈ 336 (Δ até 16) | Começa Δ amostras **mais cedo**, dentro do prefixo cíclico do próprio símbolo | Continua a descodificar. O deslocamento só dá uma rampa de fase por subportadora, igual na LTS e nos dados, que a estimativa de canal absorve. Margem real < 16, por causa do espalhamento do canal. |
| D > ≈ 336 | Apanha a cauda do símbolo anterior | Interferência entre símbolos: deve falhar. |
| D < 320 (Δ < 0) | Começa \|Δ\| amostras **mais tarde** e entra no prefixo do símbolo seguinte | Interferência entre símbolos: deve degradar logo e falhar para poucas amostras. |

A assimetria (aguenta atrasar um pouco, não aguenta adiantar) vem de o prefixo cíclico estar **antes** de cada símbolo.

---

## 5. Experiências

> **As linhas "Previsão" são hipóteses teóricas, não medições.** A linha "Resultado" é preenchida com os valores medidos.

### Montagem e scripts

| Ficheiro | O que faz |
|---|---|
| `scripts/week_03_rx.py` | Recetor completo em Python (mesma montagem do `week_02_sweep.grc` + desvio de frequência + ruído opcional + Sync Long → FFT → Frame Equalizer → Decode MAC). Imprime `DETETADAS=… DESCODIFICADAS=…`. |
| `scripts/varrimento.py` | Corre o recetor para cada valor de um parâmetro e grava `data/varrimento_<nome>.csv`. |
| `scripts/fig_week03.py` | Gera `figures/varrimento_<nome>.png` a partir do CSV, com as linhas de referência teóricas. |

Duas métricas:
- **Detetadas** = tags `wifi_start` à saída do Sync Short (na semana 2 foram 772 com threshold 0,56).
- **Descodificadas** = tramas publicadas pelo Decode MAC, que só publica tramas com CRC certo.

**N₀** = descodificadas sem perturbação (`freq_offset = 0`, D = 320). Ainda não foi medido; pode ser menor que 772, se algumas tramas falharem o CRC mesmo sem perturbação. As previsões estão escritas em relação a N₀.

Como correr (no PC com GNU Radio + gr-ieee802-11, a partir de `week03/scripts/`):

```bash
F=../../week01/data/week01_signal.bin
python3 week_03_rx.py --input-file $F            # teste: deve dar DESCODIFICADAS = N₀ > 0
python3 varrimento.py freq       --input-file $F
python3 varrimento.py delay      --input-file $F
python3 varrimento.py synclength --input-file $F
python3 fig_week03.py ../data/varrimento_freq.csv
python3 fig_week03.py ../data/varrimento_delay.csv
python3 fig_week03.py ../data/varrimento_synclength.csv
```

Para repetir com ruído: `python3 varrimento.py freq --input-file $F --saida ../data/varrimento_freq_noise3.csv -- --noise 3`.
A frequência da portadora usada pelo Frame Equalizer é `--freq` (por omissão 5,89 GHz). *Confirmar a portadora da gravação da semana 1.*

### 5.1 Desvio de frequência — `varrimento.py freq`

| `freq_offset` (kHz) | 0 | ±50 | ±100 | ±200 | ±400 | ±600 | ±650 | ±700 |
|---|---|---|---|---|---|---|---|---|
| Previsão (descodificadas) | N₀ | N₀ | N₀ | N₀ | N₀ | ≈ N₀ | cai | ≈ 0 |
| Resultado | | | | | | | | |

Previsão: falha perto de **625 kHz**, quando a estimativa grosseira dá a volta (ambiguidade de 2π). As detetadas devem manter-se perto de 772 em toda a gama, porque c[n] usa |a[n]| e não depende do desvio.
Opcional: repetir com a correção grosseira desligada no código do Sync Short. A falha deve passar para ≈ **156 kHz**.

### 5.2 Atraso da entrada `in_delayed` — `varrimento.py delay`

| D | 300 | 310 | 316 | 320 | 324 | 330 | 336 | 340 | 350 |
|---|---|---|---|---|---|---|---|---|---|
| Previsão (descodificadas) | 0 | 0 | ↓ | N₀ | N₀ | ≈ N₀ | ↓ | 0 | 0 |
| Resultado | | | | | | | | | |

### 5.3 Janela de procura `sync_length` (com D = `sync_length`) — `varrimento.py synclength`

| `sync_length` | 128 | 192 | 256 | 320 | 400 | 480 | 520 | 600 |
|---|---|---|---|---|---|---|---|---|
| Previsão (descodificadas) | ? (pode perder o 2.º pico) | ? | N₀ | N₀ | N₀ | N₀ | erro (`"wtf"`)? | erro? |
| Resultado | | | | | | | | |

O limite inferior depende de onde o Sync Short começa a copiar a trama (ainda não sabemos). O superior testa a hipótese da secção 3 (MIN_GAP = 480). Se o bloco lançar a exceção, o `varrimento.py` regista o ponto como `erro` ou `timeout`.

---

## 6. Por confirmar

- Que o Sync Short desta versão emite a tag `wifi_start` com a sua estimativa ω̂₁ como valor (pressuposto da linha 92).
- Onde o desvio total (linha 135) é usado a seguir (Frame Equalizer, semanas 5/6).
- Para registar o desvio fino estimado em cada trama, acrescentar `d_freq_offset` ao `mylog` da linha 115 e recompilar (opcional).
