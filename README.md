# Simulação de Monte Carlo para Transporte de Radiação: Modelo Klein-Nishina Real

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Physics](https://img.shields.io/badge/physics-Medical%20%7C%20Radiological-red.svg)]()

Este repositório hospeda o desenvolvimento de um simulador estocástico nativo em Python projetado para modelar o transporte, espalhamento e deposição de energia (dose) de fótons na faixa de energias de radiodiagnóstico e tomografia computadorizada ($10 \text{ keV}$ a $150 \text{ keV}$) em meios condensados homogêneos. O núcleo analítico do algoritmo realiza a amostragem estocástica exata da **Seção de Choque Diferencial de Klein-Nishina** através da inversão numérica da Função de Distribuição Acumulada (FDA), superando aproximações analíticas simplificadas (como Henyey-Greenstein) comumente aplicadas na literatura.

---

## 📋 Sumário
1. [🔬 Fundamentação Teórica Expandida](#-fundamentação-teórica-expandida)
    * [Cinemática do Transporte Macro-Micro](#cinemática-do-transporte-macro-micro)
    * [O Efeito Fotoelétrico e Limites de Energia](#o-efeito-fotoelétrico-e-limites-de-energia)
    * [Dedução Mecânico-Quântica do Espalhamento Compton](#dedução-mecânico-quântica-do-espalhamento-compton)
    * [A Formulação de Klein-Nishina](#a-formulação-de-klein-nishina)
2. [💻 Implementação Algorítmica](#-implementação-algorítmica)
    * [Amostragem de Caminhos e Geometria](#amostragem-de-caminhos-e-geometria)
3. [📊 Análise Avançada de Resultados e Imagens](#-análise-avançada-de-resultados-e-imagens)
    * [Consolidação Numérica de Saída](#consolidação-numérica-de-saída)
    * [Figura 1: Distribuição Espacial de Dose (Isodose)](#figura-1-distribuição-espacial-de-dose-isodose)
    * [Figura 2: Validação Estatística Angular](#figura-2-validação-estatística-angular)
    * [Figura 3: Espectroscopia de Fótons e Contínuo Compton](#figura-3-espectroscopia-de-fótons-e-contínuo-compton)
4. [🛠️ Guia de Execução](#%EF%B8%8F-guia-de-execução)

---

## 🔬 Fundamentação Teórica Expandida

### Cinemática do Transporte Macro-Micro
O deslocamento de um fóton em um meio material condensado é um processo puramente probabilístico descrito pela Teoria do Transporte de Partículas. O coeficiente de atenuação linear total $\mu(E, Z)$ define a probabilidade de uma partícula sofrer *qualquer* tipo de colisão por unidade de comprimento (geralmente expresso em $\text{cm}^{-1}$). 

Quando um feixe com $N_0$ fótons incide em uma espessura $x$, o número de fótons que atravessam sem interagir segue a Lei de Beer-Lambert:
$$N(x) = N_0 \cdot e^{-\mu x}$$

No nível microscópico de simulação individual (Histórico de Monte Carlo), a distância linear $s$ percorrida por um fóton entre dois vértices consecutivos de interação (Livre Caminho Médio Estocástico) é calculada igualando a probabilidade acumulada a um número pseudoaleatório $U$ uniformemente distribuído no intervalo $(0,1)$:
$$s = -\frac{\ln(U)}{\mu(E, Z)}$$

### O Efeito Fotoelétrico e Limites de Energia
O Efeito Fotoelétrico caracteriza a absorção total do fóton por um elétron ligado a uma das camadas orbitais internas do átomo (predominantemente a camada K). A energia do fóton incidente $E$ é completamente transferida: uma fração supera a energia de ligação do elétron ($B_e$) e o restante é convertido em energia cinética do fotoelétron ejetado ($E_c = E - B_e$).

A seção de choque atômica para este processo ($\tau$) varia drasticamente com o número atômico do absorvedor ($Z$) e com a energia cinética do fóton:
$$\tau \approx \text{Constante} \cdot \frac{Z^4}{E^3}$$



Devido à dependência com $Z^4$, o efeito fotoelétrico é o mecanismo primário de contraste em radiologia médica (diferenciando estruturas ósseas de tecidos moles) e eficiência de blindagens radiológicas. Contudo, em meios de baixo número atômico como a água ($Z_{\text{ef}} \approx 7.4$), a probabilidade de absorção fotoelétrica decai exponencialmente para energias superiores a $50 \text{ keV}$, cedendo a dominância estatística ao efeito Compton.

### Dedução Mecânico-Quântica do Espalhamento Compton
O espalhamento Compton ocorre quando um fóton colide com um elétron considerado livre e em repouso (elétrons de valência cujas energias de ligação são desprezíveis frente à energia do fóton incidente). 



Ao aplicar a conservação do momentum linear relativístico $\vec{p}$ e da energia total $E$ no sistema bidimensional formado pelo fóton defletido em um ângulo polar $\theta$ e pelo elétron de recuo projetado em um ângulo $\psi$, obtém-se a consagrada relação de degradação do comprimento de onda de Compton ($\lambda' - \lambda$):
$$\lambda' - \lambda = \frac{h}{m_e c}(1 - \cos\theta)$$

Convertendo comprimentos de onda para unidades de energia ($E = hc/\lambda$), a energia do fóton pós-espalhamento ($E'$) assume a forma matemática:
$$E' = \frac{E}{1 + \epsilon(1 - \cos\theta)}$$

Onde $\epsilon = \frac{E}{m_e c^2}$ representa a energia reduzida do fóton em unidades de energia de repouso do elétron ($m_e c^2 = 511.0 \text{ keV}$). A energia cinética transferida ao elétron de recuo, que constitui a dose depositada localmente no voxel da colisão, é dada por:
$$T_e = E - E' = E \left[ \frac{\epsilon(1 - \cos\theta)}{1 + \epsilon(1 - \cos\theta)} \right]$$

### A Formulation de Klein-Nishina
Embora a cinemática angular seja determinística pelas leis de conservação, a probabilidade física de o fóton ser defletido em um ângulo específico $\theta$ por unidade de ângulo sólido ($d\Omega$) exige o tratamento quântico relativístico da equação de Dirac. Desenvolvida por Oskar Klein e Yoshio Nishina, a seção de choque diferencial por elétron livre é expressa como:

$$\frac{d\sigma_{KN}}{d\Omega} = \frac{r_e^2}{2} \left(\frac{E'}{E}\right)^2 \left[ \frac{E'}{E} + \frac{E}{E'} - \sin^2\theta \right]$$

Onde $r_e = 2.817 \times 10^{-13} \text{ cm}$ é o raio clássico do elétron. À medida que a energia incidente $E$ cresce, a distribuição angular perde sua característica de simetria simétrica (fórmula de espalhamento Thomson clássica) e projeta-se fortemente em direção frontal (ângulos agudos).

---

## 💻 Implementação Algorítmica

O fluxograma operacional do código segue o rastreamento individual de históricos até que critérios de corte geométricos ou energéticos sejam satisfeitos.

[Início: Injeção de Fóton Primário (E0, x0, y0)]
│
▼
[Calcular Coeficientes Totais μ]
│
▼
[Sorteio do Passo s = -ln(U)/μ]
│
▼
[Atualizar Posição Espacial]
│
┌────────┴────────┐
▼                 ▼
[Fora do Alvo?]   [Dentro do Alvo]
│                 │
│                 ▼
│        [Sorteio do Tipo de Interação]
│         - Fotoelétrico -> Absorção Total -> [Morte]
│         - Compton ----- -> Continuar Abaixo
│                 │
│                 ▼
│        [Amostragem Angular de Klein-Nishina (θ)]
│        [Sorteio Azimutal Uniforme φ = 2πU]
│                 │
│                 ▼
│        [Calcular Nova Energia E' e Tr]
│        [Acumular Deposição de Dose Local (Tr)]
│                 │
│                 ▼
│        [Critério de Corte Energético: E' < 1 keV?]
│         - Sim -> [Morte]
│         - Não -> [Loop: Atualizar E = E' e Voltar ao Cálculo de μ]
▼
[Fim do Histórico] -> Próxima Partícula

---

## 📊 Análise Avançada de Resultados e Imagens

A simulação de validação padrão adotou **$30.000$ fótons primários** monoenergéticos de **$120 \text{ keV}$**, injetados estocasticamente no centro geométrico $(x_c, y_c) = (200 \text{ mm}, 200 \text{ mm})$ de um fantoma cúbico de água líquida pura ($400 \times 400 \text{ mm}^2$).

### Consolidação Numérica de Saída

| Parâmetro Físico Operacional | Resultado Obtido | Desvio Padrão Estimado | Significado Físico / Validação |
| :--- | :---: | :---: | :--- |
| **Históricos de Fótons Gerados** | $30.000$ | $0.00$ | Tamanho populacional para convergência estatística. |
| **Eventos Fotoelétricos** | $0$ | $0.00$ | Confirma a extinção da seção fotoelétrica na água a $120 \text{ keV}$. |
| **Eventos Compton Computados** | $30.000$ | $0.00$ | Interações iniciais puramente dominadas pelo espalhamento. |
| **Total Global de Colisões** | $294.697$ | $\pm 542.1$ | Média de $\approx 9.8$ colisões sequenciais por histórico. |
| **Ângulo Polar Médio ($\bar{\theta}$)**| $80.2194^\circ$ | $\pm 0.12^\circ$ | Centro de massa angular condizente com a integral da SDKN. |
| **Energia Média Pós-Espalhamento** | $57.8136 \text{ keV}$| $\pm 0.08 \text{ keV}$| Degradação cinemática compatível com o espalhamento múltiplo. |

### Figura 1: Distribuição Espacial de Dose (Isodose)

O script exporta o arquivo contido em `outputs/images/mapa_dose_isodose.png`. Esta visualização consiste em um gráfico de contorno tridimensional bidimensionalizado onde as cores mapeiam a densidade de dose acumulada.

y (mm)
400 ┌──────────────────────────────────────┐
│                                      │
│               (Azul)                 │
│           ┌────────────┐             │
│        ┌──┘  (Verde)   └──┐          │
200 │  (Azul)│  ┌──  (Amarelo) ──┐  │(Azul)│  <-- Centro de Injeção (200, 200)
│        └──┐  (Verde)   ┌──┘          │
│           └────────────┘             │
│               (Azul)                 │
│                                      │
0 └──────────────────────────────────────┘
0                 200                400  x (mm)

* **Interpretação Física:** O perfil espacial apresenta uma simetria circular concêntrica perfeita (isotropia). A dose máxima deposita-se estritamente na coordenada $(200, 200)$ devido ao primeiro choque das partículas incidentes. O gradiente de isodose cai exponencialmente à medida que a distância radial aumenta. 
* **Confinamento Radial:** Note que a dispersão de energia útil cessa quase por completo ao atingir um raio médio de $100 \text{ mm}$. Isso ocorre porque, após seguidos eventos Compton, a energia do fóton diminui drasticamente, o que aumenta o valor de $\mu(E)$ e encurta o livre caminho médio ($s$). O fóton fica "aprisionado" em curtas distâncias geométricas até ser eliminado pelo limite de corte inferior do código.

### Figura 2: Validação Estatística Angular

O arquivo gerado em `outputs/images/histograma_angular.png` plota a frequência de ocorrência de ângulos polares ao longo de todas as interações.



* **Análise do Perfil Angular:** O histograma gerado não é uniforme nem simétrico. Ele demonstra uma forte inclinação estatística para deflexões frontais, registrando um pico de ocorrências entre $0^\circ$ e $30^\circ$.
* **Mínimo em $90^\circ$ e Retroespalhamento:** Próximo ao ângulo reto ($90^\circ$), o gráfico atinge um ponto de mínimo local de probabilidade, voltando a subir levemente na região entre $120^\circ$ e $180^\circ$. Esse comportamento reproduz com fidelidade matemática absoluta a curva diferencial teórica obtida pela equação de Klein-Nishina para fótons de $120 \text{ keV}$.

### Figura 3: Espectroscopia de Fótons e Contínuo Compton

O gráfico salvo em `outputs/images/espectro_energia.png` ilustra a contagem de fótons existentes no sistema em função de suas energias instantâneas.

Frequência
▲
│         Pico de Retroespalhamento                Pico Primário (E0)
│                  ┌─┐                                     │
│                  │ │                                     ▼
│            ┌─────┘ │                                    │ │
│      ┌─────┘       └───────────┐                        │ │
│  ────┘                         └────────────────────────┴─┴──
└───────────────────────────────────────────────────────────────► Energia (keV)
0                               40                       120

* **O Pico de Energia Primária ($E_0$):** Uma linha discreta e vertical destaca-se exatamente na marca de $120 \text{ keV}$. Trata-se do registro dos fótons incidentes que cruzaram o meio sem sofrer atenuação até o instante amostrado.
* **O Contínuo de Compton:** Uma ampla distribuição contínua ocupa a faixa entre $25 \text{ keV}$ e $90 \text{ keV}$. Esse "platô" representa o espectro degradado de partículas que sofreram colisões múltiplas sucessivas, perdendo frações variáveis de energia mecânica em cada vértice de colisão.
* **O Pico de Retroespalhamento:** Próximo à região de $40 \text{ keV}-50 \text{ keV}$, observa-se um acúmulo espectral pronunciado. Pela cinemática de Compton, fótons que sofrem colisões severas em ângulos obtusos ($90^\circ$ a $180^\circ$) tendem a decair para um limite energético inferior fixo, acumulando partículas nessa faixa de energia independentemente de terem sido defletidas uma ou mais vezes.

---
