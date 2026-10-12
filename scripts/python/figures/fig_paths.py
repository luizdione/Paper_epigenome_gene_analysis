# -*- coding: utf-8 -*-
"""
fig_paths.py — caminhos e PERFIS DE ESTILO dos geradores de figura do projeto.

Caminhos: dados de resultado em `epigenoma/results/` (projeto-irmão, fora deste pacote); figuras
do pacote em `<pacote>/figuras/` (perfil "print", o padrão desde 2026-10-12) ou
`<pacote>/figuras_v1_pre_v2/` (perfil "v1"). Cada figura é gravada com o nome descritivo
(Fig0N_<desc>.png/.pdf) E com o nome usado pelo relatório (Fig0N.png), eliminando o passo manual
de renomear.

Perfis (variável de ambiente FIG_PROFILE; padrão = "print"):
  print   — PADRÃO. Figuras desenhadas no TAMANHO FINAL de impressão (174 mm = 6,85 pol. de
            largura, a largura de página dupla da Springer; altura <= 234 mm) com letras de
            9,5–11 pt, como a revista pede (8–12 pt no tamanho final). Motivo: o PDF de revisão
            (classe article, \textwidth = 345 pt = 4,79 pol.) e a diagramação da revista ENCOLHEM a
            imagem para a largura da coluna; uma tela de 14 pol. com letras de 8,5 pt vira ~2,9 pt
            impressa. Nenhum número, rótulo ou palavra muda em relação ao perfil v1: só tamanho de
            tela, corpo das letras, quebras de linha de título/eixo e a posição da legenda quando
            ela não cabe mais dentro do gráfico. Adotado pelo autor em 2026-10-12.
  v1      — valores EXATOS das figuras da submissão e da 1a entrega da revisão (2026-10-06), hoje
            arquivadas em `<pacote>/figuras/_archive_pre_v2_2026-10-12/`; regerar com este perfil
            reproduz aqueles PNG byte a byte (MD5 conferido 4x em 2026-10-11). Grava em
            `<pacote>/figuras_v1_pre_v2/` para não sobrescrever as figuras em uso.
            ("default", o nome antigo, é aceito como sinônimo de "v1".)

Overrides por variável de ambiente: EPIGENOMA_ROOT, EPIGENOMA_RESULTS, FIG_OUT, FIG_PROFILE, FIG_DPI.
"""
import os
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))                 # _fonte_reprodutibilidade
PKG = os.path.abspath(os.path.join(HERE, "..", "..", ".."))        # raiz do repositório
EPI = os.environ.get("EPIGENOMA_ROOT", PKG)
RESULTS = os.environ.get("EPIGENOMA_RESULTS", os.path.join(EPI, "results"))
SUB02 = os.path.join(RESULTS, "SUB02_runs_108")
DISC = os.path.join(RESULTS, "DISCOVERY_108_15kb")

PROFILE = os.environ.get("FIG_PROFILE", "print")
if PROFILE == "default":          # nome antigo do perfil v1
    PROFILE = "v1"
_DEFAULT_OUT = {"print": os.path.join(PKG, "figures"),
                "v1": os.path.join(PKG, "figures_v1")}
if PROFILE not in _DEFAULT_OUT:
    raise SystemExit("FIG_PROFILE deve ser 'print' ou 'v1' (recebido: %r)" % PROFILE)
FIG_OUT = os.environ.get("FIG_OUT", _DEFAULT_OUT[PROFILE])

PACKAGE_NAME = {
    "Fig01_discovery_dmc": "Fig01.png",
    "Fig02_celltype_dmc": "Fig02.png",
    "Fig03_bridge_concordance": "Fig03.png",
    "Fig04_module_enrichment": "Fig04.png",
    "Fig05_C6_validation": "Fig05.png",
    "Fig06_matched_atlas": "Fig06.png",
}

# Largura final de impressão (Springer, página dupla) usada pelo perfil "print".
PRINT_WIDTH_IN = 6.85

# ---------------------------------------------------------------- perfis de estilo
# Tamanhos em pontos tipográficos NA TELA DESENHADA. No perfil "default" a tela é muito maior
# que a página, por isso os pontos "parecem" grandes aqui mas saem pequenos no PDF.
STYLE = {
    "v1": {
        "dpi": 300,
        "fig1": dict(figsize=(14, 4.8), ratios=[1, 1.3, 0.8], wspace=0.42,
                     tick=8.5, gene=9.5, c_tick=10, c_xlim=(-0.6, 1.6), axis=11, c_axis=10, title=10.5,
                     legend=7.5, legend_below=False, annot=9, annot_b=8, scatter_s=85,
                     title_a='A  Positional distribution (n=%d DMCs;\n%d hyper / %d hypo, taurine)'),
        "fig2": dict(figsize=(9.2, 6.6), tick=8.6, axis=10, title=11, title_wrap=None, title_mode="axes",
                     legend=8.0, legend_outside=False, annot=8.0, xlabel_wrap=None, xlim_factor=1.16),
        "fig3": dict(figsize=(8.4, 6.6), axis=9.5, title=10.5, title_wrap=None, title_mode="axes",
                     legend=8.0, legend_loc="lower right", off=1.0, halo=2.2, label_fs=7.4, ms=9, ms_big=12),
        "fig4": dict(width=8.6, row=0.52, base=1.8, tick=9, axis=9.5, title=10.5, title_wrap=None,
                     title_mode="axes", legend=8.0, legend_loc="lower right", legend_anchor=None,
                     annot=8.4, xlim_factor=1.28),
        "fig5": dict(width=8.6, row=0.30, base=1.8, tick=7.6, axis=9.5, title=10.5, title_wrap=None,
                     title_mode="axes", legend=8.0, legend_outside=False),
        "fig6": dict(figsize=(8.2, 4.8), label=7.6, annot=8.6, axis=9.5, title=10.2, title_wrap=None,
                     title_mode="axes", xlabel_wrap=None, stagger=0.10, annot_wrap=None),
    },
    # "print": tela = tamanho final (6,85 pol. de largura); ticks/legendas 9,5–11 pt; títulos 9,5 pt
    # centrados na FIGURA (title_mode="fig") para não saírem da página quando os rótulos do eixo y
    # empurram o eixo para a direita; legendas fora do gráfico (Fig 1B, 2, 5) ou em quadrante vazio
    # (Fig 3 superior-direito, Fig 4 centro-direito) porque, maiores, deixam de caber onde estavam.
    "print": {
        "dpi": 600,
        "fig1": dict(figsize=(6.7, 4.2), ratios=[1, 1.5, 0.95], wspace=0.55,
                     tick=10, gene=10.5, c_tick=10, c_xlim=(-0.5, 1.5), axis=11, c_axis=10.5, title=9.5,
                     legend=9.5, legend_below=True, annot=9.5, annot_b=9, scatter_s=70,
                     title_a='A  Positional distribution\n(n=%d DMCs; %d hyper /\n%d hypo, taurine)'),
        "fig2": dict(figsize=(PRINT_WIDTH_IN, 5.9), tick=10.5, axis=11, title=9.5, title_wrap=76,
                     title_mode="fig", legend=10, legend_outside=True, annot=9.5, xlabel_wrap=40,
                     xlim_factor=1.24),
        "fig3": dict(figsize=(PRINT_WIDTH_IN, 6.2), axis=11, title=9.5, title_wrap=76, title_mode="fig",
                     legend=10, legend_loc="upper right", off=1.6, halo=2.6, label_fs=10, ms=9, ms_big=12),
        "fig4": dict(width=PRINT_WIDTH_IN, row=0.40, base=1.9, tick=11, axis=11, title=9.5, title_wrap=76,
                     title_mode="fig", legend=10, legend_loc="center right", legend_anchor=(1.0, 0.40),
                     annot=10, xlim_factor=1.65),
        "fig5": dict(width=PRINT_WIDTH_IN, row=0.155, base=1.7, tick=9.5, axis=11, title=9.5,
                     title_wrap=76, title_mode="fig", legend=10, legend_outside=True),
        "fig6": dict(figsize=(PRINT_WIDTH_IN, 4.4), label=10, annot=10.5, axis=11, title=9.5,
                     title_wrap=76, title_mode="fig", xlabel_wrap=48, stagger=0.17, annot_wrap=14),
    },
}
DPI = int(os.environ.get("FIG_DPI", STYLE[PROFILE]["dpi"]))


def P(fig_key):
    """Dicionário de estilo da figura (`fig1`..`fig6`) no perfil ativo."""
    return STYLE[PROFILE][fig_key]


def wrap(text, width):
    """Requebra cada linha de `text` que exceda `width` colunas (None = inalterado), equilibrando
    as duas partes para não deixar palavra órfã. Só muda quebras de linha, nunca o texto."""
    if not width:
        return text
    out = []
    for line in text.split("\n"):
        if len(line) <= width:
            out.append(line)
            continue
        parts = textwrap.wrap(line, width, break_on_hyphens=False)
        if len(parts) == 2 and len(parts[-1]) < 0.4 * width:
            parts = textwrap.wrap(line, len(line) // 2 + 4, break_on_hyphens=False)
        out.extend(parts or [""])
    return "\n".join(out)


def set_title(fig, ax, text, S):
    """Título no próprio eixo (perfil default) ou centrado na figura (perfil print).
    Devolve a fração superior da figura a reservar no tight_layout (None quando fica no eixo).
    O tight_layout do matplotlib 3.10 já desconta o suptitle; a reserva extra é só a folga
    entre o título e o gráfico (TITLE_GAP_IN, em polegadas)."""
    t = wrap(text, S.get("title_wrap"))
    if S.get("title_mode", "axes") == "axes":
        ax.set_title(t, fontsize=S["title"], fontweight="bold")
        return None
    fig.suptitle(t, fontsize=S["title"], fontweight="bold", y=0.995, va="top")
    h_in = fig.get_size_inches()[1]
    return 1.0 - TITLE_GAP_IN / h_in


TITLE_GAP_IN = float(os.environ.get("FIG_TITLE_GAP_IN", "0.08"))


def layout(fig, top=None, bottom=0.0):
    """tight_layout com reserva opcional para o suptitle (top) e para legenda externa (bottom)."""
    if top is None and not bottom:
        fig.tight_layout()
    else:
        fig.tight_layout(rect=(0, bottom, 1, 1.0 if top is None else top))


def save_figure(fig, stem, dpi=None, bbox=None):
    """Grava PNG+PDF descritivos e, se houver, a cópia com o nome do relatório."""
    dpi = dpi or DPI
    os.makedirs(FIG_OUT, exist_ok=True)
    kw = {"facecolor": "white"}
    if bbox:
        kw["bbox_inches"] = bbox
    fig.savefig(os.path.join(FIG_OUT, stem + ".png"), dpi=dpi, **kw)
    fig.savefig(os.path.join(FIG_OUT, stem + ".pdf"), **kw)
    if stem in PACKAGE_NAME:
        fig.savefig(os.path.join(FIG_OUT, PACKAGE_NAME[stem]), dpi=dpi, **kw)
        print("OK [%s, %d dpi]: %s -> %s" % (PROFILE, dpi, stem, PACKAGE_NAME[stem]))
    else:
        print("OK [%s, %d dpi]: %s" % (PROFILE, dpi, stem))
