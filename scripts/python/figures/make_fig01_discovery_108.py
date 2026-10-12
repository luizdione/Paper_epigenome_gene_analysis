# -*- coding: utf-8 -*-
"""Fig. 1 — Constitutive methylation divergence (discovery, EXPANDED 108-gene panel, +/-15 kb).
A: positional distribution of 160 DMCs.  B: top genes by DMCs (module-coloured, incl. tick + milk).
C: PC1 (MOFA2 proxy) separating subspecies. PNG+PDF.
Tamanhos de tela/letra vêm do perfil de estilo em fig_paths.py (FIG_PROFILE=default|print)."""
import pandas as pd, numpy as np, csv, collections, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from fig_paths import DISC as RES, save_figure, P  # dados em epigenoma/results; saida em FIG_OUT
S=P("fig1")
ann=RES+"/fase5_fase7/annotation_table.csv"
mat=RES+"/methylation_matrix.csv"

rows=list(csv.DictReader(open(ann,encoding="utf-8")))
lab={'upstream_distal':'Upstream distal (<−2 kb TSS)','promoter':'Proximal promoter',
     'gene_body':'Gene body','downstream_distal':'Downstream distal'}
order=['upstream_distal','promoter','gene_body','downstream_distal']
feat=collections.Counter(r['feature'].strip('"') for r in rows)
md=[float(r['meth.diff']) for r in rows]
hyper=sum(1 for x in md if x>0); hypo=sum(1 for x in md if x<0); N=len(md)

# family colours (5 groups now: ISG, HSP, tick, milk, other)
ISG={'ISG15','MX1','MX2','OAS1Z','OAS1X','OAS1Y','RSAD2','OAS2'}
HSP={'HSP90AA1','HSF1','HSPA2','HSPA4','HSPA6','HSPA8','HSPB1','HSPB11','HSPA12B','HSPA1L','DNAJB4','HSP90AB1'}
TICK={'PADI4','WNT7A','IL22','TREM1','TREM2','CD83','MYO5A','TREML1','PRSS16','BOLA-DRB3','IL2','IL21','CXCL8','CCL2','HMGB1','EIF2AK2','GRHL3','DTX1'}
MILK={'CSN1S1','CSN1S2','CSN2','CSN3','PAEP','LALBA','NCAPG','ABCG2','DGAT1','GHR','SCD1','FASN'}
def cc(g):
    if g in ISG: return '#c0392b'
    if g in HSP: return '#e67e22'
    if g in TICK: return '#8e44ad'
    if g in MILK: return '#16a085'
    return '#2e6da4'
gc=collections.Counter(r['gene'].strip('"') for r in rows)
top=gc.most_common(15)

m=pd.read_csv(mat,index_col=0); X=m.dropna(axis=1)
Xc=((X-X.mean())/X.std(ddof=0).replace(0,1)).dropna(axis=1)
U,S_,Vt=np.linalg.svd(Xc.values-Xc.values.mean(0),full_matrices=False); pc1=U[:,0]*S_[0]
if pc1[[i for i,s in enumerate(X.index) if str(s).lower().startswith('ang')][0]]<0: pc1=-pc1
var=(S_**2/np.sum(S_**2))[0]*100
breed=['Nellore' if str(i).lower().startswith('nel') else 'Angus' for i in X.index]

fig=plt.figure(figsize=S["figsize"])
gs=fig.add_gridspec(1,3,width_ratios=S["ratios"],wspace=S["wspace"])
axA=fig.add_subplot(gs[0]); axB=fig.add_subplot(gs[1]); axC=fig.add_subplot(gs[2])

vals=[feat.get(k,0) for k in order]
axA.bar(range(len(order)),vals,color=['#4C72B0','#C44E52','#55A868','#8172B3'],edgecolor='black',lw=0.5)
for i,v in enumerate(vals): axA.text(i,v+1,'%d\n(%.0f%%)'%(v,100*v/N),ha='center',fontsize=S["annot"])
axA.set_xticks(range(len(order))); axA.set_xticklabels([lab[k] for k in order],fontsize=S["tick"],rotation=45,ha='right',rotation_mode='anchor')
axA.set_ylabel('DMCs',fontsize=S["axis"]); axA.set_ylim(0,max(vals)+12)
axA.set_title(S["title_a"]%(N,hyper,hypo),
              fontsize=S["title"],loc='left',fontweight='bold')
axA.spines[['top','right']].set_visible(False)

genes=[g for g,_ in top][::-1]; cnts=[c for _,c in top][::-1]; cols=[cc(g) for g in genes]
axB.barh(range(len(genes)),cnts,color=cols,edgecolor='black',lw=0.5)
axB.set_yticks(range(len(genes))); axB.set_yticklabels(genes,fontstyle='italic',fontsize=S["gene"])
for i,v in enumerate(cnts): axB.text(v+0.2,i,str(v),va='center',fontsize=S["annot_b"])
axB.set_xlabel('DMCs',fontsize=S["axis"]); axB.set_xlim(0,max(cnts)+2.5)
axB.set_title('B  Top genes by DMC count',fontsize=S["title"],loc='left',fontweight='bold')
handlesB=[Patch(facecolor='#c0392b',ec='black',label='Interferon (ISG)'),
          Patch(facecolor='#e67e22',ec='black',label='Heat-shock (HSP)'),
          Patch(facecolor='#8e44ad',ec='black',label='Tick resistance'),
          Patch(facecolor='#16a085',ec='black',label='Milk production'),
          Patch(facecolor='#2e6da4',ec='black',label='Other panel gene')]
if S["legend_below"]:
    # perfil print: a legenda maior cobriria as contagens das barras inferiores -> vai para baixo do painel
    axB.legend(handles=handlesB,fontsize=S["legend"],loc='upper left',bbox_to_anchor=(0.0,-0.20),ncol=2,
               frameon=True,handlelength=1.0,handletextpad=0.5,borderpad=0.35,labelspacing=0.3,columnspacing=1.0)
else:
    axB.legend(handles=handlesB,fontsize=S["legend"],loc='lower right',frameon=True)
axB.spines[['top','right']].set_visible(False)

for b,cl in [('Nellore','#1f77b4'),('Angus','#2ca02c')]:
    yy=[p for p,bb in zip(pc1,breed) if bb==b]
    axC.scatter(np.full(len(yy),0 if b=='Nellore' else 1),yy,s=S["scatter_s"],color=cl,edgecolor='black',lw=0.6,label=b,zorder=3)
axC.axhline(0,color='grey',lw=0.8,ls='--')
axC.set_xticks([0,1]); axC.set_xticklabels(['Nellore','Angus'],fontsize=S["c_tick"])
axC.set_ylabel('MOFA2-proxy factor 1\n(PC1, %.0f%% var.)'%var,fontsize=S["c_axis"])
axC.set_xlim(*S["c_xlim"]); axC.set_title('C  Subspecies separation',fontsize=S["title"],loc='left',fontweight='bold')
axC.spines[['top','right']].set_visible(False)

fig.tight_layout()
save_figure(fig,"Fig01_discovery_dmc",bbox='tight')   # Fig01_discovery_dmc.{png,pdf} + Fig01.png em FIG_OUT
print(f"Fig01 (108) OK. N={N} feat={dict(feat)} PC1={var:.1f}% top={top[:5]}")
