
import pandas as pd, numpy as np, math
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
xlsx=r'postgwas_ad_pdlbd/results/22_supplement_tables_lipid8_F2_AD/metabolic_factor_triplet_supplementary_tables_v12_submission_clean.xlsx'
lipid_res=r'D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\manuscript_formal_genomicsem_fit\lipid_final8_formal_genomicsem_results.tsv'
non_res=r'D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion\manuscript_formal_genomicsem_fit\nonlipid_final8_formal_genomicsem_results.tsv'
out_png=r'figures/figure2_factor_models_disease_screen_v13_final_polished.png'; out_pdf=r'figures/figure2_factor_models_disease_screen_v13_final_polished.pdf'; out_svg=r'figures/figure2_factor_models_disease_screen_v13_final_polished.svg'
axis_order=['lipid8_F1','lipid8_F2','lipid8_F3','nonlipid8_F1','nonlipid8_F2','nonlipid8_F3']
short_labels={'lipid8_F1':'TG/VLDL','lipid8_F2':'HDL-core','lipid8_F3':'CE-lipid','nonlipid8_F1':'Ketone','nonlipid8_F2':'Amino acid','nonlipid8_F3':'Energy'}
colors={'lipid8_F1':'#C87A62','lipid8_F2':'#5B86B2','lipid8_F3':'#58A36F','nonlipid8_F1':'#64A99A','nonlipid8_F2':'#8B74B8','nonlipid8_F3':'#D9A93D'}
loading_map={'lipid8_F1':[('M_HDL_TG','M-HDL\nTG','F1'),('VLDL_size','VLDL\nsize','F1'),('MUFA','MUFA','F1'),('S_VLDL_TG','S-VLDL\nTG','F1')],'lipid8_F2':[('ApoA1','ApoA1','F2'),('HDL_CE','HDL\nCE','F2')],'lipid8_F3':[('XS_VLDL_FC','XS-VLDL\nFC','F3'),('VLDL_CE','VLDL\nCE','F3')],'nonlipid8_F1':[('Acetoacetate','Aceto-\nacetate','F1'),('bOHbutyrate','3-HB','F1')],'nonlipid8_F2':[('Val','Val','F2'),('Leu','Leu','F2'),('Phe','Phe','F2')],'nonlipid8_F3':[('Acetate','Acetate','F3'),('Glucose','Glucose','F3'),('Lactate','Lactate','F3')]}
lipid=pd.read_csv(lipid_res, sep='\t'); non=pd.read_csv(non_res, sep='\t')
for df in [lipid,non]: df['STD_All']=pd.to_numeric(df['STD_All'], errors='coerce')
loadvals={}; semcorr={}
for axis,items in loading_map.items():
    df=lipid if axis.startswith('lipid') else non
    for trait,label,factor in items:
        row=df[(df.lhs.eq(factor)) & (df.op.eq('=~')) & (df.rhs.eq(trait))]
        loadvals[(axis,trait)]=float(row.STD_All.iloc[0])
for df,axes in [(lipid,['lipid8_F1','lipid8_F2','lipid8_F3']),(non,['nonlipid8_F1','nonlipid8_F2','nonlipid8_F3'])]:
    f={axes[0]:'F1',axes[1]:'F2',axes[2]:'F3'}
    for i in range(3):
        for j in range(i+1,3):
            a,b=axes[i],axes[j]; f1,f2=f[a],f[b]
            row=df[(df.op.eq('~~')) & (((df.lhs.eq(f1)) & (df.rhs.eq(f2))) | ((df.lhs.eq(f2)) & (df.rhs.eq(f1))))]
            semcorr[(a,b)]=float(row.STD_All.iloc[0])
intv=pd.read_excel(xlsx, sheet_name='S12_IntFactorMetab_compact', header=3); intv['abs_rg']=pd.to_numeric(intv['abs_rg'], errors='coerce')
extv=pd.read_excel(xlsx, sheet_name='S18_ExtBivar_LDSC', header=3); extv['abs_rg']=pd.to_numeric(extv['abs_rg'], errors='coerce')
val=pd.merge(intv.groupby('factor',as_index=False)['abs_rg'].max().rename(columns={'abs_rg':'Internal'}), extv.groupby('factor',as_index=False)['abs_rg'].max().rename(columns={'abs_rg':'External'}), on='factor', how='outer').set_index('factor').reindex(axis_order)
ndd=pd.read_excel(xlsx, sheet_name='S23_NDD_LDSC_all18', header=3)
for c in ['rg','fdr_rg','p_rg']: ndd[c]=pd.to_numeric(ndd[c], errors='coerce')
ndd=ndd[ndd.trait1.isin(axis_order)].copy()
mix=pd.read_excel(xlsx, sheet_name='S25_MiXeR_focal_clean', header=3)
for c in ['dice','mixer_rho_beta']: mix[c]=pd.to_numeric(mix[c], errors='coerce')
mix=mix[mix.trait1.isin(axis_order)].copy().set_index('trait1').loc[['lipid8_F2','nonlipid8_F1','lipid8_F1']].reset_index(); mix['branch']=['HDL-core - AD','Ketone-body - PD','TG/VLDL - PD']
mpl.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9.2,'axes.labelsize':8,'xtick.labelsize':7.1,'ytick.labelsize':7.1,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':0.68})
fig=plt.figure(figsize=(11.8,7.05), dpi=300); gs=GridSpec(2,6,figure=fig,height_ratios=[1.38,1.0],width_ratios=[1,1,1,1,1,1],wspace=0.86,hspace=0.28)
axA=fig.add_subplot(gs[0,0:3]); axB=fig.add_subplot(gs[0,3:6])
def setup_path_ax(ax,title,panel):
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis('off'); ax.text(-0.055,1.01,panel,transform=ax.transAxes,fontsize=14,fontweight='bold',va='top'); ax.text(0.015,0.968,title,ha='left',va='top',fontsize=10.2,fontweight='bold')
def factor_node(ax,x,y,text,color):
    ax.scatter([x],[y],s=1780,c=[color],edgecolors='#222',linewidths=0.95,zorder=4); ax.text(x,y,text,ha='center',va='center',fontsize=7.75,color='white',fontweight='bold',zorder=5)
def obs_box(ax,x,y,text,width,height=0.084):
    ax.add_patch(FancyBboxPatch((x-width/2,y-height/2),width,height,boxstyle='round,pad=0.008,rounding_size=0.021',facecolor='white',edgecolor='#222',lw=0.72,zorder=4)); ax.text(x,y,text,ha='center',va='center',fontsize=5.95,zorder=5,linespacing=0.84)
def loading_arrow(ax,x1,y1,x2,y2,label):
    ax.add_patch(FancyArrowPatch((x1,y1-0.055),(x2,y2+0.045),arrowstyle='-|>',mutation_scale=7.8,lw=0.98,color='#252525',zorder=2,shrinkA=1,shrinkB=1)); mx=(x1+x2)/2; my=(y1+y2)/2+0.006; angle=np.degrees(np.arctan2((y2+0.045)-(y1-0.055),x2-x1)); angle=angle-180 if angle>90 else angle; angle=angle+180 if angle<-90 else angle; ax.text(mx,my,label,ha='center',va='center',fontsize=6.25,rotation=angle,bbox=dict(boxstyle='round,pad=0.045',facecolor='white',edgecolor='none',alpha=0.95),zorder=6)
def covariance_arc(ax,x1,x2,y,rad,label,dy=0.0):
    ax.add_patch(FancyArrowPatch((x1+0.072,y+0.068),(x2-0.072,y+0.068),connectionstyle=f'arc3,rad={rad}',arrowstyle='<->',mutation_scale=7.0,lw=0.72,color='#696969',zorder=1)); ax.text((x1+x2)/2,y+0.116+dy,label,ha='center',va='center',fontsize=6.05,bbox=dict(boxstyle='round,pad=0.045',facecolor='white',edgecolor='none',alpha=0.90),zorder=5)
def draw_model(ax,axes,xpos,offsets,title,panel):
    setup_path_ax(ax,title,panel); yF=0.655; yO=0.205
    for a1,a2,rad,dy in [(axes[0],axes[1],-0.18,0),(axes[1],axes[2],-0.18,0),(axes[0],axes[2],-0.30,0.055)]: covariance_arc(ax,xpos[a1],xpos[a2],yF,rad,f'{semcorr.get((a1,a2),semcorr.get((a2,a1))):.2f}',dy)
    for a in axes: factor_node(ax,xpos[a],yF,short_labels[a].replace(' ','\n'),colors[a])
    for a in axes:
        obs_w=0.083 if len(loading_map[a])>=4 else (0.094 if len(loading_map[a])==3 else 0.108)
        for (trait,label,factor),off in zip(loading_map[a],offsets[a]): obs_box(ax,xpos[a]+off,yO,label,obs_w); loading_arrow(ax,xpos[a],yF,xpos[a]+off,yO,f'{loadvals[(a,trait)]:.2f}')
draw_model(axA,['lipid8_F1','lipid8_F2','lipid8_F3'],{'lipid8_F1':0.235,'lipid8_F2':0.555,'lipid8_F3':0.825},{'lipid8_F1':[-0.140,-0.047,0.047,0.140],'lipid8_F2':[-0.058,0.058],'lipid8_F3':[-0.058,0.058]},'Lipid factor model','A')
draw_model(axB,['nonlipid8_F1','nonlipid8_F2','nonlipid8_F3'],{'nonlipid8_F1':0.205,'nonlipid8_F2':0.515,'nonlipid8_F3':0.820},{'nonlipid8_F1':[-0.062,0.062],'nonlipid8_F2':[-0.094,0,0.094],'nonlipid8_F3':[-0.094,0,0.094]},'Non-lipid factor model','B')
axC=fig.add_subplot(gs[1,0:2]); x=np.arange(len(axis_order)); w=0.35; axC.bar(x-w/2,val.Internal,w,color='#5A9F92',label='Internal'); axC.bar(x+w/2,val.External,w,color='#D9A93D',label='External'); axC.set_ylim(0,1.08); axC.set_ylabel('Max |rg|'); axC.set_xticks(x); axC.set_xticklabels([short_labels[a] for a in axis_order],rotation=30,ha='right'); axC.set_title('Metabolic identity validation',loc='left',fontweight='bold'); axC.legend(frameon=False,ncol=2,loc='upper left',fontsize=7,handlelength=1.5,columnspacing=1.4); axC.tick_params(axis='x', length=0); [axC.spines[s].set_visible(False) for s in ['top','right']]; axC.text(-0.17,1.08,'C',transform=axC.transAxes,fontsize=14,fontweight='bold',va='top')
axD=fig.add_subplot(gs[1,2:4]); diseases=['AD','PD','LBD']
for _,r in ndd.iterrows():
    xi=diseases.index(str(r.trait2)); yi=axis_order.index(str(r.trait1))
    evid=-math.log10(max(float(r.fdr_rg),1e-300))
    size=38+1450*abs(float(r.rg))
    sig=bool(r.get('fdr_significant_0.05',False))
    axD.scatter(xi,yi,s=size,c=[evid],cmap='YlOrRd',vmin=0,vmax=3.0,edgecolor='black' if sig else '#A8A8A8',linewidth=0.82 if sig else 0.32,zorder=3)
axD.set_xticks(range(3)); axD.set_xticklabels(diseases); axD.set_yticks(range(len(axis_order))); axD.set_yticklabels([short_labels[a] for a in axis_order]); axD.invert_yaxis(); axD.set_xlim(-0.5,2.82); axD.set_title('Factor-disease LDSC screen',loc='left',fontweight='bold'); [sp.set_visible(False) for sp in axD.spines.values()]; axD.tick_params(axis='both', length=0)
for rg_size,label in [(0.05,'0.05'),(0.10,'0.10'),(0.15,'0.15')]:
    axD.scatter([],[],s=38+1450*rg_size,c='white',edgecolor='#6F6F6F',linewidth=0.7,label=label)
leg=axD.legend(title='|rg|',frameon=False,loc='lower left',bbox_to_anchor=(0.00,-0.305),fontsize=6.35,title_fontsize=6.9,ncol=3,handletextpad=0.52,columnspacing=0.82,borderpad=0.0)
cax=axD.inset_axes([0.59,-0.285,0.37,0.044], transform=axD.transAxes)
sm=plt.cm.ScalarMappable(cmap='YlOrRd',norm=mpl.colors.Normalize(vmin=0,vmax=3.0))
cb=fig.colorbar(sm,cax=cax,orientation='horizontal')
cb.set_label('-log10(FDR)',labelpad=-1)
cb.ax.tick_params(labelsize=6.15,length=2,pad=1)
axD.text(-0.17,1.08,'D',transform=axD.transAxes,fontsize=14,fontweight='bold',va='top')
axE=fig.add_subplot(gs[1,4:6]); yy=np.arange(len(mix)); axE.barh(yy,mix.dice,color=['#5B86B2','#64A99A','#C77C5F'],height=0.46); axE.set_yticks(yy); axE.set_yticklabels(mix.branch); axE.invert_yaxis(); axE.set_xlim(0,0.86); axE.set_xlabel('Dice overlap'); axE.tick_params(axis='y', length=0); [axE.text(min(r.dice+0.025,0.76),i,f'rho={r.mixer_rho_beta:.2f}',va='center',fontsize=6.9) for i,r in mix.iterrows()]; axE.set_title('MiXeR overlap architecture',loc='left',fontweight='bold'); [axE.spines[s].set_visible(False) for s in ['top','right']]; axE.text(-0.17,1.08,'E',transform=axE.transAxes,fontsize=14,fontweight='bold',va='top')
plt.savefig(out_png,bbox_inches='tight',dpi=300); plt.savefig(out_pdf,bbox_inches='tight'); plt.savefig(out_svg,bbox_inches='tight'); print(out_png); print(out_pdf); print(out_svg)


