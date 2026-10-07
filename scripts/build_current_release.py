#!/usr/bin/env python3
"""Export current GSEA reporting and rebuild Figure 1 from the frozen v2 results.

This entry point does not rerun expression preprocessing or fgsea. It verifies
and reuses the committed NES, P values and leading edges, independently checks
BH correction, omits withdrawn historical interval columns, and redraws the
submission figure using the original publication layout.

Requirements: Python >= 3.10; reportlab, pypdf, Pillow; Poppler pdftoppm for
raster exports. Supply --font-dir for the original Arial fonts. On macOS these
are detected in /System/Library/Fonts/Supplemental; Liberation Sans is supported
as a portable substitute (the build manifest records the actual fonts).

Run from any directory: python scripts/build_current_release.py
Use --data-only to export and verify tables without graphics dependencies.
"""
from pathlib import Path
from io import BytesIO
from collections import defaultdict
import argparse, csv, hashlib, json, math, shutil, subprocess, platform
from importlib.metadata import version
import xml.sax.saxutils as xml

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-root', type=Path, default=ROOT,
                    help='Root for results/current and figures/current outputs')
parser.add_argument('--font-dir', type=Path,
                    help='Directory containing Arial or LiberationSans font files')
parser.add_argument('--pdftoppm', type=Path,
                    help='Poppler pdftoppm executable; default: search PATH')
parser.add_argument('--data-only', action='store_true',
                    help='Verify and export data only; do not render figures')
args = parser.parse_args()
RELEASE_ROOT = args.output_root.resolve()
DATA = RELEASE_ROOT/'results/current'
OUT = RELEASE_ROOT/'figures/current'
QA = DATA/'figure1_QA'
for p in [DATA, OUT]: p.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT/'results/Supp_Table_S2_v2.tsv'
CONSISTENCY = ROOT/'results/cross_cohort_consistency_v2.tsv'
PATHWAYS = [
 'HALLMARK_MTORC1_SIGNALING', 'CUSTOM_LEPTIN_MTOR_SIGNALING',
 'CUSTOM_MTOR_AUTOPHAGY_LIPOPHAGY_AXIS', 'CUSTOM_LIPOPHAGY_CORE',
 'CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION']
LABELS = ['Hallmark mTORC1 signaling','Leptin–mTOR signaling',
 'mTOR–autophagy–lipophagy axis','Lipophagy core','VIC osteogenic differentiation']
COHORTS = ['GSE51472','INTEGRATED','GSE83453']
GENES = ['RUNX2','SPP1','IBSP','ENPP1','DLX5','SLC20A1','WNT5A','COL1A1','COL1A2','ALPL']
WITHDRAWN_COLUMNS = {'NES_SD','NES_lower_95CI','NES_upper_95CI'}
sha256 = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

def read_tsv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        reader=csv.DictReader(f, delimiter='\t')
        return list(reader), reader.fieldnames

def write_table(path, rows, fields, delimiter='\t', bom=False):
    with path.open('w', encoding='utf-8-sig' if bom else 'utf-8', newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter=delimiter,
                         lineterminator='\r\n' if bom else '\n')
        w.writeheader();w.writerows(rows)

def export_data():
    rows, fields=read_tsv(SOURCE)
    if len(rows)!=68: raise ValueError('Expected 68 frozen GSEA results')
    keys=[(r['pathway'],r['cohort_short'],r['bh_pool']) for r in rows]
    if len(set(keys))!=68: raise ValueError('Duplicate frozen result key')
    groups=defaultdict(list)
    for r in rows:
        if int(r['leading_edge_size'])!=len(r['leading_edge'].split(',')):
            raise ValueError('Leading-edge size mismatch: '+r['pathway'])
        for k in ['NES','pval','padj_within_pool']:
            if not math.isfinite(float(r[k])): raise ValueError('Nonfinite '+k)
        if not 0 <= float(r['pval']) <= 1: raise ValueError('Invalid P value')
        groups[(r['cohort_short'],r['bh_pool'])].append(r)
    expected={(co,'hypothesis_5pool'):5 for co in COHORTS}
    expected.update({(co,'supplementary_solo'):1 for co in COHORTS})
    expected[('GSE51472','hallmark_50pool')]=50
    if {k:len(v) for k,v in groups.items()}!=expected:
        raise ValueError('Analysis-set or testing-pool membership mismatch')
    verification=[]
    for (co,pool), group in groups.items():
        ordered=sorted(group,key=lambda r:float(r['pval']));n=len(ordered)
        bh=[min(1.0,float(r['pval'])*n/(i+1)) for i,r in enumerate(ordered)]
        for i in range(n-2,-1,-1): bh[i]=min(bh[i],bh[i+1])
        err=max(abs(x-float(r['padj_within_pool'])) for x,r in zip(bh,ordered))
        if err>1e-12: raise ValueError('BH verification failed: '+co+' '+pool)
        verification.append({'analysis_set':co,'pool':pool,'n_tests':n,
          'maximum_absolute_error':err,'passed':True})
    clean_fields=[k for k in fields if k not in WITHDRAWN_COLUMNS]
    clean=[{k:r[k] for k in clean_fields} for r in rows]
    write_table(DATA/'GSEA_results.tsv',clean,clean_fields)
    full={(r['pathway'],r['cohort_short']):r for r in rows if r['bh_pool']=='hypothesis_5pool'}
    classes={r['pathway']:r['classification'] for r in read_tsv(CONSISTENCY)[0]}
    figure=[]
    for path,label in zip(PATHWAYS,LABELS):
        rr=[full[(path,co)] for co in COHORTS]
        concordant=len({float(r['NES'])>0 for r in rr})==1
        cls=('Strong consistency' if sum(float(r['padj_within_pool'])<.05 for r in rr)>=2 else 'Directional consistency') if concordant else 'Inconsistent'
        if classes[path]!=cls: raise ValueError('Consistency class mismatch: '+path)
        for co,r in zip(COHORTS,rr):
            figure.append(dict(pathway=path,pathway_label=label,cohort=co,
                N_input=r['N_input'],N_mapped=r['N_mapped'],NES=r['NES'],
                raw_P=r['pval'],BH_adjusted_P=r['padj_within_pool'],
                significant_BH_0_05=float(r['padj_within_pool'])<.05,
                consistency=cls,leading_edge_size=r['leading_edge_size'],
                leading_edge=r['leading_edge']))
    write_table(DATA/'Figure1_GSEA_values.csv',figure,list(figure[0]),',',True)
    leading={co:set(full[('CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION',co)]['leading_edge'].split(',')) for co in COHORTS}
    if set.union(*leading.values())!=set(GENES): raise ValueError('Leading-edge union differs from the declared figure layout')
    shared=set.intersection(*leading.values())
    if shared!=set(GENES[:6]): raise ValueError('Shared leading-edge core differs')
    matrix=[dict(gene=g,**{co:int(g in leading[co]) for co in COHORTS},shared_all_three=g in shared) for g in GENES]
    write_table(DATA/'Figure1_osteogenic_leading_edge_matrix.csv',matrix,list(matrix[0]),',',True)
    bh_report={'historical_results_recomputed':False,'BH_recalculated':True,
       'all_pools_passed':True,'pool_checks':verification,
       'withdrawn_columns_omitted':sorted(WITHDRAWN_COLUMNS),
       'full_result_rows':len(rows),'figure_rows':len(figure),'leading_edge_matrix_cells':30}
    (DATA/'BH_verification.json').write_text(json.dumps(bh_report,indent=2)+'\n')
    plot=[]
    for r in figure:
        t=dict(r)
        for k in ['NES','raw_P','BH_adjusted_P']:t[k]=float(t[k])
        for k in ['N_input','N_mapped','leading_edge_size']:t[k]=int(t[k])
        t['leading_edge']=t['leading_edge'].split(',');plot.append(t)
    return dict(pathway_order=PATHWAYS,cohort_order=COHORTS,results=plot,
                osteogenic_leading_edge={'shared_genes':GENES[:6]}),bh_report

def setup_graphics():
    global canvas,pdfmetrics,TTFont,HexColor,Paragraph,ParagraphStyle
    global PdfReader,PdfWriter,Transformation,Image,FONT_FILES,POPPLER,SVG_FONT_FAMILY
    from reportlab.pdfgen import canvas
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.lib.colors import HexColor
    from reportlab import rl_config
    from reportlab.platypus import Paragraph
    from reportlab.lib.styles import ParagraphStyle
    from pypdf import PdfReader,PdfWriter,Transformation
    from PIL import Image
    rl_config.useA85=False
    font_specs=[('Arial',['Arial.ttf','Arial Bold.ttf','Arial Italic.ttf','Arial Bold Italic.ttf']),
      ('Liberation Sans',['LiberationSans-Regular.ttf','LiberationSans-Bold.ttf','LiberationSans-Italic.ttf','LiberationSans-BoldItalic.ttf'])]
    dirs=[args.font_dir] if args.font_dir else [Path('/System/Library/Fonts/Supplemental'),Path('/usr/share/fonts/truetype/msttcorefonts'),Path('/usr/share/fonts/truetype/liberation2'),Path('/usr/share/fonts/truetype/liberation')]
    found=None
    for d in dirs:
        for family,names in font_specs:
            if all((d/n).is_file() for n in names):found=(family,[d/n for n in names]);break
        if found:break
    if not found: raise SystemExit('Arial or Liberation Sans fonts not found; supply --font-dir.')
    family,paths=found
    SVG_FONT_FAMILY='Arial, Helvetica, sans-serif' if family=='Arial' else 'Liberation Sans, Arial, Helvetica, sans-serif'
    FONT_FILES=[{'file':p.name,'family':family,'sha256':sha256(p)} for p in paths]
    for name,path in zip(['Arial','ArialBold','ArialItalic','ArialBoldItalic'],paths):
        pdfmetrics.registerFont(TTFont(name,str(path)))
    pdfmetrics.registerFontFamily('Arial',normal='Arial',bold='ArialBold',italic='ArialItalic',boldItalic='ArialBoldItalic')
    POPPLER=str(args.pdftoppm) if args.pdftoppm else shutil.which('pdftoppm')
    if not POPPLER:raise SystemExit('Poppler pdftoppm not found; install Poppler or supply --pdftoppm.')
    QA.mkdir(parents=True,exist_ok=True)

W,H,MM=170.,178.,72/25.4
INK='#252A2E'; MUTED='#64717A'; LINE='#CBD2D6'; ACCENT='#367994'
COHORT_COLOR=['#367994','#AE7159','#577E6A']
COHORT_FILL=['#AFCBD8','#E3B3A0','#BFD0C5']
COHORT_TINT=['#F0F6F8','#FAF2EE','#F1F6F3']

TITLE='Figure 1. Cross-dataset GSEA identifies consistent enrichment of the VIC osteogenic program.'
CANONICAL_LEGEND='Figure 1. Cross-dataset GSEA identifies consistent enrichment of the VIC osteogenic program. a Human valve transcriptomic datasets: GSE51472 (discovery; 5 normal and 5 calcified valves), GSE51472 + GSE12644 (ComBat-adjusted sensitivity analysis; 15 normal and 15 calcified valves), and GSE83453 (validation; 8 normal and 9 calcified tricuspid valves). The sensitivity analysis includes the discovery samples. Sclerotic valves in GSE51472 and bicuspid valves in GSE83453 were excluded. b Normalized enrichment scores (NES) for five prespecified gene sets. Positive NES indicates enrichment in calcified valves. Symbols distinguish analysis sets; filled symbols indicate Benjamini–Hochberg-adjusted P < 0.05. Adjustment was performed across the five sets within each analysis set. Strong consistency denotes concordant NES signs with adjusted P < 0.05 in at least two analysis sets; directional consistency denotes concordant signs with fewer than two significant analysis sets; inconsistent denotes opposing signs. c Membership of the ten genes in the union of osteogenic leading-edge subsets. Six genes were shared across all three subsets. Exact results and gene membership are provided in Additional file 3: Tables S1 and S2.\n'
LEGEND=CANONICAL_LEGEND[len(TITLE):].strip()

class Drawing:
    def __init__(self):
        self.c=canvas.Canvas(str(OUT/'Figure1.pdf'),pagesize=(W*MM,H*MM),pageCompression=1,initialFontName='Arial',invariant=1)
        self.c.setTitle(TITLE);self.c.setAuthor('CAVD study')
        self.svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{H}mm" viewBox="0 0 {W*MM} {H*MM}">','<rect width="100%" height="100%" fill="white"/>']
        self.bounds=[];self.points=[];self.members=[]
    def text(self,s,x,y,size=7,align='left',bold=False,color=INK,italic=False):
        font='Arial'+('Bold' if bold else '')+('Italic' if italic else '')
        tw=pdfmetrics.stringWidth(s,font,size)/MM
        dx=0 if align=='left' else -tw if align=='right' else -tw/2
        self.c.setFillColor(HexColor(color));self.c.setFont(font,size);self.c.drawString((x+dx)*MM,(H-y)*MM,s)
        anchor={'left':'start','right':'end','center':'middle'}[align]
        self.svg.append(f'<text x="{x*MM}" y="{y*MM}" text-anchor="{anchor}" font-family="{SVG_FONT_FAMILY}" font-size="{size}" font-weight="{"bold" if bold else "normal"}" font-style="{"italic" if italic else "normal"}" fill="{color}">{xml.escape(s)}</text>')
        self.bounds.append({'text':s,'box':[x+dx,y-size/MM,x+dx+tw,y+.65]})
    def line(self,x1,y1,x2,y2,color=INK,lw=.6,dash=None):
        self.c.saveState();self.c.setStrokeColor(HexColor(color));self.c.setLineWidth(lw)
        if dash:self.c.setDash(dash)
        self.c.line(x1*MM,(H-y1)*MM,x2*MM,(H-y2)*MM);self.c.restoreState()
        ds=f' stroke-dasharray="{",".join(map(str,dash))}"' if dash else ''
        self.svg.append(f'<line x1="{x1*MM}" y1="{y1*MM}" x2="{x2*MM}" y2="{y2*MM}" stroke="{color}" stroke-width="{lw}"{ds}/>')
    def rect(self,x,y,w,h,fill='white',stroke=LINE,lw=.4):
        self.c.setFillColor(HexColor(fill));self.c.setStrokeColor(HexColor(stroke));self.c.setLineWidth(lw)
        self.c.rect(x*MM,(H-y-h)*MM,w*MM,h*MM,fill=1,stroke=1)
        self.svg.append(f'<rect x="{x*MM}" y="{y*MM}" width="{w*MM}" height="{h*MM}" fill="{fill}" stroke="{stroke}" stroke-width="{lw}"/>')
    def symbol(self,x,y,cohort,r=.91,filled=True,color=None):
        color=color or COHORT_COLOR[cohort]; fill=color if filled else '#FFFFFF'
        self.c.setFillColor(HexColor(fill));self.c.setStrokeColor(HexColor(color));self.c.setLineWidth(.7)
        if cohort==0:
            self.c.circle(x*MM,(H-y)*MM,r*MM,fill=1,stroke=1)
            self.svg.append(f'<circle cx="{x*MM}" cy="{y*MM}" r="{r*MM}" fill="{fill}" stroke="{color}" stroke-width=".7"/>')
        elif cohort==1:self.rect(x-r,y-r,2*r,2*r,fill,color,.7)
        else:
            pts=[(x,y-r*1.16),(x-r*1.04,y+r*.92),(x+r*1.04,y+r*.92)]
            p=self.c.beginPath();p.moveTo(pts[0][0]*MM,(H-pts[0][1])*MM)
            for a,b in pts[1:]:p.lineTo(a*MM,(H-b)*MM)
            p.close();self.c.drawPath(p,fill=1,stroke=1)
            self.svg.append(f'<polygon points="{" ".join(f"{a*MM},{b*MM}" for a,b in pts)}" fill="{fill}" stroke="{color}" stroke-width=".7"/>')
    def finish(self):
        self.c.showPage();self.c.save();self.svg.append('</svg>')
        (OUT/'Figure1_editable.svg').write_text('\n'.join(self.svg))
        bad=[b for b in self.bounds if b['box'][0]<0 or b['box'][1]<0 or b['box'][2]>W or b['box'][3]>H]
        (QA/'text_bounds_QA.json').write_text(json.dumps({'out_of_page':bad,'text_count':len(self.bounds),'all_text_bounds':self.bounds},indent=2))
        (QA/'plotted_elements.json').write_text(json.dumps({'NES_points':self.points,'membership_cells':self.members},indent=2))
        if bad:raise ValueError(bad)

def panel(d,letter,title,y):
    d.text(letter,1.6,y,11,bold=True);d.text(title,9,y,8.2,bold=True)

def design(d,q):
    panel(d,'a','Human valve transcriptomes',5.3)
    cards=[(2,48,'Discovery','GSE51472','5 normal | 5 calcified','Affymetrix U133 Plus 2','5 sclerotic excluded'),
           (59,63,'Integrated sensitivity','GSE51472 + GSE12644','15 normal | 15 calcified','Affymetrix; ComBat-adjusted','Includes discovery samples'),
           (127,41,'Validation','GSE83453','8 normal | 9 calcified','Illumina HumanHT-12 v4','Tricuspid valves')]
    for i,(x,w,role,accession,counts,platform,note) in enumerate(cards):
        d.rect(x,9,w,27,COHORT_TINT[i],LINE,.45)
        d.line(x,9,x+w,9,COHORT_COLOR[i],1.2)
        cx=x+w/2
        d.text(role,cx,14.1,7.7,'center',True,COHORT_COLOR[i])
        d.text(accession,cx,19.0,7.2,'center',True)
        d.text(counts,cx,24.0,7,'center')
        d.text(platform,cx,28.7,6.7,'center',color=MUTED)
        d.text(note,cx,33.3,6.7,'center',color=MUTED)
    d.line(50.5,20.9,58,20.9,MUTED,.65)
    d.line(58,20.9,56.6,20,MUTED,.65);d.line(58,20.9,56.6,21.8,MUTED,.65)

def p_text(p):
    if p<.001:
        mant,exp=f'{p:.2e}'.split('e')
        sup=str(int(exp)).translate(str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
        return mant+' × 10'+sup
    return f'{p:.{max(0,2-math.floor(math.log10(p)))}f}'

def draw_p(d,p,x,y,bold):
    if p>=.001:
        d.text(p_text(p),x,y,7,'center',bold=bold)
        return
    mant,exp=f'{p:.2e}'.split('e')
    base=mant+' × 10'; exp=str(int(exp)).replace('-','−')
    font='ArialBold' if bold else 'Arial'
    wb=pdfmetrics.stringWidth(base,font,7)/MM
    we=pdfmetrics.stringWidth(exp,font,5.2)/MM
    x0=x-(wb+we)/2
    d.text(base,x0,y,7,bold=bold)
    d.text(exp,x0+wb,y-1.05,5.2,bold=bold)

def enrichment(d,q):
    panel(d,'b','Prespecified gene-set enrichment',43)
    for j,(x,label) in enumerate([(53,'Discovery'),(88,'Integrated sensitivity'),(143,'Validation')]):
        d.symbol(x,47.8,j,.8);d.text(label,x+2.1,48.65,7)
    d.text('Gene set',2,55,7.2,bold=True)
    d.text('BH-adjusted P',127.7,55,7.2,'center',True)
    d.text('Consistency',155,55,7.2,'center',True)
    x,w=49.5,61.0
    X=lambda v:x+(v+2.1)/4.2*w
    for i in range(5):
        if i%2==0:d.rect(2,57+12*i,166,12,'#F5F7F8','#F5F7F8',.35)
    for v in [-2,-1,0,1,2]:
        d.line(X(v),57.5,X(v),117,'#9DA8AF' if v==0 else '#DFE4E7',.6 if v==0 else .35,[2,2] if v==0 else None)
    labels=[['Hallmark mTORC1 signaling'],['Leptin–mTOR signaling'],['mTOR–autophagy–', 'lipophagy axis'],['Lipophagy core'],['VIC osteogenic','differentiation']]
    lookup={(r['pathway'],r['cohort']):r for r in q['results']}
    for i,path in enumerate(q['pathway_order']):
        yc=63+12*i
        for k,label in enumerate(labels[i]):d.text(label,2.7,yc+.85+(k-(len(labels[i])-1)/2)*3.2,7.2,bold=i==4)
        cls=lookup[(path,q['cohort_order'][0])]['consistency'].replace(' consistency','')
        d.text(cls,155,yc+.85,7,'center',bold=cls=='Strong',color=ACCENT if cls=='Strong' else INK)
        for j,cohort in enumerate(q['cohort_order']):
            r=lookup[(path,cohort)]; yy=yc+(j-1)*3.5;xx=X(r['NES']);sig=r['BH_adjusted_P']<.05
            d.symbol(xx,yy,j,1.,sig)
            label=p_text(r['BH_adjusted_P'])
            draw_p(d,r['BH_adjusted_P'],127.7,yy+.87,sig)
            d.points.append({'pathway':path,'cohort':cohort,'NES':r['NES'],'BH_adjusted_P':r['BH_adjusted_P'],'P_display':label,'display_xy_mm':[xx,yy],'filled':sig,'symbol':['circle','square','triangle'][j],'color':COHORT_COLOR[j],'consistency':cls})
    d.line(x,117,x+w,117,lw=.65)
    for v in [-2,-1,0,1,2]:
        d.line(X(v),117,X(v),118,lw=.65);d.text(str(v).replace('-','−'),X(v),121,7,'center')
    d.text('NES (calcified versus normal)',x+w/2,125.1,7.2,'center')
    d.symbol(46,130,0,.78,True,INK);d.text('BH P < 0.05',48,130.85,7)
    d.symbol(83,130,0,.78,False,INK);d.text('BH P ≥ 0.05',85,130.85,7)

def leading_edge(d,q):
    panel(d,'c','Osteogenic leading-edge genes',138)
    le=q['osteogenic_leading_edge']
    # Preserve a declared order: all shared genes first, then all other union genes.
    genes=['RUNX2','SPP1','IBSP','ENPP1','DLX5','SLC20A1','WNT5A','COL1A1','COL1A2','ALPL']
    lookup={(r['pathway'],r['cohort']):r for r in q['results']}
    shared=set(le['shared_genes'])
    d.rect(42.5,141.7,75,23.8,'#F0F6F8','#F0F6F8',.35)
    for k,g in enumerate(genes):
        cx=48.75+12.5*k
        d.text(g,cx,146.2,6.9,'center',bold=g in shared,italic=True)
    for j,(cohort,label) in enumerate(zip(q['cohort_order'],['Discovery','Integrated sensitivity','Validation'])):
        yy=151.2+j*6
        d.symbol(3,yy,j,.77);d.text(label,5.3,yy+.9,7)
        present=lookup[('CUSTOM_VIC_OSTEOGENIC_DIFFERENTIATION',cohort)]['leading_edge']
        for k,g in enumerate(genes):
            cx=48.75+12.5*k;hit=g in present
            d.rect(cx-5.3,yy-2.1,10.6,4.2,COHORT_FILL[j] if hit else '#FFFFFF','#B7C3CA',.4)
            d.members.append({'gene':g,'cohort':cohort,'member':hit,'display_xy_mm':[cx,yy]})
    d.line(43.45,167.0,116.55,167.0,ACCENT,.65)
    d.line(43.45,167.0,43.45,166.1,ACCENT,.65);d.line(116.55,167.0,116.55,166.1,ACCENT,.65)
    d.text('Shared in all three analysis sets',80,170.8,7,'center',color=ACCENT)
    d.text('Colored cell: present in leading edge',2,176,6.7,color=MUTED)
    d.text('White cell: absent',168,176,6.7,'right',color=MUTED)

def proof():
    style=ParagraphStyle('caption',fontName='Arial',fontSize=8,leading=9.6,textColor=HexColor(INK))
    body=xml.escape(LEGEND)
    for letter,start in [('a','Human'),('b','Normalized'),('c','Membership')]:body=body.replace(letter+' '+start,'<b>'+letter+'</b> '+start)
    p=Paragraph('<b>'+xml.escape(TITLE)+'</b> '+body,style)
    _,ph=p.wrap(166*MM,100*MM);required=H+3+ph/MM
    if required>225:raise ValueError(f'Figure and legend need {required:.2f} mm, exceeding 225 mm')
    buf=BytesIO();c=canvas.Canvas(buf,pagesize=(W*MM,225*MM),initialFontName='Arial',invariant=1)
    p.drawOn(c,2*MM,(225-H-3)*MM-ph);c.save();buf.seek(0)
    pg=PdfReader(buf).pages[0];pg.merge_transformed_page(PdfReader(OUT/'Figure1.pdf').pages[0],Transformation().translate(0,(225-H)*MM))
    wr=PdfWriter();wr.add_page(pg);wr.add_metadata({'/Title':TITLE,'/Author':'CAVD study'})
    with (OUT/'Figure1_with_legend_proof.pdf').open('wb') as f:wr.write(f)
    (QA/'legend_size_QA.json').write_text(json.dumps({'figure_mm':[W,H],'proof_mm':[W,225],'legend_height_mm':ph/MM,'combined_height_mm':required,'title_word_count':len(TITLE.split())-2,'legend_word_count':len(LEGEND.split())},indent=2))

def render():
    pp=POPPLER
    subprocess.run([pp,'-png','-r','600','-singlefile',str(OUT/'Figure1.pdf'),str(OUT/'Figure1_600dpi')],check=True)
    with Image.open(OUT/'Figure1_600dpi.png') as im:
        im.convert('RGB').save(OUT/'Figure1_600dpi.tif',compression='tiff_lzw',dpi=(600,600))
        im.thumbnail((1600,1800));im.save(OUT/'Figure1_preview.png',dpi=(240,240))
    subprocess.run([pp,'-png','-r','220','-singlefile',str(OUT/'Figure1_with_legend_proof.pdf'),str(OUT/'Figure1_with_legend_preview')],check=True)

def main():
    q,bh=export_data()
    (OUT/'Figure1_legend.txt').write_text(CANONICAL_LEGEND)
    if not args.data_only:
        setup_graphics()
        d=Drawing();design(d,{});enrichment(d,q);leading_edge(d,q);d.finish();proof();render()
    outputs=[]
    roots=[DATA] if args.data_only else [DATA,OUT]
    for directory in roots:
        for p in sorted(directory.rglob('*')):
            if p.is_file() and p.name!='build_manifest.json':
                outputs.append({'path':str(p.relative_to(RELEASE_ROOT)),'bytes':p.stat().st_size,'sha256':sha256(p)})
    manifest={'source_run':'gsea-run-v2','source_run_commit':'18f0f7f43289a0751ed17e0d6bcbcb495c0f04c6',
       'reporting_version':'v125','figure_layout_source':'Original Figure 1 publication rebuild script, 2026-10-06',
       'expression_analysis_rerun':False,'NES_intervals_reported':False,
       'data_only':args.data_only,
       'software':{'python':platform.python_version(), 'packages':{} if args.data_only else {name:version(name) for name in ['reportlab','pypdf','Pillow']}, 'pdftoppm':None if args.data_only else subprocess.run([POPPLER,'-v'],capture_output=True,text=True,check=True).stderr.splitlines()[0]},
       'input_files':[{'path':str(p.relative_to(ROOT)),'sha256':sha256(p)} for p in [SOURCE,CONSISTENCY]],
       'script_sha256':sha256(Path(__file__)),
       'fonts':[] if args.data_only else FONT_FILES,
       'figure_mm':[W,H],'raster_dpi':600,'outputs':outputs}
    (DATA/'build_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'results':'results/current','figure':'figures/current',
       'result_rows':68,'figure_rows':15,'BH_verified':True,'data_only':args.data_only}))

if __name__=='__main__':main()
