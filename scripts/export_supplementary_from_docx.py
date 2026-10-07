#!/usr/bin/env python3
"""Export only GSEA-related supplementary images and captions from Additional file 2.

Usage:
    python scripts/export_supplementary_from_docx.py Additional_file_2.docx

The DOCX is supplied locally by the caller and is not copied into the repository.
This script extracts the exact embedded PNG bytes for Figures S1, S2 and S10
(three pages), matching each image to its following caption through the DOCX
relationships. It does not regenerate the page layouts from statistical data.
The public provenance JSON maps exported panels to frozen repository sources.
Only the Python standard library is required.
"""
from pathlib import Path, PurePosixPath
from zipfile import ZipFile
from xml.etree import ElementTree as ET
import argparse, csv, hashlib, json, posixpath, re, struct

ROOT=Path(__file__).resolve().parents[1]
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'a':'http://schemas.openxmlformats.org/drawingml/2006/main',
    'r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships',
    'wp':'http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing'}
RELNS='http://schemas.openxmlformats.org/package/2006/relationships'
CAPTION=re.compile(r'^Supplementary Figure S(1|2|10)(?:\.| continued\.)')
WITHDRAWN={'NES_SD','NES_lower_95CI','NES_upper_95CI'}
SOURCE_MAP={
 ('S1',1):[{'repository_file':'results/Supp_Table_S2_v2.tsv',
     'selection':'bh_pool = hallmark_50pool; cohort_short = GSE51472; padj_within_pool < 0.05',
     'role':'NES and adjusted P values for the 20 significant discovery Hallmark sets'}],
 ('S2',1):[{'repository_file':'results/Figure1_panels/Figure1B_PCA_combat.pdf',
     'role':'Before/after ComBat PCA panel used by the published supplementary composite'}],
 ('S10',1):[
    {'repository_file':'results/Figure1_panels/Figure1A_PCA_percohort.pdf','panel':'a'},
    {'repository_file':'results/Figure1_panels/Figure1C_volcano_discovery.pdf','panel':'b'},
    {'repository_file':'results/Figure1_panels/Figure1D_volcano_validation.pdf','panel':'c'}],
 ('S10',2):[
    {'repository_file':'results/Figure1_panels/Figure1E_markers_discovery.pdf','panel':'d'},
    {'repository_file':'results/Figure1_panels/Figure1F_markers_validation.pdf','panel':'e'},
    {'repository_file':'results/Figure1_panels/Figure1G_crosscohort_forest.pdf','panel':'f'}],
 ('S10',3):[
    {'repository_file':'results/Figure1_panels/Figure1H_GO.pdf','panel':'g'},
    {'repository_file':'results/Figure1_panels/Figure1I_KEGG.pdf','panel':'h'}]
}

def sha(data):return hashlib.sha256(data).hexdigest()
def png_size(data):
    if data[:8]!=b'\x89PNG\r\n\x1a\n':raise ValueError('Expected an embedded PNG')
    return list(struct.unpack('>II',data[16:24]))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('docx',type=Path,help='Current Additional file 2 DOCX')
    parser.add_argument('--output',type=Path,default=ROOT/'figures/current/supplementary')
    args=parser.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
    source_bytes=args.docx.read_bytes();records=[];counts={'S1':0,'S2':0,'S10':0}
    with ZipFile(args.docx) as archive:
        document=ET.fromstring(archive.read('word/document.xml'))
        relationships=ET.fromstring(archive.read('word/_rels/document.xml.rels'))
        rels={r.attrib['Id']:r.attrib for r in relationships.findall('{'+RELNS+'}Relationship')}
        pending=[]
        for paragraph_index,item in enumerate(document.find('w:body',NS)):
            caption=''.join(t.text or '' for t in item.findall('.//w:t',NS)).strip()
            images=item.findall('.//a:blip',NS)
            if images:
                pending=[]
                for blip in images:
                    rid=blip.attrib.get('{'+NS['r']+'}embed')
                    if not rid:raise ValueError('Linked image found instead of embedded image')
                    crop=[x.attrib for x in item.findall('.//a:srcRect',NS)]
                    if any(any(int(v)!=0 for v in x.values()) for x in crop):
                        raise ValueError('Cropped image cannot be represented by whole-PNG extraction')
                    extent=[x.attrib for x in item.findall('.//wp:extent',NS)]
                    pending.append({'relationship_id':rid,'image_paragraph_index':paragraph_index,
                                    'crop_attributes':crop,'word_extent_EMU':extent})
            matched=CAPTION.match(caption)
            if matched:
                if len(pending)!=1:raise ValueError('Expected exactly one preceding image for '+caption[:35])
                number='S'+matched.group(1);counts[number]+=1;page=counts[number]
                if (number,page) not in SOURCE_MAP:raise ValueError('Unexpected supplementary figure/page')
                rec=dict(pending[0]);relationship=rels[rec['relationship_id']]
                if relationship.get('TargetMode')=='External':raise ValueError('External image is not supported')
                target=posixpath.normpath(posixpath.join('word',relationship['Target']))
                if not target.startswith('word/media/'):raise ValueError('Unexpected image relationship target')
                data=archive.read(target);size=png_size(data)
                stem='Supplementary_Figure_'+number+(('_page'+str(page)) if number=='S10' else '')
                image_file=stem+'.png';caption_file=stem+'_legend.txt'
                (out/image_file).write_bytes(data);(out/caption_file).write_text(caption+'\n',encoding='utf-8')
                source_map=[]
                for mapped in SOURCE_MAP[(number,page)]:
                    m=dict(mapped);path=ROOT/m['repository_file']
                    if not path.is_file():raise FileNotFoundError(m['repository_file'])
                    m['sha256']=sha(path.read_bytes());source_map.append(m)
                rec.update({'figure':number,'page':page,'image_file':image_file,
                    'caption_file':caption_file,'caption_paragraph_index':paragraph_index,
                    'docx_media_part':target,'dimensions_pixels':size,'sha256':sha(data),
                    'byte_identical_to_embedded_media':(out/image_file).read_bytes()==data,
                    'image_reencoded':False,'pixel_transform_applied':False,
                    'frozen_source_map':source_map})
                records.append(rec);pending=[]
            elif caption:
                pending=[]
    if counts!={'S1':1,'S2':1,'S10':3}:raise ValueError('Incomplete requested supplement selection: '+str(counts))
    combined='\n\n'.join((out/r['caption_file']).read_text().strip() for r in records if r['figure']=='S10')+'\n'
    (out/'Supplementary_Figure_S10_legend.txt').write_text(combined,encoding='utf-8')
    with (ROOT/'results/Supp_Table_S2_v2.tsv').open(encoding='utf-8-sig',newline='') as f:
        reader=csv.DictReader(f,delimiter='\t');fields=[x for x in reader.fieldnames if x not in WITHDRAWN]
        values=[{k:r[k] for k in fields} for r in reader if r['bh_pool']=='hallmark_50pool' and r['cohort_short']=='GSE51472' and float(r['padj_within_pool'])<.05]
    if len(values)!=20:raise ValueError('Expected 20 significant discovery Hallmark sets')
    values.sort(key=lambda r:float(r['NES']),reverse=True)
    with (out/'Supplementary_Figure_S1_data.tsv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(values)
    manifest={'source_document':{'file_name':args.docx.name,'sha256':sha(source_bytes),
       'role':'Current Additional file 2; supplied by the caller and not included in this repository'},
       'export_method':'Exact extraction of embedded PNG bytes through document.xml relationships; captions copied as text. No resizing, cropping or re-encoding.',
       'scope':'GSEA-related supplementary Figures S1, S2 and S10 only; five images in total.',
       'excluded_figures':['S3','S4','S5','S6','S7','S8','S9','S11'],
       'layout_reproduction':'These exports reproduce the embedded AF2 pixels exactly. Mapping to frozen source PDFs/TSV documents analytical provenance; this exporter does not regenerate the composite layouts from those sources.',
       'database_version_note':'S10h preserves the archived KEGG panel. A new online KEGG query may change pathway labels and ordering; the frozen panel is retained for the reported analysis.',
       'S1_data_export':{'file':'Supplementary_Figure_S1_data.tsv','rows':20,
          'order':'Descending NES','withdrawn_NES_interval_columns_included':False},
       'figures':records,
       'script_sha256':sha(Path(__file__).read_bytes())}
    (out/'export_provenance.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(out),'figures':counts,'images':len(records),'all_byte_identical':all(r['byte_identical_to_embedded_media'] for r in records)}))

if __name__=='__main__':main()
