from pathlib import Path
import json,re,time
import GEOparse, numpy as np, pandas as pd

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/"data/raw"; PROC=ROOT/"data/processed"; DOCS=ROOT/"docs"
RAW.mkdir(parents=True,exist_ok=True); PROC.mkdir(parents=True,exist_ok=True); DOCS.mkdir(parents=True,exist_ok=True)

def get_geo(acc):
    err=None
    for i in range(3):
        try: return GEOparse.get_GEO(geo=acc,destdir=str(RAW),silent=False)
        except Exception as e: err=e; time.sleep(10*(i+1))
    raise err

def first(x): return x[0] if isinstance(x,list) and x else x

gse=get_geo("GSE58331")
expr=gse.pivot_samples("VALUE").apply(pd.to_numeric,errors="coerce")
meta=pd.DataFrame([{
    "sample_id":sid,
    "title":first(g.metadata.get("title","")),
    "source_name":first(g.metadata.get("source_name_ch1","")),
    "platform_id":first(g.metadata.get("platform_id","GPL570"))
} for sid,g in gse.gsms.items()])

gpl=gse.gpls.get("GPL570")
if gpl is None:
    obj=get_geo("GPL570"); gpl=obj.gpls["GPL570"] if hasattr(obj,"gpls") else obj
ann=gpl.table.copy()
if "ID" in ann.columns: ann=ann.set_index("ID",drop=False)
ann.index=ann.index.astype(str)
sym=next(c for c in ann.columns if "gene symbol" in str(c).lower() or str(c).lower() in {"symbol","gene_symbol"})
common=expr.index.astype(str).intersection(ann.index)
expr=expr.loc[common]; symbols=ann.loc[common,sym].astype(str).str.split(" /// ").str[0].replace({"---":np.nan,"nan":np.nan})
keep=symbols.notna()
gene=expr.loc[keep].copy()
gene.insert(0,"gene_symbol",symbols.loc[keep].values)
gene=gene.groupby("gene_symbol").median(numeric_only=True)

def build(tissue,slug):
    m=meta[meta.title.str.startswith((f"TED {tissue}",f"Normal {tissue}"),na=False)].copy()
    m["class_label"]=np.where(m.title.str.startswith("TED "),"TED","Control")
    pat=re.compile(rf"^(TED|Normal) {re.escape(tissue)} (\d+)")
    m["subject_id"]=[(lambda z: f"{'TED' if z.group(1)=='TED' else 'Control'}_{z.group(2)}" if z else None)(pat.search(t)) for t in m.title]
    m=m.dropna(subset=["subject_id"])
    ids=[s for s in m.sample_id if s in gene.columns]; m=m[m.sample_id.isin(ids)]
    Xs=gene[ids].T; Xs.index.name="sample_id"
    merged=Xs.join(m.set_index("sample_id")[["subject_id","class_label","title"]])
    genes=[c for c in merged.columns if c not in {"subject_id","class_label","title"}]
    X=merged.groupby("subject_id")[genes].mean()
    sm=merged.groupby("subject_id").agg(class_label=("class_label","first"),n_source_samples=("class_label","size"),source_titles=("title",lambda x:" | ".join(x))).reset_index()
    X.to_csv(PROC/f"{slug}_subject_expression.tsv",sep="\t")
    sm.to_csv(PROC/f"{slug}_subject_metadata.csv",index=False)
    Xs.to_csv(PROC/f"{slug}_sample_expression.tsv",sep="\t")
    m[["sample_id","subject_id","class_label","title"]].to_csv(PROC/f"{slug}_sample_metadata.csv",index=False)
    return {"arrays":len(m),"subjects":len(sm),"ted_subjects":int((sm.class_label=="TED").sum()),"control_subjects":int((sm.class_label=="Control").sum()),"genes":X.shape[1]}

counts={"anterior":build("Anterior Orbit","anterior"),"lacrimal":build("Lacrimal gland","lacrimal")}
meta.to_csv(DOCS/"geo_sample_inventory.csv",index=False)
lineage={"accession":"GSE58331","platform":"GPL570","value":"submitter-processed log2 RMA","all_series_samples":int(expr.shape[1]),"probe_rows":int(expr.shape[0]),"counts":counts}
(PROC/"data_lineage.json").write_text(json.dumps(lineage,indent=2)); (DOCS/"data_lineage.json").write_text(json.dumps(lineage,indent=2))
print(json.dumps(lineage,indent=2))
