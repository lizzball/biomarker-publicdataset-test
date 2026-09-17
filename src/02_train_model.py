from pathlib import Path
import json,joblib,numpy as np,pandas as pd
from sklearn.feature_selection import SelectKBest,f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,average_precision_score,brier_score_loss,roc_curve,confusion_matrix,balanced_accuracy_score
from sklearn.model_selection import GridSearchCV,RepeatedStratifiedKFold,StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.decomposition import PCA

ROOT=Path(__file__).resolve().parents[1]; PROC=ROOT/"data/processed"; ART=ROOT/"artifacts"; ART.mkdir(exist_ok=True)

def load(slug):
    X=pd.read_csv(PROC/f"{slug}_subject_expression.tsv",sep="\t",index_col=0)
    m=pd.read_csv(PROC/f"{slug}_subject_metadata.csv").set_index("subject_id")
    X=X.loc[m.index]; y=(m.class_label=="TED").astype(int)
    return X.loc[:,X.var()>1e-8],y,m

def pipe():
    return Pipeline([("imp",SimpleImputer(strategy="median")),("sel",SelectKBest(f_classif)),("scale",StandardScaler()),("clf",LogisticRegression(max_iter=5000,class_weight="balanced",solver="liblinear",random_state=42))])

grid={"sel__k":[5,10,20,40,80],"clf__C":[.01,.1,1,10]}
X,y,meta=load("anterior")

outer=RepeatedStratifiedKFold(n_splits=5,n_repeats=10,random_state=42)
rows=[]; feats=[]
for i,(tr,te) in enumerate(outer.split(X,y)):
    inner=StratifiedKFold(3,shuffle=True,random_state=100+i)
    g=GridSearchCV(pipe(),grid,scoring="roc_auc",cv=inner,n_jobs=-1).fit(X.iloc[tr],y.iloc[tr])
    pr=g.best_estimator_.predict_proba(X.iloc[te])[:,1]
    for sid,yy,pp in zip(X.iloc[te].index,y.iloc[te],pr): rows.append({"subject_id":sid,"repeat":i//5,"fold":i%5,"y_true":int(yy),"score":float(pp)})
    sel=g.best_estimator_.named_steps["sel"]; genes=X.columns[sel.get_support()]
    co=g.best_estimator_.named_steps["clf"].coef_[0]
    for gene,c in zip(genes,co): feats.append({"gene":gene,"coefficient":float(c),"repeat":i//5,"fold":i%5})

pred=pd.DataFrame(rows)
sub=pred.groupby("subject_id").agg(y_true=("y_true","first"),score_mean=("score","mean"),score_sd=("score","std"),score_min=("score","min"),score_max=("score","max"),n_oof_predictions=("score","size")).reset_index()
sub=sub.merge(meta.reset_index(),on="subject_id")
sub.to_csv(ART/"oof_subject_predictions.csv",index=False); pred.to_csv(ART/"oof_all_predictions.csv",index=False)

repeat=[]
for r,d in pred.groupby("repeat"):
    repeat.append({"repeat":r,"auroc":roc_auc_score(d.y_true,d.score),"average_precision":average_precision_score(d.y_true,d.score),"brier":brier_score_loss(d.y_true,d.score)})
repeat=pd.DataFrame(repeat); repeat.to_csv(ART/"cv_repeat_metrics.csv",index=False)

stab=pd.DataFrame(feats).groupby("gene").agg(selection_count=("gene","size"),median_coefficient=("coefficient","median")).reset_index()
stab["selection_frequency"]=stab.selection_count/50
stab=stab.sort_values(["selection_frequency","selection_count"],ascending=False); stab.to_csv(ART/"feature_stability.csv",index=False)

yv=sub.y_true.to_numpy(); pv=sub.score_mean.to_numpy()
fpr,tpr,tt=roc_curve(yv,pv); pd.DataFrame({"fpr":fpr,"tpr":tpr,"threshold":tt}).to_csv(ART/"roc_curve.csv",index=False)
thr=[]
for t in np.arange(.05,.951,.05):
    z=(pv>=t).astype(int); tn,fp,fn,tp=confusion_matrix(yv,z,labels=[0,1]).ravel()
    thr.append({"threshold":round(float(t),2),"tn":tn,"fp":fp,"fn":fn,"tp":tp,"sensitivity":tp/(tp+fn),"specificity":tn/(tn+fp),"ppv":tp/(tp+fp) if tp+fp else np.nan,"npv":tn/(tn+fn) if tn+fn else np.nan,"balanced_accuracy":balanced_accuracy_score(yv,z)})
pd.DataFrame(thr).to_csv(ART/"threshold_metrics.csv",index=False)

Z=StandardScaler(with_std=False).fit_transform(SimpleImputer(strategy="median").fit_transform(X))
pcs=PCA(n_components=2,random_state=42).fit_transform(Z)
q=pd.DataFrame({"subject_id":X.index,"PC1":pcs[:,0],"PC2":pcs[:,1]}).merge(meta.reset_index(),on="subject_id")
q["distance"]=np.sqrt((q.PC1-q.PC1.median())**2+(q.PC2-q.PC2.median())**2); cut=q.distance.quantile(.99); q["qc_flag"]=q.distance>cut
q.to_csv(ART/"subject_pca.csv",index=False); q[["subject_id","class_label","n_source_samples","distance","qc_flag"]].to_csv(ART/"subject_qc.csv",index=False)

rng=np.random.default_rng(42); boot=[]
for _ in range(2000):
    ix=rng.integers(0,len(yv),len(yv))
    if len(np.unique(yv[ix]))==2: boot.append(roc_auc_score(yv[ix],pv[ix]))
lo,hi=np.quantile(boot,[.025,.975])

report={"primary_validation":{"dataset":"GSE58331","tissue":"Anterior Orbit","n_subjects":len(sub),"ted_n":int(yv.sum()),"control_n":int((1-yv).sum()),"oof_auroc":roc_auc_score(yv,pv),"bootstrap_auroc_95ci":[float(lo),float(hi)],"average_precision":average_precision_score(yv,pv),"brier":brier_score_loss(yv,pv),"repeat_auroc_median":float(repeat.auroc.median()),"repeat_auroc_min":float(repeat.auroc.min()),"repeat_auroc_max":float(repeat.auroc.max()),"qc_flagged_subjects":q.loc[q.qc_flag,"subject_id"].tolist(),"note":"Research tissue classifier; score is not a calibrated clinical probability."}}
(ART/"model_report.json").write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
