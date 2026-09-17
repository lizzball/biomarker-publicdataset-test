from pathlib import Path
import json,numpy as np,pandas as pd,plotly.express as px,plotly.graph_objects as go,streamlit as st
ROOT=Path(__file__).resolve().parent; ART=ROOT/"artifacts"
st.set_page_config(page_title="TED Biomarker Signal Lab",page_icon="◉",layout="wide")
st.markdown("""<style>.block-container{max-width:1320px;padding-top:2rem}.hero{background:#211B21;color:#fff;padding:2.6rem;border-radius:28px}.hero h1{font-family:Georgia,serif;font-size:4rem;font-weight:500;letter-spacing:-.04em}.summary{background:#fff;border-left:5px solid #6E1F2A;padding:1rem 1.2rem;margin:1rem 0}[data-testid=stMetric]{background:#fff;border:1px solid #ddd;padding:15px;border-radius:16px}</style>""",unsafe_allow_html=True)
st.markdown("""<div class='hero'><small>PUBLIC RARE-DISEASE BIOMARKER MODEL · GSE58331</small><h1>TED Biomarker Signal Lab</h1><p>Nested-CV gene-expression model with threshold sensitivity, biomarker stability, subject-level uncertainty and QC review.</p></div>""",unsafe_allow_html=True)
rp=ART/"model_report.json"
if not rp.exists():
    st.warning("Training artifacts are not present yet. Run the public-data workflow or GitHub Action."); st.stop()
r=json.loads(rp.read_text()); p=r["primary_validation"]
pred=pd.read_csv(ART/"oof_subject_predictions.csv"); roc=pd.read_csv(ART/"roc_curve.csv"); thr=pd.read_csv(ART/"threshold_metrics.csv"); stab=pd.read_csv(ART/"feature_stability.csv"); rep=pd.read_csv(ART/"cv_repeat_metrics.csv"); qc=pd.read_csv(ART/"subject_qc.csv"); pcs=pd.read_csv(ART/"subject_pca.csv")
st.markdown("""<div class='summary'><b>Analysis summary.</b> The model holds anatomy constant by comparing TED versus normal anterior-orbit tissue. Repeated samples are collapsed to biological subject before validation, and gene selection occurs inside each training fold. The dataset is more than a decade old and is used here as a reproducible methods case study, not as evidence for a current clinical biomarker.</div>""",unsafe_allow_html=True)
a,b,c,d=st.columns(4); a.metric("Biological subjects",p["n_subjects"]); b.metric("OOF AUROC",f"{p['oof_auroc']:.3f}"); c.metric("Bootstrap 95% CI",f"{p['bootstrap_auroc_95ci'][0]:.2f}–{p['bootstrap_auroc_95ci'][1]:.2f}"); d.metric("Avg precision",f"{p['average_precision']:.3f}")
tabs=st.tabs(["Model validation","Decision layer","Biomarker stability","Prediction stability","QC & outliers","Data & limitations"])
with tabs[0]:
    l,rcol=st.columns([1.1,1])
    with l:
        f=go.Figure([go.Scatter(x=roc.fpr,y=roc.tpr,mode="lines",name="OOF"),go.Scatter(x=[0,1],y=[0,1],mode="lines",line=dict(dash="dash"),name="Chance")]); f.update_layout(xaxis_title="False-positive rate",yaxis_title="True-positive rate",height=450); st.plotly_chart(f,use_container_width=True)
    with rcol:
        f=px.box(rep,y=["auroc","average_precision"],points="all"); f.update_layout(height=330); st.plotly_chart(f,use_container_width=True); st.caption(p["note"])
with tabs[1]:
    t=st.slider("Decision threshold",.05,.95,.50,.05); row=thr.iloc[(thr.threshold-t).abs().argsort()[:1]].iloc[0]
    c1,c2,c3,c4,c5=st.columns(5); c1.metric("TED detected",int(row.tp)); c2.metric("Controls flagged",int(row.fp)); c3.metric("Sensitivity",f"{row.sensitivity:.1%}"); c4.metric("Specificity",f"{row.specificity:.1%}"); c5.metric("Balanced accuracy",f"{row.balanced_accuracy:.1%}")
    f=px.line(thr,x="threshold",y=["sensitivity","specificity","ppv","npv"],markers=True); f.add_vline(x=t,line_dash="dash"); st.plotly_chart(f,use_container_width=True)
    st.info(f"Evaluation cohort remains fixed at {p['n_subjects']} biological subjects ({p['ted_n']} TED / {p['control_n']} controls). The slider changes only the decision rule.")
with tabs[2]:
    n=st.slider("Genes to display",5,30,15); d=stab.head(n).copy(); d["direction"]=np.where(d.median_coefficient>=0,"Higher TED score","Lower TED score")
    st.plotly_chart(px.bar(d.sort_values("selection_frequency"),x="selection_frequency",y="gene",orientation="h",color="direction",hover_data=["median_coefficient","selection_count"]),use_container_width=True)
    st.caption("Selection frequency is a stability measure under resampling, not evidence of causal disease biology.")
with tabs[3]:
    pred["correct_at_0_5"]=((pred.score_mean>=.5).astype(int)==pred.y_true)
    st.plotly_chart(px.scatter(pred,x="score_mean",y="subject_id",color="class_label",symbol="correct_at_0_5",error_x="score_sd",hover_data=["source_titles","n_source_samples","score_min","score_max"]),use_container_width=True)
    st.caption("Error bars show resampling instability across repeated held-out predictions.")
with tabs[4]:
    st.plotly_chart(px.scatter(pcs,x="PC1",y="PC2",color="class_label",text="subject_id"),use_container_width=True)
    st.dataframe(qc.sort_values("qc_flag",ascending=False),use_container_width=True,hide_index=True)
    st.info("QC flags are reviewed, not automatically deleted. Extreme biology can itself appear as a statistical outlier.")
with tabs[5]:
    st.markdown("""**Cohort path:** 175 total GEO samples → anterior-orbit TED/control subset → repeated arrays collapsed to biological subject → nested-CV model.

**Age:** GSE58331 was submitted in 2014 and released in 2015. Larger, newer and independent cohorts are required for modern biomarker claims.

**Specimen limitation:** anterior-orbit tissue is a research/surgical specimen, not a practical screening matrix.

**Public field:** `VALUE` is submitter-processed log2 RMA expression from GPL570 arrays.

**Boundary:** this is a portfolio demonstration of computational biomarker validation, not a clinical diagnostic.""")
