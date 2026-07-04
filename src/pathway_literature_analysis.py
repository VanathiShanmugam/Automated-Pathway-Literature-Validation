
import os
import time
import pandas as pd
import requests
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ================= USER SETTINGS =================

INPUT_FILE = "data/sample_pathways.txt"
EMAIL = "your_email@example.com"
NCBI_API_KEY = ""

ALPHA = 0.05
MIN_HITS = 1
SLEEP_TIME = 0.2

# ================= FILE READER =================

def read_file(fp):

    ext = os.path.splitext(fp)[1].lower()

    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(fp)

    for enc in ["utf-8", "latin1"]:
        try:
            return pd.read_csv(fp, sep="\t", encoding=enc)
        except:
            pass

    return pd.read_csv(fp)

def detect_pathway_column(df):

    for col in ["Term", "Pathway", "TERM", "term"]:
        if col in df.columns:
            return col

    return df.columns[0]

def detect_disease(fp):

    fn = os.path.basename(fp).lower()

    if any(x in fn for x in ["leuko", "leukoplakia", "opmd", "premalignant"]):
        return "LEUKOPLAKIA"

    if any(x in fn for x in ["oscc", "oral_cancer", "scc"]):
        return "OSCC"

    return "LEUKOPLAKIA"

def pubmed_count(query):

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": 0,
        "email": EMAIL
    }

    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    try:
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()
        return int(r.json()["esearchresult"]["count"])
    except Exception as e:
        print("PubMed error:", e)
        return 0

# ================= MAIN =================

disease = detect_disease(INPUT_FILE)

if disease == "LEUKOPLAKIA":
    disease_query = (
        '"oral leukoplakia"[Title/Abstract] OR '
        '"oral potentially malignant disorder"[Title/Abstract] OR '
        '"Leukoplakia, Oral"[Mesh]'
    )
else:
    disease_query = (
        '"oral squamous cell carcinoma"[Title/Abstract] OR '
        '"oral cancer"[Title/Abstract] OR '
        'OSCC[Title/Abstract] OR '
        '"Mouth Neoplasms"[Mesh]'
    )

print("Disease:", disease)

df = read_file(INPUT_FILE)
pathway_col = detect_pathway_column(df)

print("Pathway column:", pathway_col)

pathways = (
    df[pathway_col]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
)

print("Total pathways:", len(pathways))

TOTAL_PUBMED_ARTICLES = pubmed_count("all[sb]")
disease_hits = pubmed_count(disease_query)

results = []

for i, pathway in enumerate(pathways, start=1):

    print(f"[{i}/{len(pathways)}] {pathway}")

    pathway_hits = pubmed_count(f'"{pathway}"')

    combined_hits = pubmed_count(
        f'("{pathway}") AND ({disease_query})'
    )

    a = combined_hits
    b = max(pathway_hits - combined_hits, 0)
    c = max(disease_hits - combined_hits, 0)
    d = max(TOTAL_PUBMED_ARTICLES - a - b - c, 0)

    try:
        odds_ratio, pvalue = fisher_exact(
            [[a, b], [c, d]],
            alternative="greater"
        )
    except:
        odds_ratio = 0
        pvalue = 1

    results.append([
        pathway,
        pathway_hits,
        disease_hits,
        combined_hits,
        odds_ratio,
        pvalue
    ])

    time.sleep(SLEEP_TIME)

results_df = pd.DataFrame(
    results,
    columns=[
        "Pathway",
        "Pathway_Hits",
        "Disease_Hits",
        "Disease_Pathway_Hits",
        "Literature_Odds_Ratio",
        "P_Value"
    ]
)

results_df["FDR"] = multipletests(
    results_df["P_Value"],
    method="fdr_bh"
)[1]

results_df.sort_values(
    ["FDR", "Literature_Odds_Ratio"],
    ascending=[True, False],
    inplace=True
)

prefix = disease.upper()

results_df.to_csv(
    f"{prefix}_All_Pathway_Literature_Associations.csv",
    index=False
)

sig_df = results_df[results_df["FDR"] < ALPHA]

sig_df.to_csv(
    f"{prefix}_Significant_Pathway_Literature_Associations.csv",
    index=False
)

supported_df = results_df[
    results_df["Disease_Pathway_Hits"] >= MIN_HITS
]

supported_df.to_csv(
    f"{prefix}_Literature_Supported_Pathways.csv",
    index=False
)

print("\\n================================")
print("ANALYSIS COMPLETE")
print("================================")
print(f"Disease: {prefix}")
print(f"Total pathways: {len(results_df)}")
print(f"Significant pathways: {len(sig_df)}")
print(f"Literature-supported pathways: {len(supported_df)}")
