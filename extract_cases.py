#!/usr/bin/env python
#!/usr/bin/env python

import pandas
import numpy
import matplotlib.pyplot as plt
import xmltodict
from pandas.plotting import table


domain_data_dir = "temp/domain_details/"
shared_dom_color= (0,0.608,0.62,0.8)
other_dom_color= (0.6,0.6,0.6,0.8)
def visualize(p,q,shared_domain,d,p_df, q_df):
    fig = plt.figure(figsize=[12,5])
    # text_ax = fig.add_subplot(211)
    # ax = fig.add_subplot(312,)
    fig, (p_table_ax,ax, q_table_ax) = plt.subplots(3, 1, gridspec_kw={'height_ratios': [5, 1 , 5]},figsize=[12,12])

    # p_table_ax = fig.add_subplot(311,frame_on = False)
    p_table_ax.axis('off')
    p_table_ax.xaxis.set_visible(False)  # hide the x axis
    p_table_ax.yaxis.set_visible(False)  # hide the y axis
    # p_table_ax.set_title('Single-domain protein (P): ' + p + " with drug(D): " + d)
    table(p_table_ax, p_df,loc='center')

    # q_table_ax = fig.add_subplot(313,frame_on = False)
    q_table_ax.axis('off')
    q_table_ax.xaxis.set_visible(False)  # hide the x axis
    q_table_ax.yaxis.set_visible(False)  # hide the y axis
    table(q_table_ax, q_df,loc='center')
    # q_table_ax.set_title('Protein (Q) with same domain: ' + q + " with drug(D): " + d)




    # text_ax.spines["right"].set_visible(False)
    # text_ax.spines["left"].set_visible(False)
    # text_ax.spines["top"].set_visible(False)
    # text_ax.spines["bottom"].set_visible(False)
    # text_ax.get_yaxis().set_visible(False)
    # text_ax.get_xaxis().set_visible(False)

    p_path = domain_data_dir+p+".xml"
    with open(p_path) as pf:
        p_dict = xmltodict.parse(pf.read())
    p_seq  = p_dict["pfam"]["entry"]["sequence"]["#text"]
    p_len  = len (p_seq)
    ax.hlines(2, 0, p_len, linewidth=2, color="grey")
    p_domains = p_dict["pfam"]["entry"]["matches"]["match"]
    if type(p_domains) != list:
        p_domains =[p_domains]
    for dom in p_domains:
        acc = dom["@accession"]
        # type = dom["@type"]
        begin = int(dom["location"]["@start"])
        end = int(dom["location"]["@end"])
        #if pfam_A
        if acc==shared_domain:
            col = shared_dom_color
        else:
            col = other_dom_color
        ax.hlines(2, begin, end, linewidth=10, color=col)

    q_path = domain_data_dir + q + ".xml"
    with open(q_path) as qf:
        q_dict = xmltodict.parse(qf.read())
    q_seq = q_dict["pfam"]["entry"]["sequence"]["#text"]
    q_len = len(q_seq)
    ax.hlines(1, 0, q_len, linewidth=2, color="grey")
    q_domains = q_dict["pfam"]["entry"]["matches"]["match"]
    if type(q_domains) != list:
        q_domains =[q_domains]
    for dom in q_domains:
        acc = dom["@accession"]
        # type = dom["@type"]
        begin = int(dom["location"]["@start"])
        end = int(dom["location"]["@end"])
        # if pfam_A
        if acc == shared_domain:
            col = shared_dom_color
        else:
            col = other_dom_color
        ax.hlines(1, begin, end, linewidth=10, color=col)

    h_rng = float(max(p_len,q_len))
    h_margin = h_rng/10
    ax.set_xlim(-h_margin, h_rng+h_margin)
    ax.set_ylim(0.5, 2.5)
    # ax.get_yaxis().set_visible(False)
    p_label = "P ["+p+"]"
    q_label  = "Q ["+q+"]"
    ax.set_yticks([1,2])
    ax.set_yticklabels([q_label,p_label])
    fig_fname ="outputs/promissings/"+p+"_"+d+"_"+q+".svg"
    ax.set_title("Domain: "+shared_domain)
    plt.tight_layout()
    fig.savefig(fig_fname)
    plt.close(fig)


class StringConverter(dict):
    def __contains__(self, item):
        return True

    def __getitem__(self, item):
        return str

    def get(self, default=None):
        return str


col_names = pandas.read_csv("outputs/result_drug_level.csv", sep = "\t", nrows=0).columns
types_dict={col: str for col in col_names}
results_df = pandas.read_csv("outputs/result_drug_level.csv",converters=StringConverter())
results_df.fillna("",inplace=True)
Ps = results_df["onedomain-protein"]
Ms = results_df["domain"]
Ds = results_df["interacting_drug"]
Qs = results_df["neg"]
num_Qs  = results_df["num_neg"].astype(int).values

useful_rows = numpy.nonzero(num_Qs>0)[0]

import pandas as pd



col_names = pandas.read_csv("data/BindingDB_All_2021m0.tsv/BindingDB_All.tsv", sep = "\t", nrows=0).columns
# types_dict = {"Ki (nM)": float,"Kd (nM)": float,"IC50 (nM)": float,"EC50 (nM)": float}
# types_dict.update({col: str for col in col_names if col not in types_dict})
types_dict={col: str for col in col_names}

df = pandas.read_csv("data/BindingDB_All_2021m0.tsv/BindingDB_All.tsv", sep = "\t",error_bad_lines=False,converters=StringConverter())
useful_cols = ["Ki (nM)","Kd (nM)","IC50 (nM)","EC50 (nM)","kon (M-1-s-1)","koff (s-1)","pH","Temp (C)"]
mdfile = "# Potential examples for problem 1:\n"



query_table_dict = dict()

for idx, row in df.iterrows():
    uniprot_id = row["UniProt (SwissProt) Primary ID of Target Chain"]
    pubchem_cid = row["PubChem CID"]
    if (uniprot_id,pubchem_cid) in query_table_dict:
        query_table_dict[(uniprot_id,pubchem_cid)].append(idx)
    else:
        query_table_dict[(uniprot_id, pubchem_cid)] = [idx]





for i in useful_rows:
    if ";" in Qs[i]:
        current_Qs = Qs[i].split(";")
    else:
        current_Qs = [Qs[i]]
    current_Qs = [x.strip() for x in current_Qs]
    p = Ps[i]
    m = Ms[i]
    d = Ds[i]
    for q in current_Qs:

        mdfile += "## Interaction between domain "+m+" and drug "+d+"\n\n"
        # p_rows_mask = (df["UniProt (SwissProt) Primary ID of Target Chain"] == p) & (df["PubChem CID"] == d)
        # q_rows_mask = (df["UniProt (SwissProt) Primary ID of Target Chain"] == q) & (df["PubChem CID"] == d)
        #
        # p_rows = df.loc[p_rows_mask, useful_cols].copy()
        # q_rows = df.loc[q_rows_mask, useful_cols].copy()
        p_rows_idx  = query_table_dict[(p,d)]
        p_rows = df.loc[p_rows_idx,useful_cols]

        q_rows_idx  = query_table_dict[(q,d)]
        q_rows = df.loc[p_rows_idx,useful_cols]

        # mdfile += "Single-domain protein (P) interacting with the drug: " + p +"\n\n"
        # mdfile += p_rows.to_markdown(index = False)+"\n\n"
        #
        #
        # mdfile += "Another (Q) protein with the same domain: " + q +"\n\n"
        # mdfile += q_rows.to_markdown(index = False)+"\n\n"

        # mdfile += "![]("+p+"_"+q+".svg)\n\n"
        visualize(p, q, m,d,p_rows,q_rows)


