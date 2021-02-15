#!/usr/bin/env python


import matplotlib.pyplot as plt
import xmltodict



domain_data_dir = "temp/domain_details/"
shared_dom_col= (0,0.608,0.62,0.8)
other_dom_col= (0.6,0.6,0.6,0.8)
def visualize(p,q,shared_domain):
    fig = plt.figure(figsize=[8, 2])
    # text_ax = fig.add_subplot(211)
    ax = fig.add_subplot(111)
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
            col = shared_dom_col
        else:
            col = other_dom_col
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
            col = shared_dom_col
        else:
            col = other_dom_col
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
    fig_fname ="outputs/promissings/"+p+"_"+q+".svg"
    fig.suptitle("Domain: "+shared_domain)
    fig.savefig(fig_fname)
    plt.close(fig)


class StringConverter(dict):
    def __contains__(self, item):
        return True

    def __getitem__(self, item):
        return str

    def get(self, default=None):
        return str

import pandas
import numpy
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
        visualize(p,q,m)
        mdfile += "## Interaction between domain "+m+" and drug "+d+"\n\n"
        p_rows_mask = (df["UniProt (SwissProt) Primary ID of Target Chain"] == p) & (df["PubChem CID"] == d)
        q_rows_mask = (df["UniProt (SwissProt) Primary ID of Target Chain"] == q) & (df["PubChem CID"] == d)

        p_rows = df.loc[p_rows_mask, useful_cols].copy()
        q_rows = df.loc[q_rows_mask, useful_cols].copy()


        mdfile += "Single-domain protein (P) interacting with the drug: " + p +"\n\n"
        mdfile += p_rows.to_markdown(index = False)+"\n\n"


        mdfile += "Another (Q) protein with the same domain: " + q +"\n\n"
        mdfile += q_rows.to_markdown(index = False)+"\n\n"

        mdfile += "![]("+p+"_"+q+".svg)\n\n"

with open("outputs/promissings/doc.md", "w") as outf:
    outf.writelines(mdfile)

