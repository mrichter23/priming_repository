# app.py
from shiny import App, ui, render, reactive
#from shiny.types import ImgData
import pandas as pd
import matplotlib.pyplot as plt
from itertools import chain
#from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import io
from pyodide.http import open_url


def parse_uploaded_file(file):
    if file is None:
        return []

    # Shiny returns a list of uploaded files (even for single upload)
    file_info = file[0]

    df = pd.read_csv(file_info["datapath"], header=None)
    return df.iloc[:, 0].dropna().astype(str).tolist()


def build_membership_dict(lists_dict):
    all_items = set(chain.from_iterable(lists_dict.values()))
    membership = {}
    for item in all_items:
        membership[item] = {k: (item in v) for k, v in lists_dict.items()}
    return pd.DataFrame.from_dict(membership, orient="index")


def plot_upset(membership_df):
    # Custom UpSet-like plot (no upsetplot dependency)
    df = membership_df.astype(bool).copy()

    # Count intersections
    grouped = df.groupby(list(df.columns)).size().reset_index(name="count")

    # Sort by count descending
    grouped = grouped.sort_values("count", ascending=False)

    combos = grouped.drop(columns=["count"]) 
    counts = grouped["count"].values

    n_combos = len(grouped)
    n_sets = df.shape[1]

    fig, (ax_bar, ax_matrix) = plt.subplots(
        2, 1,
        figsize=(max(6, n_combos * 0.6), 8),
        gridspec_kw={"height_ratios": [3, 1]}
    )

    # --- Top bar chart ---
    ax_bar.bar(range(n_combos), counts, color="grey")
    ax_bar.set_ylabel("Intersection size")
    ax_bar.set_xticks([])

    # --- Bottom matrix ---
    for i in range(n_combos):
        for j, col in enumerate(combos.columns):
            if combos.iloc[i, j]:
                ax_matrix.scatter(i, j, s=80, color="black")

    # Improve margins so dots are not cut off
    ax_matrix.set_ylim(-0.5, n_sets - 0.5)
    ax_matrix.set_xlim(-0.5, n_combos - 0.5)

    ax_matrix.set_yticks(range(n_sets))
    ax_matrix.set_yticklabels(combos.columns)
    ax_matrix.set_xlabel("Intersections")

    plt.subplots_adjust(hspace=0.05, left=0.2)

    return fig


# ---------- UI ----------

# Predefined lists per category

rna_e14_ctx_neurons = pd.read_csv("www/DEGs_E14_Ctx_up_neurons_vs_NSC.txt", header=None, on_bad_lines='skip')
rna_e14_lge_neurons = pd.read_csv(open_url("www/DEGs_E14_LGE_up_neurons_vs_NSC.txt"), header=None, on_bad_lines='skip')
astro = pd.read_csv(open_url("www/high_expressed_genes_astrocytes_Pereira_1.csv"), header=None, on_bad_lines='skip')

atac_ctx_nsc_e14 = pd.read_csv(open_url("www/DARs_NSC_Ctx_up_E14_vs_E18.txt"), header=None, on_bad_lines='skip')
atac_ctx_nsc_e18 = pd.read_csv(open_url("www/DARs_NSC_Ctx_up_E18_vs_E14.txt"), header=None, on_bad_lines='skip')
atac_e14_nsc_ctx = pd.read_csv(open_url("www/DARs_NSC_E14_up_Ctx_vs_LGE.txt"), header=None, on_bad_lines='skip')
atac_e14_nsc_lge = pd.read_csv(open_url("www/DARs_NSC_E14_up_LGE_vs_Ctx.txt"), header=None, on_bad_lines='skip')
atac_e18_nsc_ctx = pd.read_csv(open_url("www/DARs_NSC_E18_up_Ctx_vs_LGE.txt"), header=None, on_bad_lines='skip')
atac_e18_nsc_lge = pd.read_csv(open_url("www/DARs_NSC_E18_up_LGE_vs_Ctx.txt"), header=None, on_bad_lines='skip')
atac_lge_nsc_e14 = pd.read_csv(open_url("www/DARs_NSC_LGE_up_E14_vs_E18.txt"), header=None, on_bad_lines='skip')
atac_lge_nsc_e18 = pd.read_csv(open_url("www/DARs_NSC_LGE_up_E18_vs_E14.txt"), header=None, on_bad_lines='skip')

motifs_ctx_nsc_e14 = pd.read_csv(open_url("www/Motifs_NSC_Ctx_up_E14_vs_E18.txt"), header=None, on_bad_lines='skip')
motifs_ctx_nsc_e18 = pd.read_csv(open_url("www/Motifs_NSC_Ctx_up_E18_vs_E14.txt"), header=None, on_bad_lines='skip')
motifs_e14_nsc_ctx = pd.read_csv(open_url("www/DARs_NSC_E14_up_Ctx_vs_LGE.txt"), header=None, on_bad_lines='skip')
motifs_e14_nsc_lge = pd.read_csv(open_url("www/DARs_NSC_E14_up_LGE_vs_Ctx.txt"), header=None, on_bad_lines='skip')
motifs_e18_nsc_ctx = pd.read_csv(open_url("www/DARs_NSC_E18_up_Ctx_vs_LGE.txt"), header=None, on_bad_lines='skip')
motifs_e18_nsc_lge = pd.read_csv(open_url("www/DARs_NSC_E18_up_LGE_vs_Ctx.txt"), header=None, on_bad_lines='skip')
motifs_lge_nsc_e14 = pd.read_csv(open_url("www/DARs_NSC_LGE_up_E14_vs_E18.txt"), header=None, on_bad_lines='skip')
motifs_lge_nsc_e18 = pd.read_csv(open_url("www/DARs_NSC_LGE_up_E18_vs_E14.txt"), header=None, on_bad_lines='skip')

cutnrun = pd.read_csv(open_url("www/targets_2IR.csv"), header=None, on_bad_lines='skip')

predefined_lists = {
    1: {
        "DEGs E14 Ctx up Neurons vs NSC": rna_e14_ctx_neurons[0].tolist(),
        "DEGs E14 LGE up Neurons vs NSC": rna_e14_lge_neurons[0].tolist(),
        "Genes highly expressed in Astrocytes": astro[0].tolist()
    },
    2: {
        "DARs Ctx NSC up E14 vs E18": atac_ctx_nsc_e14[0].tolist(),
        "DARs Ctx NSC up E18 vs E14": atac_ctx_nsc_e18[0].tolist(),
        "DARs E14 NSC up Ctx vs LGE": atac_e14_nsc_ctx[0].tolist(),
        "DARs E14 NSC up LGE vs Ctx": atac_e14_nsc_lge[0].tolist(),
        "DARs E18 NSC up Ctx vs LGE": atac_e18_nsc_ctx[0].tolist(),
        "DARs E18 NSC up LGE vs Ctx": atac_e18_nsc_lge[0].tolist(),
        "DARs LGE NSC up E14 vs E18": atac_lge_nsc_e14[0].tolist(),
        "DARs LGE NSC up E18 vs E14": atac_lge_nsc_e18[0].tolist()
    },
    3: {
        "Motifs Ctx NSC up E14 vs E18": motifs_ctx_nsc_e14[0].tolist(),
        "Motifs Ctx NSC up E18 vs E14": motifs_ctx_nsc_e18[0].tolist(),
        "Motifs E14 NSC up Ctx vs LGE": motifs_e14_nsc_ctx[0].tolist(),
        "Motifs E14 NSC up LGE vs Ctx": motifs_e14_nsc_lge[0].tolist(),
        "Motifs E18 NSC up Ctx vs LGE": motifs_e18_nsc_ctx[0].tolist(),
        "Motifs E18 NSC up LGE vs Ctx": motifs_e18_nsc_lge[0].tolist(),
        "Motifs LGE NSC up E14 vs E18": motifs_lge_nsc_e14[0].tolist(),
        "Motifs LGE NSC up E18 vs E14": motifs_lge_nsc_e18[0].tolist(),
    },
    4: {
        "TGIF2 cut&run": cutnrun[0].tolist()
    }
}


# Custom category names
category_names = {
    1: "DEGs",
    2: "Genes in proximity to DARs",
    3: "Enriched motifs",
    4: "Additional (e.g. cut&run)"
}


def list_input_ui(id):
    return ui.panel_well(
        ui.h4(category_names.get(id, f"Category {id}")),
        ui.input_selectize(
            f"preset_{id}",
            "Select predefined list(s)",
            choices=list(predefined_lists[id].keys()),
            multiple=True
        ),
        ui.input_file(f"file_{id}", "Or upload CSV (one column)")
    )


app_ui = ui.page_fluid(
    ui.h2("Investigating priming by differential expression and accessibility during brain development"),

    ui.div(
    {"style": "margin:0; padding:0; height:auto;"},
    
    #ui.output_image("image")
    
    ui.tags.img(src="https://github.com/mrichter23/priming_repository/tree/main/myapp/www/overview.png", height="20%")
    ),
    
    ui.p(
        "In stem cells, the primed expression of differentiation-related genes enables their progeny to rapidly initiate lineage-specific programs."
        " In our paper, we investigated transcriptional priming during neurogenesis"
        " to elucidate its role in neural stem cell (NSC) maintenance and differentiation."
        " ATAC-seq and RNA-seq of NSCs and neurons across forebrain regions allowed us to identify"
        " hundreds of neuronal differentiation genes associated to differentially open regions" 
        " while still being lowly expressed in NSCs. This is consistent with a transcriptionally primed but inactive state."
        " Our database allows determining primed gene expression across regions and stages. To that end, this app"
        " flexibly allows users to explore our differential analyses, motif searches and cut&run, or upload their own genes of interest."
        " For more information, please see our manuscript: https://doi.org/10.1101/2025.02.13.63595"
    ),

    ui.row(
        ui.column(3, list_input_ui(1)),
        ui.column(3, list_input_ui(2)),
        ui.column(3, list_input_ui(3)),
        ui.column(3, list_input_ui(4)),
    ),
    ui.hr(),
    ui.output_plot("upset_plot"),
    ui.download_button("download_plot", "Download Plot"),
    ui.hr(),
    ui.h3("Items present in ALL selected lists"),
    ui.div(
        {"style": "max-height:300px; overflow-y:auto; border:1px solid #ccc; padding:5px;"},
        ui.output_table("intersection_table")
    ),
    ui.download_button("download_table", "Download Table")
)


# ---------- Server ----------

def server(input, output, session):
    #@render.image
    #def image():
    #    img_url = "https://github.com/mrichter23/priming_repository/tree/main/myapp/overview.png"
# ---------- Helper functions ----------
    #    #from pathlib import Path
    #    img_url = "https://github.com/mrichter23/priming_repository/tree/main/myapp/overview.png"
    #    return ui.tags.img(src=img_url, height="100%", width="100%")
    #    #dir = Path(__file__).resolve().parent
    #    #img: ImgData = {"src": str(dir/"overview.png"), "width": "500px", "style": "display:block; margin:0; padding:0;"}
    #    #return {"src":img_url}
    #    #return {open_url("www/overview.png")}
    #return img_url

    @reactive.Calc
    def selected_lists():
        lists = {}
        for i in range(1, 5):
            preset_list = input[f"preset_{i}"]()
            uploaded = input[f"file_{i}"]()

            items = []

            # Combine multiple predefined lists if selected
            if preset_list:
                for p in preset_list:
                    items.extend(predefined_lists[i][p])

            # If file uploaded, extend (not overwrite)
            if uploaded is not None:
                items.extend(parse_uploaded_file(uploaded))

            if items:
                # Use more informative names for lists
                selected_names = []
                if preset_list:
                    selected_names.extend(preset_list)
                if uploaded is not None:
                    selected_names.append(f"Upload_{i}")

                label = " + ".join(selected_names) if selected_names else f"List_{i}"
                lists[label] = list(set(items))

        return lists


    @reactive.Calc
    def membership_df():
        lists = selected_lists()
        if len(lists) < 2:
            return None
        return build_membership_dict(lists)


    @output
    @render.plot
    def upset_plot():
        df = membership_df()
        if df is None:
            return None
        fig = plot_upset(df)
        return fig


    @output
    @render.table
    def intersection_table():
        lists = selected_lists()
        if len(lists) < 2:
            return pd.DataFrame()

        common = set.intersection(*map(set, lists.values()))
        return pd.DataFrame({"Items in all lists": list(common)})


    @output
    @render.download(filename="upset_plot.png")
    def download_plot():
        df = membership_df()
        if df is None:
            return None
        fig = plot_upset(df)
        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        buf.seek(0)
        return buf


    @output
    @render.download(filename="intersection.csv")
    def download_table():
        lists = selected_lists()
        if len(lists) < 2:
            return None

        common = set.intersection(*map(set, lists.values()))
        df = pd.DataFrame({"Items": list(common)})
        buf = io.StringIO()
        df.to_csv(buf, index=False)
        buf.seek(0)
        return buf


# ---------- Run App ----------

app = App(app_ui, server)


# Requirements:
# pip install shiny pandas matplotlib upsetplot
