from dash import html, dcc, register_page, Input, Output, State, callback
import dash_bootstrap_components as dbc

register_page(__name__, path="/resources", name="Resources")

# Define your resources data
RESOURCES = [
    # Guidelines
    dict(title="WHO – Hepatitis Overview", virus="General", category="Guidelines",
         href="https://www.who.int/health-topics/hepatitis",
         blurb="Global overview, facts, and strategy from WHO."),
    dict(title="WHO – Hepatitis B Fact Sheet (2025)", virus="HBV", category="Guidelines",
         href="https://www.who.int/news-room/fact-sheets/detail/hepatitis-b",
         blurb="Up-to-date HBV facts (transmission, prevention, care)."),
    dict(title="WHO – Hepatitis C Fact Sheet (2025)", virus="HCV", category="Guidelines",
         href="https://www.who.int/news-room/fact-sheets/detail/hepatitis-c",
         blurb="Key HCV facts, prevention, and care."),
    dict(title="CDC – Viral Hepatitis", virus="General", category="Guidelines",
         href="https://www.cdc.gov/hepatitis/index.html",
         blurb="US surveillance, basics, and prevention resources."),
    dict(title="EASL Clinical Practice Guidelines: HBV (2025)", virus="HBV", category="Guidelines",
         href="https://www.hepb.org/assets/Uploads/EASL-guidelines-May-2025.pdf",
         blurb="Latest European guidance on HBV diagnosis, treatment, monitoring."),
    dict(title="AASLD-IDSA HCV Guidance (Dec 2023, web)", virus="HCV", category="Guidelines",
         href="https://www.hcvguidelines.org/",
         blurb="Continuously updated HCV testing and treatment guidance."),

    # Databases
    dict(title="HBVdb (INSERM Lyon)", virus="HBV", category="Databases",
         href="https://hbvdb.lyon.inserm.fr/",
         blurb="Curated HBV sequences, tools (genotyper, resistance)."),
    dict(title="Los Alamos HCV Sequence & Immunology DB", virus="HCV", category="Databases",
         href="https://hcv.lanl.gov/",
         blurb="Canonical HCV sequence database and tools."),

    # Pipelines (analysis/workflows)
    dict(title="nf-core/viralrecon", virus="General", category="Pipelines",
         href="https://nf-co.re/viralrecon/2.6.0/",
         blurb="End-to-end viral WGS pipeline (Illumina/ONT), assembly & low-freq variants."),
    dict(title="V-pipe (ETH Zürich)", virus="General", category="Pipelines",
         href="https://github.com/cbg-ethz/V-pipe",
         blurb="Within-sample diversity, consensus, SNVs, haplotypes for viral NGS."),
    dict(title="Galaxy – Public Server", virus="General", category="Pipelines",
         href="https://galaxy-main.usegalaxy.org/",
         blurb="No-code workflows for viral NGS; publish & reuse pipelines."),

    # Tools (genotyping / resistance / phylo)
    dict(title="HCV-GLUE (CVR Glasgow)", virus="HCV", category="Tools",
         href="https://hcv-glue.cvr.gla.ac.uk/",
         blurb="Genotyping & DAA resistance summaries; report builder (research use)."),
    dict(title="geno2pheno[HCV]", virus="HCV", category="Tools",
         href="https://hcv.geno2pheno.org/",
         blurb="Subtype and DAA resistance interpretation (web app)."),
    dict(title="HBVdb Genotyper", virus="HBV", category="Tools",
         href="https://hbvdb.lyon.inserm.fr/HBVdb/HBVdbGenotype",
         blurb="Paste/upload sequences for HBV genotype calls."),
    dict(title="Stanford HBVseq", virus="HBV", category="Tools",
         href="https://hivdb.stanford.edu/HBV/HBVseq/development/HBVseq.html",
         blurb="Genotype and RT mutation interpretation (research use)."),
    dict(title="Genome Detective – HBV Typing Tool", virus="HBV", category="Tools",
         href="https://www.genomedetective.com/app/typingtool/hbv/",
         blurb="Batch phylogenetic HBV genotyping."),
    dict(title="Genome Detective – HCV Typing Tool", virus="HCV", category="Tools",
         href="https://www.genomedetective.com/app/typingtool/hcv/",
         blurb="Batch phylogenetic HCV genotyping."),
    dict(title="MAFFT (MSA)", virus="General", category="Tools",
         href="https://mafft.cbrc.jp/",
         blurb="Multiple sequence alignment (web & CLI)."),
    dict(title="IQ-TREE (Phylogenetics)", virus="General", category="Tools",
         href="https://iqtree.github.io/",
         blurb="Maximum-likelihood trees, ModelFinder, ultrafast bootstrap."),
    dict(title="BEAST 2 (Bayesian phylogenetics)", virus="General", category="Tools",
         href="https://www.beast2.org/",
         blurb="Time-scaled trees, phylodynamics (packages & tutorials)."),

    # Protocols (wet lab & WGS)
    dict(title="HBV WGS (HEP-TILE ONT) – protocols.io", virus="HBV", category="Protocols",
         href="https://www.protocols.io/view/hep-tile-hbv-whole-genome-sequencing-nanopore-prot-dpn45mgw.html",
         blurb="Tiled amplicon ONT protocol, pan-genotypic primer scheme."),
    dict(title="HBV ONT Whole Genome – protocols.io (PDF)", virus="HBV", category="Protocols",
         href="https://www.protocols.io/view/complete-hepatitis-b-virus-sequencing-using-an-ont-cu9uwz6w.pdf",
         blurb="Step-by-step HBV complete genome sequencing on ONT."),
    dict(title="HCV Near-whole-genome (Illumina MiSeq)", virus="HCV", category="Protocols",
         href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8473162/",
         blurb="Genotype-independent HCV assay validated on MiSeq."),
]

CATEGORIES = ["Guidelines", "Databases", "Pipelines", "Tools", "Protocols"]
VIRUSES = ["General", "HBV", "HCV"]

# Define color scheme for different categories
CATEGORY_COLORS = {
    "Guidelines": "primary",
    "Databases": "success",
    "Pipelines": "warning",
    "Tools": "info",
    "Protocols": "secondary"
}

VIRUS_COLORS = {
    "General": "dark",
    "HBV": "danger",
    "HCV": "secondary"  # Using secondary instead of custom purple for simplicity
}

def _resource_card(item):
    # Use appropriate colors for badges
    virus_badge = dbc.Badge(
        item["virus"], 
        color=VIRUS_COLORS.get(item["virus"], "secondary"), 
        className="me-1"
    )
    category_badge = dbc.Badge(
        item["category"], 
        color=CATEGORY_COLORS.get(item["category"], "info"), 
        className="me-1"
    )
    
    return dbc.Card(
        dbc.CardBody([
            html.Div([virus_badge, category_badge], className="mb-2"),
            html.H5(item["title"], className="card-title mb-2 text-primary"),
            html.P(item["blurb"], className="card-text mb-3 text-muted"),
            dbc.Button(
                "Visit Resource", 
                href=item["href"], 
                target="_blank", 
                color="primary", 
                size="sm",
                className="d-flex align-items-center"
            )
        ]),
        className="h-100 shadow-sm border-0 resource-card",
        style={"borderRadius": "12px"}
    )

layout = dbc.Container([
    # Header section
    dbc.Row([
        dbc.Col([
            html.H1("Hepatitis Resources", className="mb-2 fw-bold"),
            html.P("A curated collection of tools, guidelines, and references for hepatitis research and clinical practice.", 
                  className="text-muted mb-4")
        ], width=12)
    ]),
    
    # Filters section
    dbc.Card([
        dbc.CardBody([
            html.H5("Filter Resources", className="mb-3"),
            dbc.Row([
                dbc.Col([
                    dbc.InputGroup([
                        dbc.InputGroupText(html.I(className="fas fa-search")),
                        dcc.Input(
                            id="res-search", 
                            type="text", 
                            placeholder="Search resources…",
                            className="form-control",
                            debounce=True
                        ),
                    ]),
                ], md=6, className="mb-2"),
                dbc.Col([
                    dcc.Dropdown(
                        id="res-virus", 
                        options=[{"label": v, "value": v} for v in ["All"] + VIRUSES],
                        value="All", 
                        clearable=False,
                        placeholder="Filter by virus"
                    ),
                ], md=3, className="mb-2"),
                dbc.Col([
                    dcc.Dropdown(
                        id="res-category", 
                        options=[{"label": c, "value": c} for c in ["All"] + CATEGORIES],
                        value="All", 
                        clearable=False,
                        placeholder="Filter by category"
                    ),
                ], md=3, className="mb-2")
            ], className="g-2 align-items-end"),
        ])
    ], className="shadow-sm mb-4 border-0", style={"borderRadius": "12px"}),
    
    # Results counter
    html.Div(id="res-count", className="text-muted mb-3 fw-medium"),
    
    # Resources grid
    dbc.Row(id="res-grid", className="g-4"),
    
    # Loading component for better UX
    dcc.Loading(
        id="loading-resources",
        type="circle",
        children=html.Div(id="loading-output")
    )
], fluid=True, className="py-4 px-4 px-md-5")

# Add custom CSS through an assets folder or inline style
# Create a file in assets/style.css with the following content:
"""
.resource-card {
    transition: transform 0.2s, box-shadow 0.2s;
}

.resource-card:hover {
    transform: translateY(-5px);
    box-shadow: 0 0.5rem 1rem rgba(0, 0, 0, 0.15) !important;
}
"""

@callback(
    Output("res-grid", "children"),
    Output("res-count", "children"),
    Input("res-search", "value"),
    Input("res-virus", "value"),
    Input("res-category", "value"),
)
def _filter_resources(q, virus, category):
    q = (q or "").strip().lower()
    def keep(x):
        virus_ok = (virus == "All") or (x["virus"] == virus)
        cat_ok = (category == "All") or (x["category"] == category)
        text = " ".join([x["title"], x["blurb"], x["virus"], x["category"]]).lower()
        return virus_ok and cat_ok and (q in text if q else True)

    items = [x for x in RESOURCES if keep(x)]
    
    if not items:
        return [
            dbc.Col([
                html.Div(
                    html.P("No resources match your filters. Try adjusting your search criteria.", 
                          className="text-center text-muted p-5"),
                    className="d-flex justify-content-center align-items-center"
                )
            ], width=12)
        ], f"No resources found"
    
    cards = [dbc.Col(_resource_card(x), md=6, lg=4, className="mb-3") for x in items]
    
    count_text = f"Showing {len(items)} resource{'s' if len(items) != 1 else ''}"
    if virus != "All" or category != "All" or q:
        count_text += " (filtered)"
    
    return cards, count_text