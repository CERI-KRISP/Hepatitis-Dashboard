from dash import html, dcc, register_page, callback, Output, Input
import dash_bootstrap_components as dbc
from dash import ctx
from dash import callback_context as ctx
from datetime import datetime
from weasyprint import HTML
from dash.exceptions import PreventUpdate
import pandas as pd

register_page(__name__, path="/about", name="About")

# -----------------------
# Reusable components with enhanced styling
# -----------------------

def _fact_card(value: str, caption: str, icon: str = None) -> dbc.Card:
    return dbc.Card(
        dbc.CardBody(
            [
                html.Div(
                    [
                        html.I(className=icon, style={"fontSize": "2rem", "color": "#3a8ca7"}) if icon else None,
                        html.H3(
                            value,
                            className="text-center mb-2",
                            style={"fontWeight": 800, "color": "#3a8ca7", "fontSize": "1.8rem"}
                        ),
                    ],
                    className="text-center"
                ),
                html.P(
                    caption,
                    className="text-center mb-0 text-muted",
                    style={"fontSize": "0.95rem", "lineHeight": "1.4", "fontWeight": 400}
                ),
            ]
        ),
        className="border-0 shadow-sm mb-3 bg-white",
        style={"transition": "transform 0.2s", "borderRadius": "12px"}
    )

def key_facts_card():
    """Enhanced 'Key facts at a glance' for sidebar."""
    return html.Div(
        [
            html.H4("📊 Key Facts at a Glance", 
                   className="text-center mb-4", 
                   style={"color": "#2c3e50", "fontWeight": 600, "paddingBottom": "10px", "borderBottom": "2px solid #4ca1af"}),
            _fact_card("14,529", "Whole genomes analyzed", "bi bi-dna"),
            _fact_card("141", "Countries covered", "bi bi-globe"),
            _fact_card("18+", "Genotypes tracked", "bi bi-diagram-3"),
            _fact_card("Real-time", "Analytics platform", "bi bi-lightning"),
        ],
        className="h-100 p-4",
        style={"backgroundColor": "#f8fafc", "borderRadius": "15px", "border": "1px solid #e9ecef"}
    )

def _section_header(title: str, emoji: str) -> html.H2:
    return html.H2(
        f"{emoji} {title}",
        style={
            "color": "#2c3e50",
            "textAlign": "center",
            "marginBottom": "40px",
            "fontWeight": 700,
            "paddingBottom": "15px",
            "borderBottom": "3px solid #4ca1af",
            "position": "relative"
        }
    )

def _feature_card(title: str, items: list, icon: str, color: str) -> dbc.Card:
    return dbc.Card(
        dbc.CardBody([
            html.Div(
                html.I(className=icon, style={"fontSize": "2.5rem", "color": color}),
                className="text-center mb-3"
            ),
            html.H4(title, className="text-center mb-3", style={"color": "#2c3e50", "fontWeight": 600}),
            html.Ul([html.Li(item, style={"marginBottom": "8px"}) for item in items], className="ps-3")
        ]),
        className="h-100 shadow-sm border-0",
        style={"borderRadius": "12px", "transition": "all 0.3s ease"}
    )

# -----------------------
# New Enhanced Sections
# -----------------------

def dashboard_mission_section():
    return dbc.Container(
        [
            html.H1("About This Dashboard", style={"color": "#2c3e50", "marginBottom": "30px"}),
            
            dbc.Row([
                dbc.Col([
                    html.H2("🎯 Our Mission", style={"color": "#2c3e50", "marginBottom": "20px"}),
                    html.P(
                        "This platform transforms global hepatitis surveillance by integrating genomic sequencing data "
                        "with epidemiological insights to identify critical gaps and guide public health interventions.",
                        className="lead mb-4"
                    ),
                    html.P(
                        "We bridge the divide between vast genomic datasets and actionable public health intelligence, "
                        "enabling researchers and policymakers to optimize resource allocation and track progress toward "
                        "WHO elimination goals."
                    ),
                    html.Div([
                        html.H5("Core Objectives", style={"color": "#2c3e50", "marginTop": "25px", "marginBottom": "15px"}),
                        html.Ul([
                            html.Li("Identify geographic and genotypic surveillance gaps"),
                            html.Li("Monitor emerging drug resistance patterns"),
                            html.Li("Facilitate data-driven public health interventions"),
                            html.Li("Support global hepatitis elimination goals"),
                        ])
                    ])
                ], md=8),
                dbc.Col(key_facts_card(), md=4)
            ])
        ],
        fluid=True,
        className="my-5"
    )

def user_types_section():
    return dbc.Container(
        [
            _section_header("👥 Designed For Multiple Stakeholders", "👥"),
            
            dbc.Row([
                dbc.Col([
                    _feature_card(
                        "Researchers & Academics",
                        [
                            "Identify sequencing gaps and priorities",
                            "Track viral evolution and diversity",
                            "Access curated, analysis-ready datasets",
                            "Generate publication-ready visualizations"
                        ],
                        "bi bi-mortarboard",
                        "#3182bd"
                    )
                ], md=4, className="mb-3"),
                dbc.Col([
                    _feature_card(
                        "Public Health Agencies",
                        [
                            "Monitor elimination progress metrics",
                            "Target interventions to high-need areas",
                            "Allocate sequencing resources efficiently",
                            "Generate automated surveillance reports"
                        ],
                        "bi bi-building",
                        "#d94801"
                    )
                ], md=4, className="mb-3"),
                dbc.Col([
                    _feature_card(
                        "Global Health Partners",
                        [
                            "Coordinate international surveillance",
                            "Identify priority regions for support",
                            "Measure impact of interventions",
                            "Inform policy and funding decisions"
                        ],
                        "bi bi-globe",
                        "#2c3e50"
                    )
                ], md=4, className="mb-3"),
            ])
        ],
        fluid=True,
        className="my-5"
    )

def technical_specs_section():
    return dbc.Container(
        [
            _section_header("🛠️ Technical Implementation", "🛠️"),
            
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader(html.H5("Architecture & Stack", className="mb-0")),
                        dbc.CardBody([
                            html.Ul([
                                html.Li("Frontend: Plotly Dash with Bootstrap components"),
                                html.Li("Backend: Python with Pandas, NumPy, Scikit-learn"),
                                html.Li("Visualization: Plotly, Mapbox for interactive mapping"),
                                html.Li("Deployment: Docker containerization, cloud-ready"),
                                html.Li("Data Processing: Automated ETL pipelines")
                            ])
                        ])
                    ], className="h-100 shadow-sm")
                ], md=6, className="mb-3"),
                
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader(html.H5("Data Pipeline", className="mb-0")),
                        dbc.CardBody([
                            html.Ul([
                                html.Li("Genomic data: GenBank API integration"),
                                html.Li("Epidemiological data: IHME GBD, WHO reports"),
                                html.Li("Geographic data: Standardized ISO-3 mapping"),
                                html.Li("Quality control: Automated validation checks"),
                                html.Li("Update frequency: Monthly data refreshes")
                            ])
                        ])
                    ], className="h-100 shadow-sm")
                ], md=6, className="mb-3")
            ])
        ],
        fluid=True,
        className="my-5"
    )

def hepatitis_context_card():
    return dbc.Card([
        dbc.CardHeader(
            html.H5("🧬 Hepatitis Context & Significance", className="mb-0")
        ),
        dbc.CardBody([
            html.P(
                "This dashboard focuses on Hepatitis B and C, which collectively affect 325 million people worldwide "
                "and cause over 1.1 million annual deaths. These viruses present significant global health challenges "
                "due to their potential for chronic infection, progression to cirrhosis and hepatocellular carcinoma, "
                "and uneven global distribution.",
                className="mb-3"
            ),
            html.P(
                "Genomic surveillance is critical for understanding transmission patterns, monitoring drug resistance, "
                "and guiding elimination efforts. This platform addresses the urgent need for integrated, accessible "
                "genomic epidemiology tools.",
                className="mb-3"
            ),
            dbc.Button(
                "Learn More About Hepatitis",
                href="/resources#hepatitis-background", 
                color="outline-primary",
                size="sm"
            )
        ])
    ], className="mb-4")

def data_methodology_section():
    return dbc.Container(
        [
            _section_header("📊 Data & Methodology", "📊"),
            
            dbc.Row([
                dbc.Col(
                    _feature_card(
                        "Data Sources",
                        [
                            "10,996 HBV + 3,533 HCV whole genome sequences",
                            "IHME Global Burden of Disease estimates",
                            "World Bank & UN population data",
                            "Standardized geographic coordinates"
                        ],
                        "bi bi-database",
                        "#3182bd"
                    ),
                    md=6, className="mb-4"
                ),
                dbc.Col(
                    _feature_card(
                        "Analytical Methods",
                        [
                            "Genomic data integration & standardization",
                            "Burden-adjusted coverage metrics",
                            "Geographic disparity analysis",
                            "Drug resistance mutation tracking"
                        ],
                        "bi bi-gear",
                        "#d94801"
                    ),
                    md=6, className="mb-4"
                )
            ]),
        ],
        fluid=True,
        className="my-5"
    )

def surveillance_gaps_section():
    return dbc.Container(
        [
            _section_header("🔍 Global Surveillance Insights", "🔍"),
            
            dbc.Row([
                dbc.Col(
                    dbc.Card([
                        dbc.CardHeader(html.H5("🌍 Geographic Disparities", className="mb-0", style={"color": "#2c3e50"})),
                        dbc.CardBody([
                            html.P("Sequencing concentrated in China/US while high-endemicity regions like Africa are severely underrepresented.",
                                  className="text-muted mb-2"),
                            html.Div([
                                html.Strong("Key finding: ", style={"color": "#dc3545"}),
                                "Tens of thousands of additional sequences needed from underrepresented regions"
                            ], style={"backgroundColor": "#fff3cd", "padding": "10px", "borderRadius": "5px"})
                        ])
                    ], className="h-100 shadow-sm border-0"),
                    width=4, className="mb-3"
                ),
                dbc.Col(
                    dbc.Card([
                        dbc.CardHeader(html.H5("🧬 Genotypic Gaps", className="mb-0", style={"color": "#2c3e50"})),
                        dbc.CardBody([
                            html.P("Critical genotypes (HBV-E, HCV-5/8) are significantly under-sampled, limiting understanding of viral diversity.",
                                  className="text-muted mb-2"),
                            html.Div([
                                html.Strong("Impact: ", style={"color": "#fd7e14"}),
                                "Hinders personalized treatment approaches and vaccine development"
                            ], style={"backgroundColor": "#e7f4e4", "padding": "10px", "borderRadius": "5px"})
                        ])
                    ], className="h-100 shadow-sm border-0"),
                    width=4, className="mb-3"
                ),
                dbc.Col(
                    dbc.Card([
                        dbc.CardHeader(html.H5("⚖️ Burden Mismatch", className="mb-0", style={"color": "#2c3e50"})),
                        dbc.CardBody([
                            html.P("Critical misalignment between disease burden and sequencing efforts across regions.",
                                  className="text-muted mb-2"),
                            html.Div([
                                html.Strong("Recommendation: ", style={"color": "#0d6efd"}),
                                "Targeted sequencing in high-burden, low-sequence regions"
                            ], style={"backgroundColor": "#cce7ff", "padding": "10px", "borderRadius": "5px"})
                        ])
                    ], className="h-100 shadow-sm border-0"),
                    width=4, className="mb-3"
                )
            ], className="mb-4"),
        ],
        fluid=True,
        className="my-5"
    )

def interactive_dashboard_section():
    return dbc.Container(
        [
            html.H3("📈 Interactive Dashboard Preview", 
                   style={"textAlign": "center", "color": "#2c3e50", "marginBottom": "30px", "fontWeight": 700}),
            
            # Key statistics
            dbc.Row([
                dbc.Col(_fact_card("10,996", "HBV Genomes", "bi bi-dna"), width=3, className="mb-3"),
                dbc.Col(_fact_card("3,533", "HCV Genomes", "bi bi-dna"), width=3, className="mb-3"),
                dbc.Col(_fact_card("89 + 52", "Countries", "bi bi-globe"), width=3, className="mb-3"),
                dbc.Col(_fact_card("10 + 8", "Genotypes", "bi bi-diagram-3"), width=3, className="mb-3"),
            ], className="mb-4"),
            
            # Dual maps
            html.Div([
                html.H4("Global Sequence Distribution", 
                       style={"textAlign": "center", "color": "#2c3e50", "marginBottom": "25px", "fontWeight": 600}),
                
                dbc.Row([
                    dbc.Col(
                        html.Div([
                            html.H5("Hepatitis B (HBV)", 
                                   style={"textAlign": "center", "color": "#3182bd", "marginBottom": "15px", "fontWeight": 600}),
                            html.Img(
                                src="/assets/HBV_map.png",
                                style={
                                    "width": "100%", 
                                    "height": "auto",
                                    "borderRadius": "12px",
                                    "boxShadow": "0 6px 12px rgba(0,0,0,0.1)",
                                    "border": "3px solid #3182bd"
                                },
                                alt="HBV global sequence distribution"
                            ),
                            html.P("10,996 sequences across 89 countries",
                                  className="text-center text-muted mt-3",
                                  style={"fontSize": "0.9rem"})
                        ]),
                        width=6, className="mb-4"
                    ),
                    dbc.Col(
                        html.Div([
                            html.H5("Hepatitis C (HCV)", 
                                   style={"textAlign": "center", "color": "#d94801", "marginBottom": "15px", "fontWeight": 600}),
                            html.Img(
                                src="/assets/HCV_map.png",
                                style={
                                    "width": "100%", 
                                    "height": "auto",
                                    "borderRadius": "12px",
                                    "boxShadow": "0 6px 12px rgba(0,0,0,0.1)",
                                    "border": "3px solid #d94801"
                                },
                                alt="HCV global sequence distribution"
                            ),
                            html.P("3,533 sequences across 52 countries",
                                  className="text-center text-muted mt-3",
                                  style={"fontSize": "0.9rem"})
                        ]),
                        width=6, className="mb-4"
                    ),
                ]),
            ], style={"marginBottom": "40px"}),
            
            # Feature highlights
            dbc.Row([
                dbc.Col(
                    _feature_card(
                        "Interactive Visualizations",
                        ["Real-time surveillance maps", "Genotype distribution", "Temporal trends", "Coverage metrics"],
                        "bi bi-bar-chart",
                        "#4ca1af"
                    ),
                    width=6, className="mb-4"
                ),
                dbc.Col(
                    _feature_card(
                        "Advanced Analytics",
                        ["Drug resistance tracking", "Geographic analysis", "Regional comparisons", "Custom filtering"],
                        "bi bi-graph-up",
                        "#2c3e50"
                    ),
                    width=6, className="mb-4"
                ),
            ]),
            
            # CTA buttons
            html.Div(
                id="report-section",
                children=[
                    html.H4("Download Genomic Surveillance Report", style={
                        "textAlign": "center",
                        "color": "#2c3e50",
                        "marginBottom": "20px",
                        "fontWeight": 700
                    }),
                
                    html.Div([
                        dbc.Label("Select Report Type:", html_for="report-type-dropdown",
                                  style={"fontWeight": "bold", "color": "#2c3e50", "marginBottom": "10px"}),
                
                        dcc.Dropdown(
                            id="report-type-dropdown",
                            options=[
                                {"label": "📘 Full Comprehensive Report", "value": "full"},
                                {"label": "📄 Executive Summary", "value": "summary"},
                                {"label": "🔍 Gap Analysis Only", "value": "gaps"},
                                {"label": "⚙️ Technical Methodology", "value": "methods"},
                            ],
                            value="full",
                            style={"width": "100%", "maxWidth": "350px", "margin": "0 auto"},
                            clearable=False,
                        ),
                
                        dbc.Button(
                            "Download Selected Report",
                            id="btn-download-report",
                            color="secondary",
                            size="lg",
                            className="mt-3",
                            style={
                                "backgroundColor": "#2c3e50",
                                "border": "none",
                                "borderRadius": "8px",
                                "fontWeight": 600,
                                "padding": "12px 24px"
                            }
                        ),
                
                        dcc.Download(id="download-data-report"),
                        
                        dbc.Button(
                            "Explore Full Dashboard",
                            color="primary",
                            href="/dashboard",
                            size="lg",
                            className="mt-3 mx-2 px-4",
                            style={
                                "backgroundColor": "#4ca1af",
                                "border": "none",
                                "borderRadius": "8px",
                                "fontWeight": 600,
                                "padding": "12px 24px"
                            }
                        ),
                
                    ], className="text-center")
                ],
                className="text-center mt-5 py-4",
                style={
                    "backgroundColor": "#f8fafc",
                    "borderRadius": "15px",
                    "padding": "30px",
                    "maxWidth": "600px",
                    "margin": "0 auto"
                }
            ),
        ],
        fluid=True,
        className="my-5 py-4",
        style={"backgroundColor": "#f8fafc", "borderRadius": "15px"}
    )

def acknowledgments_section():
    return dbc.Container(
        [
            dbc.Card(
                dbc.CardBody(
                    [
                        html.H2("Acknowledgments & Resources", style={"color": "#2c3e50", "marginBottom": "16px"}),
                        html.P(
                            "Developed by [Your Lab / Institution] in collaboration with [Collaborators].",
                            style={"marginBottom": "6px"},
                        ),
                        html.P(
                            html.Span(
                                [
                                    html.Strong("Data Sources: "),
                                    "GenBank, IHME Global Burden of Disease, World Bank, UN Population Division.",
                                ]
                            )
                        ),
                        html.P(
                            html.Span(
                                [
                                    html.Strong("Technologies: "),
                                    "Plotly Dash, Python, Pandas, NumPy, Mapbox, Docker.",
                                ]
                            )
                        ),
                        html.Hr(),
                        html.P(
                            html.Span(
                                [
                                    html.Strong("Suggested Citation: "),
                                    "[Author(s)], Hepatitis Genomic Surveillance Dashboard, Version [X.X], [Year]. DOI: [if applicable].",
                                ]
                            ),
                            className="mb-2",
                        ),
                        html.P(
                            html.Span(
                                [
                                    html.Strong("Contact: "),
                                    "[email@institution.edu] for technical support or collaboration inquiries.",
                                ]
                            ),
                            className="mb-3",
                        ),
                        dbc.Button(
                            "Explore Full Dashboard",
                            href="/dashboard",
                            color="primary",
                            className="me-2",
                            style={"backgroundColor": "#4ca1af", "border": "none"},
                        ),
                        dbc.Button(
                            "Access Research Library", 
                            href="/resources",
                            color="secondary",
                            style={"backgroundColor": "#2c3e50", "border": "none"},
                        ),
                    ]
                ),
                className="shadow-sm border-0",
                style={"borderRadius": "12px", "backgroundColor": "#f8f9fa"},
            )
        ],
        fluid=True,
        className="my-5",
    )

# Keep your existing callback and PDF generation functions (they remain the same)
@callback(
    Output("download-data-report", "data"),
    Input("btn-download-report", "n_clicks"),
    Input("report-type-dropdown", "value"),
    prevent_initial_call=True
)
def download_selected_report(n_clicks, report_type):
    if n_clicks:
        pdf_content = generate_data_report(report_type)
        filename = f"Hepatitis_Report_{report_type}_{datetime.now().strftime('%Y%m%d')}.pdf"
        return dcc.send_bytes(pdf_content, filename)
    raise PreventUpdate
    
def generate_data_report(report_type="full") -> bytes:
    # ... keep your existing PDF generation code unchanged ...
    pass

# -----------------------
# Enhanced Page Layout
# -----------------------

layout = html.Div(
    [
        # Enhanced Banner
        dbc.Container(
            fluid=True,
            children=[
                html.Div(
                    className="hepatitis-banner",
                    style={
                        "background": "linear-gradient(rgba(0, 0, 0, 0.7), rgba(0, 0, 0, 0.7)), url('/assets/Hepatitis_banner.png') center/cover no-repeat",
                        "color": "white",
                        "padding": "5rem 1rem",
                        "textAlign": "center",
                        "marginBottom": "3rem",
                        "position": "relative",
                    },
                    children=[
                        html.H1("HEPATITIS GENOMIC SURVEILLANCE DASHBOARD", style={
                            "fontSize": "3.5rem", "fontWeight": "bold",
                            "textShadow": "2px 2px 8px rgba(0,0,0,0.7)",
                            "marginBottom": "1rem"
                        }),
                        html.Hr(style={
                            "width": "120px", "height": "4px", "backgroundColor": "#4ca1af",
                            "margin": "0 auto", "border": "none", "marginBottom": "1.5rem"
                        }),
                        html.P("Open-access analytics integrating global HBV/HCV genomic and epidemiological data to inform elimination strategies.",
                            style={"fontSize": "1.5rem", "fontWeight": 300, "opacity": 0.9}
                        ),
                    ],
                )
            ],
        ),

        # Main content container
        dbc.Container(
            [
                # New dashboard-focused sections
                dashboard_mission_section(),
                user_types_section(),
                hepatitis_context_card(),  # Streamlined hepatitis context
                technical_specs_section(),
                data_methodology_section(),
                surveillance_gaps_section(), 
                interactive_dashboard_section(),
                acknowledgments_section(),
            ],
            fluid=True,
            style={'maxWidth': '1440px'}
        ),
    ],
    style={"backgroundColor": "#ffffff", "minHeight": "100vh"}
)
