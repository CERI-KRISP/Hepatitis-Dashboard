from dash import html, dcc, register_page, Input, Output, State, callback
import dash_bootstrap_components as dbc
from urllib.parse import quote_plus

register_page(__name__, path="/contact", name="Contact")

# Multiple organizations data - replace with your actual details
ORGANIZATIONS = [
    {
        "name": "Centre for Epidemic Response and Innovation (CERI)",
        "logo": "/assets/ceri_logo.png",  # Replace with your actual logo path
        "description": "We build in our successful experience of the Network for Genomic Surveillance in South Africa (NGS-SA) with like minded people who do generate ideas and realize them.",
        "address": "Van Der Byl Rd, Stellenbosch Central, Stellenbosch, 7600, South Africa",
        "phone": "+27 (0) 21 808 3815",
        "email": "ceri@sun.ac.za",
        "hours": "Monday-Friday: 8:00 AM - 4:30 PM"
    },
    {
        "name": "Francis Crick Institute (CRICK)",
        "logo": "/assets/CRICK.png",  # Replace with your actual logo path
        "description": "The Francis Crick Institute is an independent charity, established to be a UK flagship for discovery research in biomedicine.",
        "address": "1 Midland Rd, London NW1 1AT, United Kingdom",
        "phone": "+44 20 3796 0000",
        "email": "contact@viralgenomics.org",
        "hours": "Monday-Friday: 10:00 AM - 4:00 PM, Wednesday: 10:00 AM - 8:00 PM"
    },
    {
        "name": "Africa Health Research Institute (AHRI)",
        "logo": "/assets/AHRI_logo.png",  # Replace with your actual logo path
        "description": "AHRI’s research combines population, basic and translational, social, and clinical sciences to understand and intervene in the health and well-being of South African communities.",
        "address": "Nelson R. Mandela School of Medicine, 719 Umbilo Road, Durban, 4013, South Frica",
        "phone": "+27 (0) 31 521 0038",
        "email": "info@ahri.org",
        "hours": ""
    }
]

# Team members data - replace with your actual team
TEAM_MEMBERS = [
    {
        "name": "Prof. Tulio De Oliveira",
        "title": "Director of Research",
        "organization": "Centre for Epidemic Response and Innovation (CERI)",
        "email": "tulio@sun.ac.za",
        "phone": "+27 (0) 31 260 4898",
        "expertise": "Professor of Bioinformatics at the School for Data Science and Computational Thinking, Stellenbosch University and at the College of Health Sciences at University of KwaZulu-Natal.",
        "image": "/assets/Tulio.jpeg"  # Replace with actual image paths
    },
    {
        "name": "Prof. Michael Chen",
        "title": "Lead Bioinformatician",
        "organization": "Viral Genomics Center",
        "email": "mchen@viralgenomics.org",
        "phone": "+1 (555) 234-5002",
        "expertise": "Viral Sequencing & Analysis",
        "image": "/assets/team2.jpg"
    },
    {
        "name": "Prof. Willem Hanekom",
        "title": "Executive Director",
        "organization": "Africa Health Research Institute (AHRI)",
        "email": "info@ahri.org",
        "phone": "+27 (0)31 521 0038",
        "expertise": "Clinician Scientist",
        "image": "/assets/Willem_Hanekom_website.png"
    },
]

# Department contacts by organization
DEPARTMENTS = {
    "Centre for Epidemic Research and Innovation (CERI)": [
        {
            "name": "Research Collaboration",
            "email": "cbaxter@sun.ac.za",
            "phone": "27 (0) 83 515 6929"
        },
        {
            "name": "Sequencing Services",
            "email": "lavanya@sun.ac.za",
            "phone": "+27 (0) 31 240 1887"
        },
        {
            "name": "Bioinformatics Support",
            "email": "ewilkinson@sun.ac.za",
            "phone": "+27 (0) 21 808 3472"
        },
    ],
    "Viral Genomics Center": [
        {
            "name": "General Inquiries",
            "email": "contact@viralgenomics.org",
            "phone": "+1 (555) 234-5678"
        },
        {
            "name": "Sequencing Services",
            "email": "sequencing@viralgenomics.org",
            "phone": "+1 (555) 234-6000"
        },
        {
            "name": "Bioinformatics Support",
            "email": "bioinfo@viralgenomics.org",
            "phone": "+1 (555) 234-7000"
        },
    ],
    "Global Health Partnership": [
        {
            "name": "General Inquiries",
            "email": "info@globalheppartnership.org",
            "phone": "+1 (555) 345-6789"
        },
        {
            "name": "Partnership Opportunities",
            "email": "partnerships@globalheppartnership.org",
            "phone": "+1 (555) 345-7000"
        },
        {
            "name": "Training Programs",
            "email": "training@globalheppartnership.org",
            "phone": "+1 (555) 345-8000"
        },
    ]
}

def maps_embed_src(org, zoom=15):
    # Prefer lat/lng if present; else fall back to address
    if "lat" in org and "lng" in org and org["lat"] and org["lng"]:
        q = f"{org['lat']},{org['lng']}"
    else:
        q = quote_plus(org.get("address", org.get("name", "")))

    # Works without an API key
    return f"https://www.google.com/maps?q={q}&z={zoom}&output=embed"

def maps_open_link(org):
    # For the “Open in Google Maps” button (optional)
    if "lat" in org and "lng" in org and org["lat"] and org["lng"]:
        q = f"{org['lat']},{org['lng']}"
        return f"https://www.google.com/maps?q={q}"
    q = quote_plus(org.get("address", org.get("name", "")))
    return f"https://www.google.com/maps/search/?api=1&query={q}"


def _organization_card(org):
    return dbc.Card([
        dbc.CardBody([
            html.Div([
                html.Img(
                    src=org["logo"], 
                    className="mb-3",
                    style={"height": "60px"}
                ),
                html.H4(org["name"], className="card-title"),
                html.P(org["description"], className="card-text text-muted"),
            ], className="text-center mb-4"),
            
            html.Div([
                html.Div([
                    html.I(className="fas fa-map-marker-alt text-primary me-2"),
                    html.Span(org["address"], className="text-muted")
                ], className="mb-2"),
                html.Div([
                    html.I(className="fas fa-phone text-primary me-2"),
                    html.Span(org["phone"], className="text-muted")
                ], className="mb-2"),
                html.Div([
                    html.I(className="fas fa-envelope text-primary me-2"),
                    html.A(org["email"], href=f"mailto:{org['email']}", 
                          className="text-muted text-decoration-none")
                ], className="mb-2"),
                html.Div([
                    html.I(className="fas fa-clock text-primary me-2"),
                    html.Span(org["hours"], className="text-muted")
                ])
            ])
        ])
    ], className="h-100")

def _team_member_card(member):
    return dbc.Card(
        [
            dbc.CardBody([
                # avatar
                html.Div(
                    html.Img(
                        src=member["image"],
                        alt=f"{member['name']} portrait"
                    ),
                    className="avatar-ring mx-auto mb-3"  # <- rings + circle crop
                ),

                # text
                html.H5(member["name"], className="card-title mb-1"),
                html.H6(member["title"], className="card-subtitle text-muted mb-2"),
                html.P(member["expertise"], className="card-text"),

                # actions
                html.Div([
                    html.A(
                        html.I(className="fas fa-envelope me-2"),
                        href=f"mailto:{member['email']}",
                        className="text-decoration-none text-dark me-3"
                    ),
                    html.A(
                        html.I(className="fas fa-phone me-2"),
                        href=f"tel:{member['phone']}",
                        className="text-decoration-none text-dark"
                    ),
                ], className="mt-3"),
                dbc.Button("Contact", href=f"mailto:{member['email']}",
                           color="primary", size="sm", className="mt-2"),
            ])
        ],
        className="h-100 text-center team-card"
    )


def _department_contact(dept, org_name):
    return dbc.Card(
        dbc.CardBody([
            html.H6(org_name, className="card-subtitle text-primary mb-2"),
            html.H5(dept["name"], className="card-title"),
            html.Div([
                html.Div([
                    html.I(className="fas fa-envelope me-2 text-primary"),
                    html.A(dept["email"], href=f"mailto:{dept['email']}", className="text-decoration-none")
                ], className="mb-2"),
                html.Div([
                    html.I(className="fas fa-phone me-2 text-primary"),
                    html.A(dept["phone"], href=f"tel:{dept['phone']}", className="text-decoration-none")
                ])
            ])
        ]),
        className="h-100 department-card"
    )

layout = dbc.Container([
    # Header
    dbc.Row([
        dbc.Col([
            html.H1("Contact Our Organizations", className="mb-3 fw-bold text-center"),
            html.P("Connect with our collaborative network of hepatitis research institutions", 
                  className="lead text-muted text-center mb-5"),
        ], width=12)
    ]),
    
    # Organization cards
    dbc.Row([
        html.H2("Our Organizations", className="mb-4 text-center")
    ]),
    dbc.Row([
        dbc.Col(_organization_card(org), md=4, className="mb-4") for org in ORGANIZATIONS
    ], className="mb-5"),
    
    # Department contacts
    dbc.Row([
        dbc.Col([
            html.H2("Department Contacts", className="mb-4 text-center")
        ], width=12)
    ]),
    dbc.Row([
        dbc.Col(_department_contact(dept, org_name), md=4, className="mb-3") 
        for org_name, departments in DEPARTMENTS.items() 
        for dept in departments
    ], className="mb-5"),
    
    # Team section
    dbc.Row([
        dbc.Col([
            html.H2("Our Team", className="mb-4 text-center")
        ], width=12)
    ]),
    dbc.Row([
        dbc.Col(_team_member_card(member), md=4, className="mb-4") for member in TEAM_MEMBERS
    ], className="mb-5"),
    
    # Map section for each organization (optional)
    dbc.Row([
        dbc.Col([
            html.H2("Find Our Locations", className="mb-4 text-center")
        ], width=12)
    ]),
    dbc.Row([
		dbc.Col([
			dbc.Card([
				dbc.CardHeader(org["name"]),
				dbc.CardBody([
					html.Iframe(
						src=maps_embed_src(org),           # 👈 dynamic per org
						width="100%",
						height="300",
						style={"border": "0"},
						allow="fullscreen",
						referrerPolicy="no-referrer-when-downgrade",
						sandbox="allow-same-origin allow-scripts allow-popups allow-forms"
					),
					dbc.Button("Open in Google Maps",
							   href=maps_open_link(org),
							   target="_blank",
							   color="link",
							   className="mt-2 p-0")
				])
			])
    ], md=4, className="mb-4")
    for org in ORGANIZATIONS
])
], fluid=True, className="py-5 px-4 px-md-5")