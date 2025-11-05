import os
import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from flask import Flask, redirect
from utils.data_loader import load_and_preprocess_data

# Create Flask app first
server = Flask(__name__)
server.secret_key = os.environ.get('SECRET_KEY', 'default-secret-key')

# Add redirect before creating Dash app - ONLY ONE REDIRECT
@server.route('/')
def redirect_to_dashboard():
    return redirect('/dashboard', code=302)

app = dash.Dash(
    __name__,
    server=server,
    use_pages=True,
    external_stylesheets=[
        dbc.themes.BOOTSTRAP,
        "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css"
    ],
    suppress_callback_exceptions=True,
    pages_folder="pages"  # Explicitly specify pages folder
)

def navbar():
    # Get page registry and create navigation
    reg = {p["name"]: p["path"] for p in dash.page_registry.values()}
    order = ["Dashboard", "About", "Resources", "Contact"]
    items = [dbc.NavItem(dbc.NavLink(name, href=reg[name], active="exact"))
             for name in order if name in reg]
    return dbc.Navbar(
        dbc.Container([
            dbc.NavbarBrand("Hepatitis Dashboard", className="fw-bold"),
            dbc.Nav(items, pills=True, navbar=True)
        ]),
        color="primary", dark=True, sticky="top", className="mb-4",
    )

app.layout = dbc.Container(
    [
        navbar(),
        dash.page_container,
    ],
    fluid=True,
)

# Load data and store in server config
app.server.config["DATA_STORE"] = load_and_preprocess_data()

if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=8050)
