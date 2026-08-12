#!/usr/bin/env python3
"""
SNPster - Genomic Analysis Platform
Main application entry point for Posit Cloud deployment
"""

import sys
import os
from pathlib import Path

# Add the web_app directory to Python path for imports
current_dir = Path(__file__).parent
web_app_dir = current_dir / "web_app"
sys.path.insert(0, str(web_app_dir))

# Import the Shiny app
try:
    from web_app.shiny_ui import app
except ImportError as e:
    print(f"Error importing the Shiny app: {e}")
    # Create a minimal error app
    from shiny import App, ui
    
    error_ui = ui.page_fluid(
        ui.h1("Import Error"),
        ui.p(f"Failed to import app: {str(e)}"),
        ui.p("Check that all dependencies are installed.")
    )
    
    def error_server(input, output, session):
        pass
    
    app = App(error_ui, error_server)

# This is the main app object that Posit Cloud will use
if __name__ == "__main__":
    app.run()
