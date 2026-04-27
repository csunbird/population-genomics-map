# =============================================================================
# app.R — Entry point
# Run with: shiny::runApp() from the project root
#           or deploy via rsconnect::deployApp()
# =============================================================================

source("global.R")
source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)
