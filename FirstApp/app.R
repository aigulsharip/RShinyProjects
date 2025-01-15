#Simple Shiny App with inputs and output
#Install packages


# Load R packages
library(shiny)
library(shinythemes)

# Define UI
ui <- fluidPage(theme = shinytheme("cosmo"),
  navbarPage(
    # theme = "cerulean",  # <--- To use a theme, uncomment this
    "My first app",
    tabPanel("Navbar 1",
             sidebarPanel(
               #h1("Sidebar"),
               tags$h3("Input:"),
               textInput("txt1", "First Name:", ""),
               textInput("txt2", "Last Name:", ""),
               
             ), # sidebarPanel
             mainPanel(
                          #h1("MainPanel"),
                          h3("Full Name"),
                          verbatimTextOutput("txtout"),

             ) # mainPanel
             
    ), # Navbar 1, tabPanel
    tabPanel("Navbar 2", "This panel is intentionally left blank"),
    tabPanel("Navbar 3", "This panel is intentionally left blank")

  ) # navbarPage
) # fluidPage


# Define server function  
server <- function(input, output) {
  
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = " " )
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)
