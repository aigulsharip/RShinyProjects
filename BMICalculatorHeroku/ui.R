############################################
# A BMI Calculator is a tool that computes the Body Mass Index (BMI) of an 
# individual. BMI is a measure derived from a person’s weight and height and 
# is used as a quick and simple way to assess whether a person is underweight, 
# normal weight, overweight, or obese.

############################################
library(shiny)
library(shinythemes)
library(markdown)

####################################
# User Interface                   #
####################################
fluidPage(theme = shinytheme("spacelab"),
              navbarPage("BMI Calculator:",
                         
                         tabPanel("Home",
                                  # Input values
                                  sidebarPanel(
                                    HTML("<h3>Input parameters</h3>"),
                                    sliderInput("height", 
                                                label = "Height", 
                                                value = 175, 
                                                min = 40, 
                                                max = 250),
                                    sliderInput("weight", 
                                                label = "Weight", 
                                                value = 70, 
                                                min = 20, 
                                                max = 100),
                                    
                                    actionButton("submitbutton", 
                                                 "Submit", 
                                                 class = "btn btn-primary")
                                  ),
                                  
                                  mainPanel(
                                    tags$label(h3('Status/Output')), # Status/Output Text Box
                                    verbatimTextOutput('contents'),
                                    tableOutput('tabledata') # Results table
                                  ) # mainPanel()
                                  
                         ), #tabPanel(), Home
                         
                         tabPanel("About", 
                                  titlePanel("About"), 
                                  div(includeMarkdown("about.md"), 
                                      align="justify")
                         ) #tabPanel(), About
                         
              ) # navbarPage()
) # fluidPage()

