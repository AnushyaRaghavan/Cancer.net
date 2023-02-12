
library(leaflet)
library(dplyr)
library(BCRA)
library(flexdashboard)
library(stringr)
library(png)
library(shinyjs)
library(DT)
library(visNetwork)
library(rintrojs)

char_race <- c("Asian","Non-Asian")

deter_race <- function(x) {
  return(char_race[as.numeric(x)])
}


# Choices for drop-downs
vars <- c(
    "All cancer site"="all_cancer_sites",
    "Tongue"="tongue",
    "Mouth"="mouth",
    "Breast"="breast",
    "Lungs"="lungs",
    "Stomach"="stomach",
    "Bladder"="bladder"
)
vars1<- c(
  "Male"="M",
  "Female"="F"
  
)


source("carouselPanel.R")

# Panel div for visualization
# override the currently broken definition in shinyLP version 1.1.0
panel_div <- function(class_type, content) {
  div(class = sprintf("panel panel-%s", class_type),
      div(class = "panel-body", content)
  )
}
navbarPage(title = img(src="logo.jpg", height = "30px"), id = "navBar",
                   theme = "paper.css",
                   collapsible = TRUE,
                   inverse = TRUE,
                   windowTitle = "Cancer Data",
                   position = "fixed-top",
                   footer = includeHTML("./www/include_footer.html"),
                   header = tags$style(
                     ".navbar-right {
                       float: right !important;
                       }",
                     "body {padding-top: 75px;}"),
           
           
                   tabPanel("HOME", value = "home",
                    
                    shinyjs::useShinyjs(),
                    
                    tags$head(tags$script(HTML('
                                                       var fakeClick = function(tabName) {
                                                       var dropdownList = document.getElementsByTagName("a");
                                                       for (var i = 0; i < dropdownList.length; i++) {
                                                       var link = dropdownList[i];
                                                       if(link.getAttribute("data-value") == tabName) {
                                                       link.click();
                                                       };
                                                       }
                                                       };
                                                       '))),
                    fluidRow(
                      HTML("
                                     
                                     <section class='banner'>
                                     
                                     <br>  <br>  <br>  <br> <br> <br> <br>
                                     <br> <br>
                                     <br>  <br> <br>  <br> <br>  <br> <br>
                                     <p class='parallax_description'>You Are Strong Enough To fight This Through.</p>
                                     </section>
                                     ")
                    ),
                    
                    #                                      <h2 class='parallax'>CANCER</h2>

                    tags$hr(),
                    
                    # WHAT
                    fluidRow(
                      column(3),
                      column(6,
                             shiny::HTML("<br><br><center> <h1>What is Cancer</h1> </center><br>"),
                             shiny::HTML("<center><h5>Our body has different types of cells which grow and
                                       multiply several times in a day to maintain normal function.
                                       Occasionally this cell division is affected and creates a wrong
                                       cell or an abnormal cell, these divide further in an uncontrolled
                                       way affecting the other parts of the body.
                                       These cells which divide abnormally are called cancers.
                                       These cells can locally infiltrate the adjacent organs
                                       or invade the blood and lymphatic system and cause spread
                                                   or metastasis</h5><center>"),
                             shiny::HTML("<center><h5>Cancer can occur in almost
                                                   all parts of the body and they are categorized mainly from
                                                   the cell in which they begin.</h5> </center><br>"),
                             shiny::HTML("<h6><li>Cancer arising from epithelium or covering of an organ - CARCINOMA.</li>
                                                         <li>Cancer arising from connective or supportive cells (E.g.: bone, muscle,) - SARCOMA.</li>
                                                         <li>Cancer arising from cells which produce blood(bone marrow) - LEUKAEMIA.</li>
                                                         <li>Cancer arising from cells of immune system - LYMPHOMA and MYELOMA</li>
                                                         <li>Cancer arising from nervous system - ASTROCYTOMA, GLIOMA.</li></h6>")
                      ),
                      column(3)
                    ),
                    
                    fluidRow(
                      
                      style = "height:50px;"),
                    
                    # PAGE BREAK
                    tags$hr(),
                    
                    # HOW
                    fluidRow(
                      column(3),
                      column(6,
                             shiny::HTML("<br><br><center> <h1>How Does Cancer Start</h1> </center><br>"),
                             shiny::HTML("<h5><li>By damage to DNA - the genetic material in the cell.</li>
                                                         <li>By interfering with cell division and other controlling mechanisms of the cell.</li>
                                                         <li>Damage to DNA can be caused by certain hydrocarbons, chemicals which are present
                                                         in tobacco for example, radiation or due to certain chemical action in our body when
                                                         the built in protection mechanism fails in our body.</li>
                                                         <li>Certain substances for example stimulate cell division causing damage to DNA and causing mutation.</li>
                                                         <li>Any factor which leads to the above can produce cancer.</li></h5>")
                      ),
                      column(3)
                    ),
                    tags$hr(),
                    
                    fluidRow(
                      column(3),
                      column(6,
                             
                             shiny::HTML("<br><br><center> <h1>What Do We Understand By Stages of Cancer
                              </h1> </center>"),
                             shiny::HTML("<br><br><center> <h5>Staging is done when the patient is first 
                                          seen with the particular cancer.Staging is important because it
                                          decides the treatment to be given and also gives an indication 
                                          of the outcome.</h5> </center><br>"))),
                    
                    fluidRow(
                      column(2),
                      
                      column(2,
                             div(class="panel panel-default", 
                                 div(class="panel-body",  width = "600px",
                                     align = "center",
                                     div(
                                       tags$img(src = "zero.png", 
                                                width = "50px", height = "50px")
                                     ),
                                     div(
                                       h5(
                                         "Stage 0 - pre cancerous stage."
                                       )
                                     )
                                 )
                             )
                      ),
                      column(2,
                             div(class="panel panel-default",
                                 div(class="panel-body",  width = "600px", 
                                     align = "center",
                                     div(
                                       tags$img(src = "one.svg", 
                                                width = "50px", height = "50px")                                             ),
                                     div(
                                       h5(
                                         "Stage 1 - localized cancer spread."
                                       )
                                     )
                                 )
                             )
                      ),
                      column(2,
                             div(class="panel panel-default",
                                 div(class="panel-body",  width = "600px", 
                                     align = "center",
                                     div(
                                       tags$img(src = "twoth.png",
                                                width = "100px", height = "50px")),
                                     
                                     div(
                                       
                                       h5(
                                         "Stage 2 and 3 - Regional spread."
                                       )
                                     )
                                 )
                             )
                      ),
                      column(2,
                             div(class="panel panel-default",
                                 div(class="panel-body",  width = "600px", 
                                     align = "center",
                                     div(
                                       tags$img(src = "four.svg", 
                                                width = "50px", height = "50px")),
                                     div(
                                       h5(
                                         "Stage 4 -Distant Spread Cancer."
                                       )
                                     )
                                 )
                             )
                      ),
                      column(2)
                      
                    ),
                    
                    
                    
                    fluidRow(
                      
                      style = "height:50px;"),
                    
                    # PAGE BREAK
                    tags$hr(),
                    
                    # WHERE
                    fluidRow(
                      column(3),
                      column(6,
                             shiny::HTML("<br><br><center> <h1>Can We Prevent Cancer</h1> </center><br>"),
                             shiny::HTML("<h5>Yes, some common cancers are preventable.Leading a healthy
                                                   life certainly takes you a long way.Avoid smoking, tobacco chewing
                                                   which can reduce the risk of oral cavity and lung cancers.Limit 
                                                   taking alcohol, exercise regularly and maintaining healthy body
                                                   weight as overweight and obesity are risk factors for breast and
                                                   endometrial cancers.Cervical smear and mammogram after the age of 
                                                   35 years or earlier if needed as per your doctor's advice.</h5>")
                      ),
                      column(3)
                    ),
                    
                    fluidRow(
                      
                      style = "height:50px;"),
                    
                    # PAGE BREAK
                    tags$hr(),
                    
                    # HOW TO START
                    fluidRow(
                      column(3),
                      column(6,
                             shiny::HTML("<br><br><center> <h1>What is Genetic Risk</h1> </center><br>"),
                             shiny::HTML("<h5>It is the likelihood of developing cancer
                                       depending upon your family history.There are certain genes
                                       in our body which protect us from developing cancer, when 
                                       these genes are faulty that individual is at higher than 
                                       average risk.Only a small percentage (up to 5 per cent) of
                                       certain types of cancer is due to people being born with a
                                       faulty gene.Gene with the fault is inherited and sometimes
                                       it occurs for unknown reasons.Not all those who have a mutated
                                       cancer protection gene will develop cancer, but they are at a 
                                       higher risk.</h5>")
                      ),
                      column(3)
                    ),
                    
                    # BUTTONS TO START
                    
                    fluidRow(
                      
                      style = "height:50px;"),
                    
                    # PAGE BREAK
                    
                    # INSTRUCTIONAL SECTION
                    
                    # Embedded Video from Vimeo on how to use this tool
                    # fluidRow(
                    #     column(3),
                    #     column(6,
                    #            tags$embed(src = "https://player.vimeo.com/video/8419440",
                    #                       width = "640", height = "360") 
                    #     ),
                    #     column(3)
                    # ),
                    
                    fluidRow(
                      
                      style = "height:50px;"),
                    
                    # PAGE BREAK
                    tags$hr(),
                    
                    # AFTERWARD
                    
                    
                    fluidRow(
                      
                      style = "height:50px;"),
      

                    fluidRow(style = "height:25px;"
                    )
                    
           ),
           
           
           #interactive mapp================
           
           
           tabPanel("INTERACTIVE MAP",
                    div(class="outer",
                        
                        tags$head(
                            # Include our custom CSS
                            includeCSS("style.css"),
                            includeScript("gomap.js")
                        ),
                        
                        # If not using custom CSS, set height of leafletOutput to a number instead of percent
                        leafletOutput("map", width="100%", height="100%"),
                        
                        # Shiny versions priosr to 0.11 should use class = "modal" instead.
                        absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                                      draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
                                      width = 330, height = "auto",
                                      
                                      h2("Cancer Data"),
                                      
                                      selectInput("color", "Data", vars),
                                      selectInput("color", "Sex", vars1)
                                     
                                      
                                     
                        ),
                        
                        
                    )
           ),
           # Application title
           tabPanel("PREDICT TOOL",
                    
                    # Sidebar with a slider input for number of bins
                    sidebarLayout(
                      sidebarPanel(
                        sliderInput("Split",
                                    "Minimum Split:",
                                    min = 2,
                                    max = 60,
                                    value = 20),
                        sliderInput("Clump","Clump Size:",
                                    min = 1,max = 10,value = 5),
                        sliderInput("UniCell_Size","Uniformity of Cell Size:",
                                    min = 1,max = 10,value = 5),
                        sliderInput("Uni_CellShape","Uniformity of CellShape:",
                                    min = 1,max = 10,value = 5),
                        sliderInput("MargAdh","Marginal Adhesion:",
                                    min = 1,max = 10,value = 5),
                        sliderInput("SEpith","SEpith:",
                                    min = 1,max = 10,value = 5),
                        sliderInput("BareN","Bare Nuclei:",
                                    min = 1,max = 10,value = 5),
                        sliderInput("BChromatin","Bland Chromatin:",
                                    min = 1,max = 10,value = 5),
                        sliderInput("NoemN","NoemN:",
                                    min = 1,max = 10,value = 5),
                        sliderInput("Mitoses","Mitoses:",
                                    min = 1,max = 10,value = 5)
                        
                      ),
                      
                      # Show a plot of the generated distribution
                      mainPanel(
                        plotOutput("distPlot"),
                        h4(textOutput("Prediction"))
                      )
                    )),
           
           tabPanel("CANCER RISK CALCULATOR",
                   
                      fluidPage(theme = shinytheme("slate"),
                                headerPanel("Breast Cancer Risk Assessment"),
                                fluidRow(
                                  sidebarPanel(
                                    fluidRow(column(10, h2('Questionaire')),
                                             column(2, tags$br(), 
                                                    imageOutput("image", height = 25, width = 25,
                                                                click = clickOpts(id = "image_click")),
                                                    tags$br(), tags$br())),
                                    
                                    numericInput('age',"Current Age?", 35, min=35, max=85, step = 1),
                                    
                                    selectInput('menstruation',"Age of first menstruation?",
                                                choices = list("Less than 12 years old" = 11,
                                                               "12 through 13 years old" = 12,
                                                               "Greater than 13 years old" = 14,
                                                               "Unknown" = 99), selected = 99),
                                    
                                    selectInput('first_birth',"Age at first birth?",
                                                choices = list("No births" = 98,
                                                               "Less than 20 years old"= 19,
                                                               "20 through 24 years old" = 20,
                                                               "25 through 29 years old" = 25,
                                                               "30 years old and greater"= 30,
                                                               "Unknown" = 99), selected = 99),
                                    
                                    selectInput('relatives',"Number of 1st degree relatives that have had breast cancer?",
                                                choices = list("Unknown" = 99,
                                                               "0 relatives" = 0,
                                                               "1 relatives" = 1,
                                                               "2 or more relatives" = 2), selected = 99),
                                    
                                    selectInput('biopsies', "Number of breast biopsies?",
                                                choices = list("Unknown" = 99,
                                                               "no biopsies" = 0,
                                                               "1 biopsies" = 1,
                                                               "2 or more biopsies" = 2), selected = 99),
                                    
                                    conditionalPanel(
                                      condition = "input.biopsies == 1 || input.biopsies == 2",
                                      selectInput("hyperplasia", "Did the biopsy display hyperplasia?",
                                                  choices = list("No" = 0,
                                                                 "Yes" = 1,
                                                                 "Unknown" = 99), selected = 99)),
                                    
                                    selectInput('race',"Race/Ethnicity?",
                                                choices = list("Asian" = 1,
                                                               "Non-Asian" = 2
                                                ), selected = 1),
                                    
                                    actionButton("do","Calculate Risk", class = "btn-primary")
                                    
                                  ),
                                  
                                  mainPanel(
                                    h3('Your Results'),
                                    tags$hr(),
                                    span(textOutput("five_yr_title"), style="font-size: 20px"),
                                    textOutput("five_yr_text"),
                                    tags$br(),
                                    span(textOutput("lifetime_title"), style="font-size: 20px"),
                                    textOutput("lifetime_text"),
                                    tags$hr(),
                                    fluidRow(
                                      column(7
                                             , fluidRow(
                                               column(12, span(textOutput("advice_title"), style="font-size: 20px"),
                                                      tags$hr(width = "50%", align = "left"))
                                             )
                                             , fluidRow(
                                               column(12, textOutput("advice_text1"), tags$br(),
                                                      textOutput("advice_text2"))
                                             )
                                      )
                                      , column(5, tags$br(), htmlOutput("short_five"), htmlOutput("short_life"), tags$br(), gaugeOutput("plt1"))
                                    )
                                  )), ),
                   
         ),
                      
                      
                      
           tabPanel("CANCER DATA EXPLORER",
                    fluidRow(),
                    
                    hr(),
                    DT::dataTableOutput("ziptable")
           ),
           
           conditionalPanel("false", icon("crosshair"))
)
