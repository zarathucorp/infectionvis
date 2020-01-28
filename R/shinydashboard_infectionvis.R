
#' @title fluODE equation
#' @description fluODE equation
#' @param time time
#' @param state initial states
#' @param params parameters
#' @return fluODE output
#' @details fluODE
#' @examples 
#' init    <- c(S = 49656776, E = 0, I = 1, A = 0, R = 0)
#' params  <- c(beta = 0.4726, delta = 1/2, p = 2/3, kappa = 1.9, alpha = 1/6, eta = 1/6, q = 1/2)
#' times   <- seq(0, 141, by = 1)
#' deSolve::ode(y = init, times = times, func = fluODE, parms = params)
#' @rdname fluODE
#' @export 

fluODE <- function(time, state, params) {
  with(as.list(c(state, params)), {
    N <- S + E + I + A + R
    dS <- -beta*S*(q*I + delta*A)/N
    dE <- beta*S*(q*I + delta*A)/N - kappa*E
    dI <- kappa*(p)*E - alpha*I
    dA <- kappa*(1-p)*E - eta*A
    dR <- alpha*I + eta*A
    return(list(c(dS, dE, dI, dA, dR)))
  })
}




#' @title 'Shinydashboard' for Infection Simulation
#' @description RStudio addin- 'Shinydashboard' for Infection Simulation 
#' @param max.filesize Maximum file size to upload (MB), Default: 2048 (2 GB)
#' @return RStudio addin- 'Shinydashboard' for Infection Simulation
#' @details RStudio addin- 'Shinydashboard' for Infection Simulation
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[deSolve]{ode}}
#' @rdname infectionvis
#' @export 
#' @importFrom deSolve ode
#' @importFrom shinycustomloader withLoader
#' @importFrom shinydashboard box
#' @importFrom data.table fread
#' @import shiny 
#' @import shinydashboard
#' @import highcharter

infectionvis <- function(max.filesize = 2048){
  options(shiny.maxRequestSize = max.filesize * 1024^2)

  
  ui <- dashboardPage(
    dashboardHeader(
      title = "Infection Dashboard",
      dropdownMenu(type = "messages",
                   messageItem(
                     from = "Sales Dept",
                     message = "Sales are steady this month."
                   ),
                   messageItem(
                     from = "New User",
                     message = "How do I register?",
                     icon = icon("question"),
                     time = "13:45"
                   ),
                   messageItem(
                     from = "Support",
                     message = "The new server is ready.",
                     icon = icon("life-ring"),
                     time = "2014-12-01"
                   )
      )),
    
    dashboardSidebar(
      withMathJax(),
      sidebarMenu(
        fileInput("file", "Upload zip file"),
        helpText(a(h4("Example data", class = "btn btn-default action-button" , 
                      style = "fontweight:600"), href="https://zarathu.s3.ap-northeast-2.amazonaws.com/exampledata/fludata09.csv")),
        uiOutput("inicond"),
        menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
        menuItem("Other", icon = icon("th"), tabName = "other",
                 badgeLabel = "new", badgeColor = "green")
      )
    ),
    dashboardBody(
      tabItems(
        # First tab content
        tabItem(tabName = "dashboard",
                fluidRow(
                  box(title = "Parameters", status = "primary", solidHeader = TRUE, collapsible = T, uiOutput("parameter"), width = 12) 
                ),
                fluidRow(
                  box(title = "Estimation", status = "primary", solidHeader = TRUE, collapsible = T,
                      withLoader(highchartOutput(paste0("chart", "Estimation")), type="html", loader="loader6"))
                )
                
        )
      )
    )
  )
  

  
  server <- function(input, output, session){
    

    userFile <- eventReactive(input$file, {
      # If no file is selected, don't do anything
      #validate(need(input$file, message = FALSE))
      input$file
    })
    
    incidata <- eventReactive(input$file, {
      validate(need((grepl("csv", userFile()$name) == T), message = "Please upload zip file"))
      out <- fread(userFile()$datapath)
      return(out)
    })
    
    
    output$inicond <- renderUI({
      tagList(
        h4("Initial conditions"),
        numericInput("ini_S", "S", value = 49656776, min = 0, max = 100000000),
        sliderInput("ini_E", "E", value = 0, min = 0, max = 1),
        sliderInput("ini_I", "I", value = 1, min = 0, max = 1),
        sliderInput("ini_A", "A", value = 0, min = 0, max = 1),
        sliderInput("ini_R", "R", value = 0, min = 0, max = 1)
      )
    })
    
    
    output$parameter <- renderUI({
      tagList(
        withMathJax(),
        column(12/7, numericInput("param_beta", '$$\\beta$$', value = 0.4726, min = 0, step = 0.01)),
        column(12/7, numericInput("param_delta", "$$\\delta$$", value = 1/2, min = 0)),
        column(12/7, numericInput("param_p", "$$p$$", value = 2/3, min = 0)),
        column(12/7, numericInput("param_kappa", "$$\\kappa$$", value = 1.9, min = 0)),
        column(12/7, numericInput("param_alpha", "$$\\alpha$$", value = 1/6, min = 0)),
        column(12/7, numericInput("param_eta", "$$\\eta$$", value = 1/6, min = 0)),
        column(12/7, numericInput("param_q", "$$q$$", value = 1/2, min = 0))
        
      )
    })
    
    
    
    out.ode <- reactive({
      req(incidata())
      params  <- c(beta = input$param_beta ,delta = input$param_delta, p = input$param_p, kappa = input$param_kappa, alpha = input$param_alpha, eta = input$param_eta, q = input$param_q)
      init <-  c(S = input$ini_S, E = input$ini_E, I = input$ini_I, A = input$ini_A, R = input$ini_R)
      times   <- seq(0, 141, by = 1)
      output <- deSolve::ode(y = init, times = times, func = fluODE, parms = params)
      output <- as.data.frame(output)
      ## adding cumulative incidence
      output$cuminci <- cumsum(output$E * input$param_p * input$param_kappa)
      return(output)
    })
    
    
    output$chartEstimation <- renderHighchart({
      highchart() %>% 
        hc_xAxis(time = out.ode()$time) %>% 
        hc_add_series(name = "Data", data = incidata()$V1, type = "scatter") %>% 
        hc_add_series(name = "Cumulative incidence", data = out.ode()$cuminci) %>% 
        hc_exporting(enabled = T) %>% 
        hc_tooltip(valueDecimals = 0)
      
    })
    
    
    session$onSessionEnded(function() {
      stopApp()
    })
    
  }
  
  
  
  viewer <- browserViewer(browser = getOption("browser"))
  runGadget(ui, server, viewer = viewer)
}

