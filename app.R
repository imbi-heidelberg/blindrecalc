if (!require("shiny"))
  install.packages("shiny")
if (!require("devtools"))
  install.packages("devtools")
if (!require("blindrecalc"))
  devtools::install_github("imbi-heidelberg/blindrecalc")

library(shiny)
library(blindrecalc)

ui <- fluidPage(
   titlePanel("Blinded Sample Size Recalculation With the Package blindrecalc"),

   fluidRow(

     column(3,
            sliderInput("alpha",
                        "Maximal Type One Error Rate",
                               min = .01,
                               max = .1,
                               step = .005,
                               value = .025),
                   sliderInput("beta",
                               "Maximal Type Two Error Rate",
                               min = .05,
                               max = .4,
                               step = .01,
                               value = .2),
                   sliderInput("r",
                               "Allocation Ratio",
                               min = 1/5,
                               max = 5,
                               step = .1,
                               value = 1),
                   sliderInput("delta",
                               "Alternative Effect Size",
                               min = .05,
                               max = 1,
                               value = .5,
                               step = .01),
                   sliderInput("delta_NI",
                               "Non-Inferiority Margin",
                               min = 0,
                               max = 1,
                               value = 0,
                               step = .01),
                   sliderInput("n_max",
                               "Maximal Sample Size",
                               min = 10,
                               max = 10000,
                               value = 10000,
                               step = 2)
     ),
     column(3,
                   sliderInput("n1",
                               "Sample Size of Internal Pilot Study",
                               min = 1,
                               max = 1000,
                               value = 10,
                               step = 1),
                   sliderInput("nuisance",
                               "Assumed Variance",
                               min = 0.5,
                               max = 20,
                               value = 1,
                               step = .5),
                   sliderInput("iters",
                               "Number of Simulation Iterations",
                               min = 1000,
                               max = 100000,
                               value = 10000,
                               step = 1000),
                   submitButton("Compute!")
     ),
     column(6,
            "sample size distribution",
            plotOutput(outputId = "dist_plot"),
            verbatimTextOutput("toer"),
            verbatimTextOutput("pow"),
            verbatimTextOutput("n_fix"))


      ),

   p("This shiny app was developed by:"),
   a(href="https://www.klinikum.uni-heidelberg.de/medizinische-biometrie/wir-
ueber-uns/wir-ueber-uns", img(src="Logo-IMBI_ENGL.jpg", height=100))
)




server <- function(input, output) {



  output$dist_plot <- renderPlot({

    d <- blindrecalc::setupStudent(alpha = input$alpha,
                                   beta = input$beta,
                                   r = input$r,
                                   delta = input$delta,
                                   delta_NI = input$delta_NI,
                                   alternative = "greater",
                                   n_max = input$n_max)


    withProgress(message = "Computing characteristics", {

    return(blindrecalc::n_dist(d, input$n1, input$nuisance, summary = FALSE,
                               plot = TRUE, iters = input$iters, seed = 2020))
    })

  })

  output$toer <- renderText({
    d <- blindrecalc::setupStudent(alpha = input$alpha,
                                   beta = input$beta,
                                   r = input$r,
                                   delta = input$delta,
                                   delta_NI = input$delta_NI,
                                   alternative = "greater",
                                   n_max = input$n_max)

    return(paste("type I error rate:", blindrecalc::toer(d, input$n1, input$nuisance,
                                                         recalculation = TRUE,
                                                         iters = input$iters, seed = 2020)))
  })

  output$pow <- renderText({
    d <- blindrecalc::setupStudent(alpha = input$alpha,
                                   beta = input$beta,
                                   r = input$r,
                                   delta = input$delta,
                                   delta_NI = input$delta_NI,
                                   alternative = "greater",
                                   n_max = input$n_max)

    return(paste("power:", blindrecalc::pow(d, input$n1, input$nuisance,
                                                         recalculation = TRUE,
                                                         iters = input$iters, seed = 2020)))
  })


  output$n_fix <- renderText({
    d <- blindrecalc::setupStudent(alpha = input$alpha,
                                   beta = input$beta,
                                   r = input$r,
                                   delta = input$delta,
                                   delta_NI = input$delta_NI,
                                   alternative = "greater",
                                   n_max = input$n_max)

    return(paste("fixed sample size:", ceiling(blindrecalc::n_fix(d, input$nuisance))))
  })



}

# Run the application
shinyApp(ui = ui, server = server)

