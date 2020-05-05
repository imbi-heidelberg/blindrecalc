if (!require("shiny"))
  install.packages("shiny")
if (!require("devtools"))
  install.packages("devtools")
if (!require("blindrecalc"))
  devtools::install_github("imbi-heidelberg/blindrecalc")

library(shiny)
library(blindrecalc)

max_sliders <- 10

ui <- navbarPage("blindrecalc",
                 tabPanel("t-test",
                          fluidPage(
                            titlePanel("Blinded Sample Size Recalculation for Student's t-Test"),
                            fluidRow(
                              column(3,
                                     sliderInput("t_alpha",
                                                 "Maximal Type One Error Rate",
                                                 min = .01,
                                                 max = .1,
                                                 step = .005,
                                                 value = .025),
                                     sliderInput("t_beta",
                                                 "Maximal Type Two Error Rate",
                                                 min = .05,
                                                 max = .4,
                                                 step = .01,
                                                 value = .2),
                                     sliderInput("t_r",
                                                 "Allocation Ratio",
                                                 min = 1/5,
                                                 max = 5,
                                                 step = .1,
                                                 value = 1),
                                     sliderInput("t_delta",
                                                 "Alternative Effect Size",
                                                 min = .05,
                                                 max = 1,
                                                 value = .5,
                                                 step = .01),
                                     sliderInput("t_delta_NI",
                                                 "Non-Inferiority Margin",
                                                 min = 0,
                                                 max = 1,
                                                 value = 0,
                                                 step = .01),
                                     sliderInput("t_n_max",
                                                 "Maximal Sample Size",
                                                 min = 10,
                                                 max = 10000,
                                                 value = 10000,
                                                 step = 2)
                                     ),
                              column(3,
                                     sliderInput("t_n1",
                                                 "Sample Size of Internal Pilot Study",
                                                 min = 1,
                                                 max = 1000,
                                                 value = 10,
                                                 step = 1),
                                     sliderInput("t_nuisance",
                                                 "Assumed Variance",
                                                 min = 0.5,
                                                 max = 20,
                                                 value = 1,
                                                 step = .5),
                                     sliderInput("t_iters",
                                                 "Number of Simulation Iterations",
                                                 min = 1000,
                                                 max = 100000,
                                                 value = 10000,
                                                 step = 1000),
                                     submitButton("Compute!")
                                     ),
                              column(6,
                                     "sample size distribution",
                                     plotOutput(outputId = "t_dist_plot"),
                                     verbatimTextOutput("t_toer"),
                                     verbatimTextOutput("t_pow"),
                                     verbatimTextOutput("t_n_fix"))
                              ),
                            p("This shiny app was developed by:"),
                            a(href="https://www.klinikum.uni-heidelberg.de/medizinische-biometrie/wir-ueber-uns/wir-ueber-uns",
                              img(src="Logo-IMBI_ENGL.jpg", height=100))
                            )
                          ),
                 tabPanel("Chi-Square-test",
                          fluidPage(
                            titlePanel("Blinded Sample Size Recalculation for the Chi-Square-Test"),
                            fluidRow(
                              column(3,
                                     sliderInput("cs_alpha",
                                                 "Maximal Type One Error Rate",
                                                 min = .01,
                                                 max = .1,
                                                 step = .005,
                                                 value = .025),
                                     sliderInput("cs_beta",
                                                 "Maximal Type Two Error Rate",
                                                 min = .05,
                                                 max = .4,
                                                 step = .01,
                                                 value = .2),
                                     sliderInput("cs_r",
                                                 "Allocation Ratio",
                                                 min = 1/5,
                                                 max = 5,
                                                 step = .1,
                                                 value = 1),
                                     sliderInput("cs_delta",
                                                 "Alternative Effect Size",
                                                 min = .05,
                                                 max = 1,
                                                 value = .5,
                                                 step = .01),
                                     sliderInput("cs_n_max",
                                                 "Maximal Sample Size",
                                                 min = 10,
                                                 max = 10000,
                                                 value = 10000,
                                                 step = 2)
                              ),
                              column(3,
                                     sliderInput("cs_n1",
                                                 "Sample Size of Internal Pilot Study",
                                                 min = 1,
                                                 max = 1000,
                                                 value = 10,
                                                 step = 1),
                                     sliderInput("cs_nuisance",
                                                 "Assumed Overall Response Rate",
                                                 min = 0.05,
                                                 max = 0.95,
                                                 value = 0.5,
                                                 step = .05),
                                     submitButton("Compute!")
                              ),
                              column(6,
                                     "sample size distribution",
                                     plotOutput(outputId = "cs_dist_plot"),
                                     verbatimTextOutput("cs_toer"),
                                     verbatimTextOutput("cs_pow"),
                                     verbatimTextOutput("cs_n_fix"))
                            ),
                            p("This shiny app was developed by:"),
                            a(href="https://www.klinikum.uni-heidelberg.de/medizinische-biometrie/wir-ueber-uns/wir-ueber-uns",
                              img(src="Logo-IMBI_ENGL.jpg", height=100))
                          )
                 ),
                 tabPanel("Farrington Manning test",
                          fluidPage(
                            titlePanel("Blinded Sample Size Recalculation for the Farrington Manning Test"),
                            fluidRow(
                              column(3,
                                     sliderInput("fm_alpha",
                                                 "Maximal Type One Error Rate",
                                                 min = .01,
                                                 max = .1,
                                                 step = .005,
                                                 value = .025),
                                     sliderInput("fm_beta",
                                                 "Maximal Type Two Error Rate",
                                                 min = .05,
                                                 max = .4,
                                                 step = .01,
                                                 value = .2),
                                     sliderInput("fm_r",
                                                 "Allocation Ratio",
                                                 min = 1/5,
                                                 max = 5,
                                                 step = .1,
                                                 value = 1),
                                     sliderInput("fm_delta",
                                                 "Alternative Effect Size",
                                                 min = .05,
                                                 max = 1,
                                                 value = .5,
                                                 step = .01),
                                     sliderInput("fm_delta_NI",
                                                 "Non-Inferiority Margin",
                                                 min = 0.05,
                                                 max = 1,
                                                 value = 0.1,
                                                 step = .01),
                                     sliderInput("fm_n_max",
                                                 "Maximal Sample Size",
                                                 min = 10,
                                                 max = 10000,
                                                 value = 10000,
                                                 step = 2)
                              ),
                              column(3,
                                     sliderInput("fm_n1",
                                                 "Sample Size of Internal Pilot Study",
                                                 min = 1,
                                                 max = 1000,
                                                 value = 10,
                                                 step = 1),
                                     sliderInput("fm_nuisance",
                                                 "Assumed Overall Response Rate",
                                                 min = 0.05,
                                                 max = 0.95,
                                                 value = 0.5,
                                                 step = .05),
                                     submitButton("Compute!")
                              ),
                              column(6,
                                     "sample size distribution",
                                     plotOutput(outputId = "fm_dist_plot"),
                                     verbatimTextOutput("fm_toer"),
                                     verbatimTextOutput("fm_pow"),
                                     verbatimTextOutput("fm_n_fix"))
                            ),
                            p("This shiny app was developed by:"),
                            a(href="https://www.klinikum.uni-heidelberg.de/medizinische-biometrie/wir-ueber-uns/wir-ueber-uns",
                              img(src="Logo-IMBI_ENGL.jpg", height=100))
                          )
                 )

                 )




server <- function(input, output) {

  # outpur for t-test

  output$t_dist_plot <- renderPlot({

    d <- blindrecalc::setupStudent(alpha = input$t_alpha,
                                   beta = input$t_beta,
                                   r = input$t_r,
                                   delta = input$t_delta,
                                   delta_NI = input$t_delta_NI,
                                   alternative = "greater",
                                   n_max = input$t_n_max)


    withProgress(message = "Computing characteristics", {

    return(blindrecalc::n_dist(d, input$t_n1, input$t_nuisance, summary = FALSE,
                               plot = TRUE, iters = input$t_iters, seed = 2020))
    })

  })

  output$t_toer <- renderText({
    d <- blindrecalc::setupStudent(alpha = input$t_alpha,
                                   beta = input$t_beta,
                                   r = input$t_r,
                                   delta = input$t_delta,
                                   delta_NI = input$t_delta_NI,
                                   alternative = "greater",
                                   n_max = input$t_n_max)

    return(paste("type I error rate:", blindrecalc::toer(d, input$t_n1, input$t_nuisance,
                                                         recalculation = TRUE,
                                                         iters = input$t_iters, seed = 2020)))
  })

  output$t_pow <- renderText({
    d <- blindrecalc::setupStudent(alpha = input$t_alpha,
                                   beta = input$t_beta,
                                   r = input$t_r,
                                   delta = input$t_delta,
                                   delta_NI = input$t_delta_NI,
                                   alternative = "greater",
                                   n_max = input$t_n_max)

    return(paste("power:", blindrecalc::pow(d, input$t_n1, input$t_nuisance,
                                                         recalculation = TRUE,
                                                         iters = input$t_iters, seed = 2020)))
  })


  output$t_n_fix <- renderText({
    d <- blindrecalc::setupStudent(alpha = input$t_alpha,
                                   beta = input$t_beta,
                                   r = input$t_r,
                                   delta = input$t_delta,
                                   delta_NI = input$t_delta_NI,
                                   alternative = "greater",
                                   n_max = input$t_n_max)

    return(paste("fixed sample size:", ceiling(blindrecalc::n_fix(d, input$t_nuisance))))
  })


  # output for chi-square test

  output$cs_dist_plot <- renderPlot({

    d <- blindrecalc::setupChiSquare(alpha = input$cs_alpha,
                                     beta = input$cs_beta,
                                     r = input$cs_r,
                                     delta = input$cs_delta,
                                     n_max = input$cs_n_max)


    withProgress(message = "Computing characteristics", {

      return(blindrecalc::n_dist(d, input$cs_n1, input$cs_nuisance, summary = FALSE,
                                 plot = TRUE))
    })

  })

  output$cs_toer <- renderText({
    d <- blindrecalc::setupChiSquare(alpha = input$cs_alpha,
                                     beta = input$cs_beta,
                                     r = input$cs_r,
                                     delta = input$cs_delta,
                                     n_max = input$cs_n_max)

    return(paste("type I error rate:", blindrecalc::toer(d, input$cs_n1, input$cs_nuisance,
                                                         recalculation = TRUE,
                                                         iters = input$cs_iters, seed = 2020)))
  })

  output$cs_pow <- renderText({
    d <- blindrecalc::setupChiSquare(alpha = input$cs_alpha,
                                     beta = input$cs_beta,
                                     r = input$cs_r,
                                     delta = input$cs_delta,
                                     n_max = input$cs_n_max)

    return(paste("power:", blindrecalc::pow(d, input$cs_n1, input$cs_nuisance,
                                            recalculation = TRUE,
                                            iters = input$cs_iters, seed = 2020)))
  })


  output$cs_n_fix <- renderText({
    d <- blindrecalc::setupChiSquare(alpha = input$cs_alpha,
                                     beta = input$cs_beta,
                                     r = input$cs_r,
                                     delta = input$cs_delta,
                                     n_max = input$cs_n_max)

    return(paste("fixed sample size:", ceiling(blindrecalc::n_fix(d, input$cs_nuisance))))
  })





  # output for Farrington Manning test

  output$fm_dist_plot <- renderPlot({

    d <- blindrecalc::setupChiSquare(alpha = input$fm_alpha,
                                     beta = input$fm_beta,
                                     r = input$fm_r,
                                     delta = input$fm_delta,
                                     delta_NI = input$fm_delta_NI,
                                     n_max = input$fm_n_max)


    withProgress(message = "Computing characteristics", {

      return(blindrecalc::n_dist(d, input$fm_n1, input$fm_nuisance, summary = FALSE,
                                 plot = TRUE))
    })

  })

  output$fm_toer <- renderText({
    d <- blindrecalc::setupChiSquare(alpha = input$fm_alpha,
                                     beta = input$fm_beta,
                                     r = input$fm_r,
                                     delta = input$fm_delta,
                                     delta_NI = input$fm_delta_NI,
                                     n_max = input$fm_n_max)

    return(paste("type I error rate:", blindrecalc::toer(d, input$fm_n1, input$fm_nuisance,
                                                         recalculation = TRUE,
                                                         iters = input$fm_iters, seed = 2020)))
  })

  output$fm_pow <- renderText({
    d <- blindrecalc::setupChiSquare(alpha = input$fm_alpha,
                                     beta = input$fm_beta,
                                     r = input$fm_r,
                                     delta = input$fm_delta,
                                     delta_NI = input$fm_delta_NI,
                                     n_max = input$fm_n_max)

    return(paste("power:", blindrecalc::pow(d, input$fm_n1, input$fm_nuisance,
                                            recalculation = TRUE,
                                            iters = input$fm_iters, seed = 2020)))
  })


  output$fm_n_fix <- renderText({
    d <- blindrecalc::setupChiSquare(alpha = input$fm_alpha,
                                     beta = input$fm_beta,
                                     r = input$fm_r,
                                     delta = input$fm_delta,
                                     delta_NI = input$fm_delta_NI,
                                     n_max = input$fm_n_max)

    return(paste("fixed sample size:", ceiling(blindrecalc::n_fix(d, input$fm_nuisance))))
  })



}

# Run the application
shinyApp(ui = ui, server = server)

