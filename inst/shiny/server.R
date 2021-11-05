server <- function(input, output, session) {

    ## file uploads ----------------------
    data_num <- reactive({
        file <- input$file_num
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == "xyz", "Please upload an xyz file for numerator"))
        read_xyz(file$datapath)
    })
    data_den <- reactive({
        file <- input$file_den
        ext <- tools::file_ext(file$datapath)
        req(file)
        validate(need(ext == "xyz", "Please upload an xyz file for denominator"))
        read_xyz(file$datapath)
    })
#    data_surf <- reactive({
#        file <- input$file_surf
#        ext <- tools::file_ext(file$datapath)
#        req(file)
#        validate(need(ext == "surf", "Please upload a surf file"))
#        read.hxsurf(file$datapath)
#    })

    ## numerator --------------------
    output$num_summary1 <- renderPrint({
        req(data_num())
        data_num()
    })
    output$num_summary2 <- renderPrint({
        req(data_num())
        landmarks(data_num())
    })
    output$num_ord <- renderPlot({
        req(data_num())
        plot_ord(data_num())
    })
    num_edma_fit <- reactive({
        req(data_num(), input$B)
        op <- pboptions(label="Numerator FM")
        on.exit(pboptions(op))
        edma_fit(data_num(), B=input$B)
    })
    output$num_fm <- renderReactable({
        req(num_edma_fit())
        reactable(Meanform(num_edma_fit()), filterable=TRUE)
    })
    output$num_3d <- renderRglwidget({
        req(num_edma_fit())
        plot_3d(num_edma_fit(), labels=TRUE)
        rglwidget()
    })

    ## denominator --------------------
    output$den_summary1 <- renderPrint({
        req(data_den())
        data_den()
    })
    output$den_summary2 <- renderPrint({
        req(data_den())
        landmarks(data_den())
    })
    output$den_ord <- renderPlot({
        req(data_den())
        plot_ord(data_den())
    })
    den_edma_fit <- reactive({
        req(data_den(), input$B)
        op <- pboptions(label="Denominator FM")
        on.exit(pboptions(op))
        edma_fit(data_den(), B=input$B)
    })
    output$den_fm <- renderReactable({
        req(den_edma_fit())
        reactable(Meanform(den_edma_fit()), filterable=TRUE)
    })
    output$den_3d <- renderRglwidget({
        req(den_edma_fit())
        plot_3d(den_edma_fit(), labels=TRUE)
        rglwidget()
    })

    ## FDM --------------------
    fdm <- reactive({
        req(num_edma_fit(), den_edma_fit(), input$B)
        op <- pboptions(label="FDM")
        on.exit(pboptions(op))
        edma_fdm(numerator=num_edma_fit(),
                denominator=den_edma_fit(),
                ref_denom = input$refdenom == "Denominator",
                mix=input$mix == "Mixed",
                B=input$B)
    })
    output$fdm_ord <- renderPlot({
        req(fdm())
        plot_ord(fdm())
    })
    output$fdm_tobs <- renderPlot({
        req(fdm())
        plot_test(fdm())
    })
    output$fdm_ci <- renderReactable({
        req(fdm())
        reactable(get_fdm(fdm(), sort=TRUE), filterable=TRUE)
    })
    output$fdm_3d <- renderRglwidget({
        req(fdm())
        plot_3d(fdm(), labels=TRUE, signif_only=FALSE, q=0.1)
#        if (isTruthy(data_surf())) {
#            xyz001 <- as.matrix(data_num()[,,1])
#            rownames(xyz001) <- landmarks(data_num())
#            plot3d(data_surf(), col="#FFF8DC", add=TRUE, alpha=0.2)
#        }
        rglwidget()
    })


}

