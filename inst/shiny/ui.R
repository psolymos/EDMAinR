## note: the ui_* UI pieces can be stored & sourced as separate file
## or added as modules

## home, docs, settings ---------------

ui_home <- fluidRow(
  column(width=12,
    h2("Welcome!")
  )
)

ui_docs <- fluidRow(
  column(width=12,
    h2("Docs")
  )
)

ui_settings <- fluidRow(
  column(width=12,
    h2("Data sets"),
    fileInput("file_num", "Choose xyz file for numerator",
                     multiple = FALSE,
                     accept = c("text/plain", ".xyz")),
    fileInput("file_den", "Choose xyz file for denominator",
                     multiple = FALSE,
                     accept = c("text/plain", ".xyz")),
#    fileInput("file_surf", "Choose surf file (optional)",
#                     multiple = FALSE,
#                     accept = c("text/plain", ".xyz")),
    h2("Settings"),
    column(width = 6,
        sliderInput("B", "Number of bootstrap replicates", 50, 1000, 100, 50),
        sliderInput("level", "Confidence level (%)", 80, 100, 90, 1)
    ),
    column(width = 6,
        radioButtons("refdenom", "Reference sample",
                     c("Numerator", "Denominator"), "Denominator"),
        radioButtons("mix", "Resampling type",
                     c("Mixed", "Fixed reference"), "Fixed reference")
    )
  )
)

## numerator --------------------

ui_tab_data_xyz_num <- fluidRow(
  column(width=12,
    h2("Numerator"),
    h3("Summary"),
    verbatimTextOutput("num_summary1"),
    h3("Landmarks"),
    verbatimTextOutput("num_summary2")
  )
)
ui_tab_data_ord_num <- fluidRow(
  column(width=12,
    h2("Numerator"),
    h3("Ordination"),
    plotOutput("num_ord")
  )
)
ui_tab_fm_num <- fluidRow(
  column(width=12,
    h2("Numerator"),
    h3("Form matrix"),
    reactableOutput("num_fm")
  )
)
ui_tab_fm_3d_num <- fluidRow(
  column(width=12,
    h2("Numerator"),
    h3("3D plot"),
    rglwidgetOutput("num_3d")
  )
)

## denominator --------------------

ui_tab_data_xyz_den <- fluidRow(
  column(width=12,
    h2("Denumerator"),
    h3("Summary"),
    verbatimTextOutput("den_summary1"),
    h3("Landmarks"),
    verbatimTextOutput("den_summary2")
  )
)
ui_tab_data_ord_den <- fluidRow(
  column(width=12,
    h2("Denumerator"),
    h3("Ordination"),
    plotOutput("den_ord")
  )
)
ui_tab_fm_den <- fluidRow(
  column(width=12,
    h2("Denumerator"),
    h3("Form matrix"),
    reactableOutput("den_fm")
  )
)
ui_tab_fm_3d_den <- fluidRow(
  column(width=12,
    h2("Denumerator"),
    h3("3D plot"),
    rglwidgetOutput("den_3d")
  )
)

## FDM --------------------

ui_tab_fdm_tobs <- fluidRow(
  column(width=12,
    h2("FDM"),
    h3("Global testing"),
    plotOutput("fdm_tobs")
  )
)
ui_tab_fdm_ci <- fluidRow(
  column(width=12,
    h2("FDM"),
    h3("Local testing / CI"),
    reactableOutput("fdm_ci")
  )
)
ui_tab_fdm_ord <- fluidRow(
  column(width=12,
    h2("FDM"),
    h3("Ordination"),
    plotOutput("fdm_ord")
  )
)
ui_tab_fdm_3d <- fluidRow(
  column(width=12,
    h2("FDM"),
    h3("3D plot"),
    rglwidgetOutput("fdm_3d")
  )
)

## Layout

dashboardPage(
  dashboardHeader(title = paste("EDMA in R", ver[1])),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon=icon("home")),
      menuItem("Settings", tabName = "settings", icon=icon("cog")),
      menuItem("Numerator", tabName = "numerator",
               icon=icon("chevron-circle-up"),
        menuSubItem("XYZ data", tabName = "tab_data_xyz_num"),
        menuSubItem("Ordination", tabName = "tab_data_ord_num"),
        menuSubItem("Form matrix", tabName = "tab_fm_num"),
        menuSubItem("3D plot", tabName = "tab_fm_3d_num")
      ),
      menuItem("Denominator", tabName = "denominator",
               icon=icon("chevron-circle-down"),
        menuSubItem("XYZ data", tabName = "tab_data_xyz_den"),
        menuSubItem("Ordination", tabName = "tab_data_ord_den"),
        menuSubItem("Form matrix", tabName = "tab_fm_den"),
        menuSubItem("3D plot", tabName = "tab_fm_3d_den")
      ),
      menuItem("Form difference", tabName = "fm",
               icon=icon("border-all"),
        menuSubItem("Global testing", tabName = "tab_fdm_tobs"),
        menuSubItem("Local testing", tabName = "tab_fdm_ci"),
        menuSubItem("Ordination", tabName = "tab_fdm_ord"),
        menuSubItem("3D plot", tabName = "tab_fdm_3d")
      ),
      menuItem("Documentation", tabName = "docs", icon=icon("book"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem("home", ui_home),
      tabItem("settings", ui_settings),

      tabItem("tab_data_xyz_num", ui_tab_data_xyz_num),
      tabItem("tab_data_ord_num", ui_tab_data_ord_num),
      tabItem("tab_fm_num", ui_tab_fm_num),
      tabItem("tab_fm_3d_num", ui_tab_fm_3d_num),

      tabItem("tab_data_xyz_den", ui_tab_data_xyz_den),
      tabItem("tab_data_ord_den", ui_tab_data_ord_den),
      tabItem("tab_fm_den", ui_tab_fm_den),
      tabItem("tab_fm_3d_den", ui_tab_fm_3d_den),

      tabItem("tab_fdm_tobs", ui_tab_fdm_tobs),
      tabItem("tab_fdm_ci", ui_tab_fdm_ci),
      tabItem("tab_fdm_ord", ui_tab_fdm_ord),
      tabItem("tab_fdm_3d", ui_tab_fdm_3d),

      tabItem("docs", ui_docs)
    )
  )
)

