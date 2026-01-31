#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
})

haplotypes <- c("MF", "Mf", "Mz", "mF", "mf", "mz")

hap_F_allele <- function(h) substr(h, 2, 2)
hap_M_allele <- function(h) substr(h, 1, 1)

normalize_freqs <- function(p) {
  total <- sum(p)
  if (total <= 0) {
    stop("Haplotype frequencies must sum to a positive value.")
  }
  p / total
}

validate_params <- function(p, s_M, s_m, u_m, u_f, u_z, z_reduction, max_generations, tol) {
  if (any(!is.finite(c(s_M, s_m, u_m, u_f, u_z, z_reduction, max_generations, tol)))) {
    stop("Parameters must be numeric.")
  }
  if (max_generations < 1) {
    stop("max_generations must be >= 1.")
  }
  if (s_M < 0 || s_M > 1 || s_m < 0 || s_m > 1 ||
      u_m < 0 || u_m > 1 || u_f < 0 || u_f > 1 || u_z < 0 || u_z > 1 ||
      z_reduction < 0 || z_reduction > 1) {
    stop("Parameters s_M, s_m, u_m, u_f, u_z, z_reduction must be in [0, 1].")
  }
  if (u_f + u_z > 1) {
    stop("u_f + u_z must be <= 1.")
  }
  normalize_freqs(p)
}

mutation_step <- function(p, u_m, u_f, u_z) {
  mutated <- setNames(rep(0, length(haplotypes)), haplotypes)
  for (h in haplotypes) {
    m_allele <- hap_M_allele(h)
    f_allele <- hap_F_allele(h)

    m_targets <- if (m_allele == "M") c(M = 1 - u_m, m = u_m) else c(m = 1)
    f_targets <- if (f_allele == "F") {
      c(F = 1 - u_f - u_z, f = u_f, z = u_z)
    } else {
      setNames(1, f_allele)
    }

    for (m_key in names(m_targets)) {
      for (f_key in names(f_targets)) {
        target <- paste0(m_key, f_key)
        mutated[target] <- mutated[target] + p[h] * m_targets[[m_key]] * f_targets[[f_key]]
      }
    }
  }
  mutated
}

simulate_recursion <- function(cfg, progress_cb = NULL) {
  p <- cfg$p
  s_M <- cfg$s_M
  s_m <- cfg$s_m
  u_m <- cfg$u_m
  u_f <- cfg$u_f
  u_z <- cfg$u_z
  z_reduction <- cfg$z_reduction
  max_generations <- cfg$max_generations
  tol <- cfg$tol

  records <- list()
  records[[1]] <- c(gen = 0, setNames(p, haplotypes))

  for (gen in seq_len(max_generations)) {
    if (!is.null(progress_cb) && (gen == 1 || gen %% 50 == 0 || gen == max_generations)) {
      progress_cb(gen, max_generations)
    }
    next_p <- setNames(rep(0, length(haplotypes)), haplotypes)

    for (i in seq_along(haplotypes)) {
      for (j in i:length(haplotypes)) {
        fg <- if (i == j) p[i]^2 else 2 * p[i] * p[j]
        if (fg == 0) next

        h_i <- haplotypes[i]
        h_j <- haplotypes[j]
        maternal_f <- c(hap_F_allele(h_i), hap_F_allele(h_j))
        maternal_has_F <- any(maternal_f == "F")
        maternal_is_ff <- all(maternal_f == "f")
        maternal_is_Fz <- all(sort(maternal_f) == c("F", "z"))

        if (maternal_has_F) {
          s_m_eff <- if (maternal_is_Fz) s_m * (1 - z_reduction) else s_m
          weights <- ifelse(sapply(haplotypes, function(h) hap_M_allele(h) == "M"), 1, 1 - s_m_eff)
        } else {
          weights <- rep(1, length(haplotypes))
        }

        if (maternal_is_ff) {
          weights <- ifelse(sapply(haplotypes, function(h) hap_M_allele(h) == "M"),
                            (1 - s_M), weights)
        }

        pollen <- p * weights
        pollen_total <- sum(pollen)
        if (pollen_total <= 0) {
          stop("All pollen weights are zero; adjust selection coefficients.")
        }
        pollen <- pollen / pollen_total

        if (i == j) {
          maternal_gametes <- setNames(1, i)
        } else {
          maternal_gametes <- setNames(c(0.5, 0.5), c(i, j))
        }

        for (m_idx in as.integer(names(maternal_gametes))) {
          m_prob <- maternal_gametes[[as.character(m_idx)]]
          for (k in seq_along(haplotypes)) {
            contrib <- fg * m_prob * pollen[k] * 0.5
            next_p[m_idx] <- next_p[m_idx] + contrib
            next_p[k] <- next_p[k] + contrib
          }
        }
      }
    }

    if (u_m > 0 || u_f > 0 || u_z > 0) {
      next_p <- mutation_step(next_p, u_m, u_f, u_z)
    }

    next_p <- normalize_freqs(next_p)
    p <- next_p

    records[[gen + 1]] <- c(gen = gen, setNames(p, haplotypes))

    m_freq <- sum(p[sapply(haplotypes, function(h) hap_M_allele(h) == "M")])
    f_freq <- sum(p[sapply(haplotypes, function(h) hap_F_allele(h) == "F")])
    if (min(m_freq, 1 - m_freq) <= tol && min(f_freq, 1 - f_freq) <= tol) {
      if (!is.null(progress_cb)) {
        progress_cb(gen, max_generations)
      }
      break
    }
  }

  df <- as.data.frame(do.call(rbind, records), stringsAsFactors = FALSE)
  df$gen <- as.integer(df$gen)
  for (h in haplotypes) {
    df[[h]] <- as.numeric(df[[h]])
  }
  df
}

make_long_df <- function(df) {
  do.call(
    rbind,
    lapply(haplotypes, function(h) {
      data.frame(gen = df$gen, haplotype = h, freq = df[[h]], stringsAsFactors = FALSE)
    })
  )
}

ui <- fluidPage(
  titlePanel("Two-Locus Haplotype Recursion (M/F/z)") ,
  sidebarLayout(
    sidebarPanel(
      actionButton("run", "Run simulation"),
      tags$hr(),
      h4("Starting haplotype frequencies"),
      numericInput("p_MF", "MF", value = 1, min = 0, step = 0.01),
      numericInput("p_Mf", "Mf", value = 0, min = 0, step = 0.01),
      numericInput("p_Mz", "Mz", value = 0, min = 0, step = 0.01),
      numericInput("p_mF", "mF", value = 0, min = 0, step = 0.01),
      numericInput("p_mf", "mf", value = 0, min = 0, step = 0.01),
      numericInput("p_mz", "mz", value = 0, min = 0, step = 0.01),
      tags$hr(),
      h4("Selection"),
      numericInput("s_M", "s_M (M pollen on ff)", value = 0.1, min = 0, max = 1, step = 0.01),
      numericInput("s_m", "s_m (m pollen on any F)", value = 0.9, min = 0, max = 1, step = 0.01),
      numericInput("z_reduction", "z reduction (Fz on s_M)", value = 0.5, min = 0, max = 1, step = 0.01),
      tags$hr(),
      h4("Mutation"),
      numericInput("u_m", "u_m (M -> m)", value = 1e-4, min = 0, max = 1, step = 1e-5),
      numericInput("u_f", "u_f (F -> f)", value = 1e-4, min = 0, max = 1, step = 1e-5),
      numericInput("u_z", "u_z (F -> z)", value = 1e-6, min = 0, max = 1, step = 1e-7),
      tags$hr(),
      h4("Run settings"),
      numericInput("max_generations", "max_generations", value = 20000, min = 1, step = 1),
      numericInput("tol", "tol", value = 1e-8, min = 0, step = 1e-8),
      checkboxGroupInput(
        "show_haplotypes",
        "Show haplotypes",
        choices = haplotypes,
        selected = haplotypes
      ),
      checkboxInput("show_mutation", "Show mutation-only trajectories (M, F)", value = TRUE)
    ),
    mainPanel(
      plotOutput("freq_plot", height = 420),
      tableOutput("summary")
    )
  )
)

server <- function(input, output, session) {
  sim_data <- reactiveVal(NULL)

  observeEvent(input$run, {
    p <- c(
      MF = input$p_MF,
      Mf = input$p_Mf,
      Mz = input$p_Mz,
      mF = input$p_mF,
      mf = input$p_mf,
      mz = input$p_mz
    )

    result <- tryCatch({
      p <- validate_params(p, input$s_M, input$s_m, input$u_m, input$u_f,
                           input$u_z, input$z_reduction, input$max_generations, input$tol)
      cfg <- list(
        p = p,
        s_M = input$s_M,
        s_m = input$s_m,
        u_m = input$u_m,
        u_f = input$u_f,
        u_z = input$u_z,
        z_reduction = input$z_reduction,
        max_generations = input$max_generations,
        tol = input$tol
      )
      withProgress(message = "Running simulation", value = 0, {
        simulate_recursion(cfg, progress_cb = function(gen, total) {
          setProgress(gen / total, detail = paste("Generation", gen, "of", total))
        })
      })
    }, error = function(e) e)

    if (inherits(result, "error")) {
      showNotification(result$message, type = "error")
      return(NULL)
    }

    sim_data(result)
  })

  output$freq_plot <- renderPlot({
    df <- sim_data()
    if (is.null(df)) return(NULL)

    show <- input$show_haplotypes
    if (is.null(show) || length(show) == 0) return(NULL)

    long_df <- make_long_df(df)
    long_df <- long_df[long_df$haplotype %in% show, , drop = FALSE]

    p_plot <- ggplot(long_df, aes(x = gen, y = freq, color = haplotype)) +
      geom_line(linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = "Generation", y = "Haplotype frequency") +
      theme_minimal() +
      theme(
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14)
      )

    if (isTRUE(input$show_mutation)) {
      p0 <- c(
        MF = input$p_MF,
        Mf = input$p_Mf,
        Mz = input$p_Mz,
        mF = input$p_mF,
        mf = input$p_mf,
        mz = input$p_mz
      )
      p0 <- normalize_freqs(p0)
      m0 <- sum(p0[c("MF", "Mf", "Mz")])
      f0 <- sum(p0[c("MF", "mF")])
      gens <- df$gen
      m_mut <- m0 * (1 - input$u_m) ^ gens
      f_mut <- f0 * (1 - input$u_f - input$u_z) ^ gens
      mut_df <- data.frame(
        gen = c(gens, gens),
        freq = c(m_mut, f_mut),
        label = c(rep("M (mutation only)", length(gens)),
                  rep("F (mutation only)", length(gens)))
      )
      p_plot <- p_plot +
        geom_line(
          data = mut_df,
          aes(x = gen, y = freq, linetype = label),
          inherit.aes = FALSE,
          color = "grey50",
          linewidth = 1,
          alpha = 0.6
        ) +
        scale_linetype_manual(
          name = "Mutation only",
          values = c("M (mutation only)" = "dotted", "F (mutation only)" = "dotdash")
        )
    }

    p_plot
  })

  output$summary <- renderTable({
    df <- sim_data()
    if (is.null(df)) return(NULL)
    tail(df, 1)
  })
}

shinyApp(ui = ui, server = server)
