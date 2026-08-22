#!/usr/bin/env Rscript

# Deterministic Nature-style renderer for CMM schema-v2 production runs.
# Scientific decisions are made upstream. This script only reads tidy artifacts declared in
# 00_manifest.json and maps their existing columns to publication graphics.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: render_publication_figures.R RUN_DIR MANIFEST OUTPUT_DIR FIGURE_MANIFEST")
}

minimum_r_version <- "4.3.2"
minimum_package_versions <- c(
  jsonlite = "1.8.8",
  ggplot2 = "3.5.2",
  ggrepel = "0.9.6",
  patchwork = "1.2.0",
  svglite = "2.2.1",
  ragg = "1.2.7"
)
required_packages <- names(minimum_package_versions)
if (getRversion() < numeric_version(minimum_r_version)) {
  stop(sprintf(
    paste0(
      "incompatible R version: actual %s; required >= %s. ",
      "Restore the checked-in renv.lock or install a compatible R runtime."
    ),
    as.character(getRversion()),
    minimum_r_version
  ))
}
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(sprintf("missing required R packages: %s", paste(missing_packages, collapse = ", ")))
}
incompatible_packages <- required_packages[vapply(
  required_packages,
  function(package) utils::packageVersion(package) < numeric_version(minimum_package_versions[[package]]),
  logical(1)
)]
if (length(incompatible_packages) > 0L) {
  details <- vapply(
    incompatible_packages,
    function(package) sprintf(
      "%s actual %s; required >= %s",
      package,
      as.character(utils::packageVersion(package)),
      minimum_package_versions[[package]]
    ),
    character(1)
  )
  stop(sprintf(
    "incompatible R package versions: %s. Restore the checked-in renv.lock.",
    paste(details, collapse = "; ")
  ))
}
# A publication render must not silently discard rows or use an invalid scale/device.
options(warn = 2)
zero_reference_flux_tolerance <- 1e-9

run_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
manifest_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
figure_manifest_path <- normalizePath(args[[4]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
if (!identical(manifest$schema_version, 2L)) {
  stop("manifest schema_version must be 2")
}

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0L) right else left
}

plural_noun <- function(count, singular, plural = paste0(singular, "s")) {
  if (count == 1L) singular else plural
}

artifact_entry <- function(role) {
  primary <- manifest$artifacts[[role]]
  if (!is.null(primary)) return(primary)
  supplementary <- manifest$supplementary_artifacts %||% list()
  matches <- Filter(
    function(entry) identical(as.character(entry$role %||% ""), role),
    supplementary
  )
  if (length(matches) == 1L) return(matches[[1L]])
  NULL
}

artifact_status <- function(role) {
  entry <- artifact_entry(role)
  if (is.null(entry)) return("not_declared")
  if (is.character(entry)) return("complete")
  as.character(entry$status %||% "complete")
}

artifact_relative <- function(role) {
  entry <- artifact_entry(role)
  if (is.null(entry)) return(NULL)
  if (is.character(entry)) return(as.character(entry))
  if (is.null(entry$path)) return(NULL)
  as.character(entry$path)
}

artifact_available <- function(role) {
  artifact_status(role) %in% c("complete", "partial") && !is.null(artifact_relative(role))
}

artifact_path <- function(role) {
  relative <- artifact_relative(role)
  if (is.null(relative)) stop(sprintf("artifact role %s has no path", role))
  file.path(run_dir, relative)
}

read_artifact_csv <- function(role) {
  utils::read.csv(
    artifact_path(role),
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NaN")
  )
}

read_artifact_json <- function(role) {
  jsonlite::fromJSON(artifact_path(role), simplifyVector = FALSE)
}

relative_to_run <- function(path) {
  absolute <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(run_dir, "/")
  if (!startsWith(absolute, prefix)) stop(sprintf("output escaped run directory: %s", absolute))
  substring(absolute, nchar(prefix) + 1L)
}

as_flag <- function(value) {
  tolower(trimws(as.character(value))) %in% c("true", "1", "yes", "y")
}

font_requested <- Sys.getenv("CMM_FIGURE_FONT", unset = "Helvetica")
font_family <- font_requested
# The base PDF device provides editable vector text without a Cairo/X11 runtime dependency.
# Restrict its family to the registered PostScript set; Helvetica is Nature-compatible and is
# available in a stock headless R installation on every supported platform.
if (!(font_family %in% names(grDevices::pdfFonts()))) {
  font_family <- "Helvetica"
}
font_resolved <- "not inspected"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  match <- tryCatch(systemfonts::match_fonts(font_requested), error = function(e) NULL)
  if (!is.null(match) && nrow(match) > 0L && nzchar(match$path[[1]])) {
    font_resolved <- basename(match$path[[1]])
  } else {
    font_family <- "Helvetica"
    fallback <- tryCatch(systemfonts::match_fonts(font_family), error = function(e) NULL)
    if (!is.null(fallback) && nrow(fallback) > 0L) {
      font_resolved <- basename(fallback$path[[1]])
    }
  }
}

blue <- "#1769AA"
orange <- "#D55E00"
green <- "#00876C"
purple <- "#8E5EA2"
grey <- "#767676"
light_grey <- "#D9D9D9"
palette_safe <- c(blue, orange, green, purple, "#C58A00", "#4B9CD3", "#444444", "#CC79A7")

theme_nature <- function() {
  ggplot2::theme_classic(base_family = font_family, base_size = 6) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.title = ggplot2::element_text(size = 6.5, colour = "black"),
      axis.text = ggplot2::element_text(size = 6, colour = "black"),
      legend.title = ggplot2::element_text(size = 6),
      legend.text = ggplot2::element_text(size = 5.5),
      legend.key.height = grid::unit(3, "mm"),
      legend.key.width = grid::unit(4, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 6.5, face = "bold", colour = "black"),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "mm")
    )
}

draw_on_device <- function(open_device, plot) {
  open_device()
  device <- grDevices::dev.cur()
  tryCatch(
    print(plot),
    finally = {
      if (grDevices::dev.cur() == device) grDevices::dev.off()
    }
  )
}

save_triplet <- function(plot, figure_id, width_mm, height_mm) {
  paths <- list(
    png = file.path(output_dir, paste0(figure_id, ".png")),
    pdf = file.path(output_dir, paste0(figure_id, ".pdf")),
    svg = file.path(output_dir, paste0(figure_id, ".svg"))
  )
  draw_on_device(
    function() ragg::agg_png(
      filename = paths$png,
      width = width_mm,
      height = height_mm,
      units = "mm",
      res = 300,
      background = "white"
    ),
    plot
  )
  draw_on_device(
    function() grDevices::pdf(
      file = paths$pdf,
      width = width_mm / 25.4,
      height = height_mm / 25.4,
      family = font_family,
      onefile = TRUE,
      useDingbats = FALSE,
      paper = "special"
    ),
    plot
  )
  draw_on_device(
    function() svglite::svglite(
      file = paths$svg,
      width = width_mm / 25.4,
      height = height_mm / 25.4,
      bg = "white"
    ),
    plot
  )
  lapply(paths, relative_to_run)
}

figures <- list()
required_failed <- FALSE

record_skipped <- function(id, label, caption, roles, reason, required = FALSE) {
  figures[[length(figures) + 1L]] <<- list(
    id = id,
    label = label,
    status = if (required) "failed" else "skipped",
    caption = caption,
    alt = caption,
    sources = Filter(Negate(is.null), lapply(roles, artifact_relative)),
    reason = reason
  )
  if (required) required_failed <<- TRUE
}

record_rendered <- function(id, label, caption, alt, roles, width_mm, height_mm, plot) {
  outputs <- save_triplet(plot, id, width_mm, height_mm)
  figures[[length(figures) + 1L]] <<- list(
    id = id,
    label = label,
    status = "rendered",
    caption = caption,
    alt = alt,
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = 300,
    sources = unname(lapply(roles, artifact_relative)),
    outputs = outputs
  )
}

render_panel <- function(id, label, caption, alt, roles, width_mm, height_mm, required, build) {
  unavailable <- roles[!vapply(roles, artifact_available, logical(1))]
  if (length(unavailable) > 0L) {
    reason <- sprintf("Unavailable artifact roles: %s", paste(unavailable, collapse = ", "))
    record_skipped(id, label, caption, roles, reason, required = required)
    return(invisible(NULL))
  }
  tryCatch(
    {
      plot <- build()
      rendered_height_mm <- attr(plot, "cmm_height_mm", exact = TRUE) %||% height_mm
      rendered_caption <- attr(plot, "cmm_caption", exact = TRUE) %||% caption
      rendered_alt <- attr(plot, "cmm_alt", exact = TRUE) %||% alt
      record_rendered(
        id,
        label,
        rendered_caption,
        rendered_alt,
        roles,
        width_mm,
        rendered_height_mm,
        plot
      )
    },
    error = function(error) {
      record_skipped(id, label, caption, roles, conditionMessage(error), required = required)
    }
  )
}

product_label <- as.character(manifest$report$product_label %||% "Target product")

flux_axis_label <- function(prefix) {
  bquote(.(prefix) ~ "(" * mmol ~ gDW^{-1} ~ h^{-1} * ")")
}

growth_axis_label <- function(prefix = "Growth rate") {
  bquote(.(prefix) ~ "(" * h^{-1} * ")")
}

render_panel(
  "fig01_yield_envelope",
  "Figure 1",
  "Feasible growth range across enforced target-product flux (mmol gDW⁻¹ h⁻¹); growth is in h⁻¹ and the annotated molar yield is in mol product per mol substrate.",
  "Production envelope showing minimum and maximum growth against target-product flux.",
  c("theoretical_yield", "production_envelope"),
  89,
  74,
  TRUE,
  function() {
    envelope <- read_artifact_csv("production_envelope")
    yield <- read_artifact_csv("theoretical_yield")
    envelope$product_flux <- as.numeric(envelope$product_flux)
    envelope$growth_min <- as.numeric(envelope$growth_min)
    envelope$growth_max <- as.numeric(envelope$growth_max)
    annotation <- if (nrow(yield) > 0L && is.finite(as.numeric(yield$molar_yield[[1]]))) {
      formatted_yield <- format(
        signif(as.numeric(yield$molar_yield[[1]]), digits = 4),
        trim = TRUE,
        scientific = FALSE
      )
      sprintf("'Yield ceiling:'~%s~mol~mol^{-1}", formatted_yield)
    } else {
      "Yield ceiling not returned"
    }
    ggplot2::ggplot(envelope, ggplot2::aes(x = product_flux)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = growth_min, ymax = growth_max),
        fill = blue,
        alpha = 0.16
      ) +
      ggplot2::geom_line(ggplot2::aes(y = growth_max), colour = blue, linewidth = 0.55) +
      ggplot2::geom_line(
        ggplot2::aes(y = growth_min),
        colour = blue,
        linewidth = 0.4,
        linetype = "22"
      ) +
      ggplot2::annotate(
        "label",
        x = Inf,
        y = Inf,
        label = annotation,
        parse = TRUE,
        hjust = 1.03,
        vjust = 1.25,
        size = 1.9,
        family = font_family,
        label.size = 0.2,
        fill = "white"
      ) +
      ggplot2::labs(
        x = flux_axis_label(sprintf("%s flux", product_label)),
        y = growth_axis_label()
      ) +
      theme_nature()
  }
)

amplification_roles <- c("fseof_tidy", "fvseof_tidy")
for (ranking_role in c(
  "amplification_target_ranking",
  "variability_supported_amplification_targets"
)) {
  if (artifact_available(ranking_role)) {
    amplification_roles <- c(amplification_roles, ranking_role)
  }
}
if (artifact_available("amplification_loop_diagnostic")) {
  amplification_roles <- c(amplification_roles, "amplification_loop_diagnostic")
}

single_knockout_roles <- c("single_knockout_moma", "single_knockout_room", "summary")
if (artifact_available("reproduction_config")) {
  single_knockout_roles <- c(single_knockout_roles, "reproduction_config")
}
if (artifact_available("single_knockout_consensus")) {
  single_knockout_roles <- c(single_knockout_roles, "single_knockout_consensus")
}
if (artifact_available("recommendations")) {
  single_knockout_roles <- c(single_knockout_roles, "recommendations")
}

render_panel(
  "fig02_single_knockout",
  "Figure 2",
  "Single-knockout phenotypes predicted independently by MOMA and ROOM. D1-D5 use the method-specific workflow display rank when declared (current exports retain one representative per blocked-reaction signature); legacy runs fall back to the five highest-product feasible, viability-qualified rows within each method. Every D1-D5 row is a forward-validation candidate, not a recommendation. Grey denotes other feasible screens, blue denotes D1-D5 validation candidates, orange denotes final single-gene support in recommendations.csv, and the black star is the wild-type reference. Infeasible outcomes remain in the source tables.",
  "Faceted scatter plots of predicted product flux against growth rate for MOMA and ROOM single knockouts.",
  single_knockout_roles,
  180,
  96,
  TRUE,
  function() {
    moma <- read_artifact_csv("single_knockout_moma")
    room <- read_artifact_csv("single_knockout_room")
    moma$method <- "MOMA"
    room$method <- "ROOM"
    screen <- rbind(moma, room)
    screen$objective <- suppressWarnings(as.numeric(screen$objective))
    screen$product_flux <- suppressWarnings(as.numeric(screen$product_flux))
    summary <- read_artifact_json("summary")
    wild_type_growth <- as.numeric(summary$wild_type_growth %||% NA_real_)
    viability_fraction <- 0
    if (artifact_available("reproduction_config")) {
      configuration <- read_artifact_json("reproduction_config")
      configured_viability <- suppressWarnings(as.numeric(configuration$viability_fraction %||% NA_real_))
      if (is.finite(configured_viability)) viability_fraction <- configured_viability
    }
    if ("growth_fraction" %in% names(screen)) {
      screen$growth_fraction <- suppressWarnings(as.numeric(screen$growth_fraction))
    } else {
      screen$growth_fraction <- if (is.finite(wild_type_growth) && wild_type_growth > 0) {
        screen$objective / wild_type_growth
      } else {
        rep(NA_real_, nrow(screen))
      }
    }
    declared_display_rank <- if ("display_rank" %in% names(screen)) {
      suppressWarnings(as.numeric(screen$display_rank))
    } else {
      rep(NA_real_, nrow(screen))
    }
    screen$display_rank <- NA_integer_
    for (method_name in c("MOMA", "ROOM")) {
      method_rows <- which(screen$method == method_name)
      declared <- method_rows[
        is.finite(declared_display_rank[method_rows]) &
          declared_display_rank[method_rows] >= 1 &
          declared_display_rank[method_rows] <= 5
      ]
      if (length(declared) > 0L) {
        screen$display_rank[declared] <- as.integer(declared_display_rank[declared])
        next
      }
      eligible <- method_rows[
        screen$status[method_rows] == "optimal" &
          is.finite(screen$objective[method_rows]) &
          is.finite(screen$product_flux[method_rows]) &
          is.finite(screen$growth_fraction[method_rows]) &
          screen$growth_fraction[method_rows] >= viability_fraction
      ]
      eligible <- eligible[order(
        -screen$product_flux[eligible],
        -screen$objective[eligible],
        as.character(screen$target_id[eligible])
      )]
      displayed <- utils::head(eligible, 5L)
      screen$display_rank[displayed] <- seq_along(displayed)
    }
    supported_targets <- character(0)
    if (artifact_available("recommendations")) {
      recommendations <- read_artifact_csv("recommendations")
      supported_targets <- as.character(recommendations$target[
        recommendations$type == "single_gene_knockout" &
          tolower(as.character(recommendations$verdict)) == "support"
      ])
    }
    screen$category <- "Other feasible"
    screen$category[is.finite(screen$display_rank)] <- "Validation candidate (D1-D5)"
    screen$category[screen$target_id %in% supported_targets] <- "Supported recommendation"
    screen$category <- factor(
      screen$category,
      levels = c(
        "Other feasible",
        "Validation candidate (D1-D5)",
        "Supported recommendation",
        "Wild-type reference"
      )
    )
    wild_type <- data.frame(
      objective = rep(wild_type_growth, 2L),
      product_flux = rep(as.numeric(summary$wild_type_product %||% NA_real_), 2L),
      method = c("MOMA", "ROOM"),
      category = factor(
        rep("Wild-type reference", 2L),
        levels = levels(screen$category)
      ),
      stringsAsFactors = FALSE
    )
    wild_type <- wild_type[
      is.finite(wild_type$objective) & is.finite(wild_type$product_flux),
      ,
      drop = FALSE
    ]
    feasible <- screen[
      screen$status == "optimal" & is.finite(screen$objective) & is.finite(screen$product_flux),
      ,
      drop = FALSE
    ]
    if (nrow(feasible) == 0L) {
      return(
        ggplot2::ggplot() +
          ggplot2::annotate(
            "text", x = 0, y = 0, label = "No feasible single-knockout outcome", size = 2.1,
            family = font_family
          ) +
          ggplot2::theme_void(base_family = font_family, base_size = 6)
      )
    }
    feasible <- feasible[order(feasible$category), , drop = FALSE]
    plot <- ggplot2::ggplot(
      feasible,
      ggplot2::aes(x = objective, y = product_flux, colour = category, shape = category)
    ) +
      ggplot2::geom_point(size = 1.35, alpha = 0.82) +
      ggplot2::geom_point(
        data = wild_type,
        size = 2,
        stroke = 0.55,
        alpha = 1
      ) +
      ggplot2::facet_wrap(~method, nrow = 1) +
      ggplot2::scale_colour_manual(
        values = c(
          "Other feasible" = grey,
          "Validation candidate (D1-D5)" = blue,
          "Supported recommendation" = orange,
          "Wild-type reference" = "black"
        ),
        name = NULL
      ) +
      ggplot2::scale_shape_manual(
        values = c(
          "Other feasible" = 16,
          "Validation candidate (D1-D5)" = 17,
          "Supported recommendation" = 18,
          "Wild-type reference" = 8
        ),
        name = NULL
      ) +
      ggplot2::labs(
        x = growth_axis_label("Predicted growth rate"),
        y = flux_axis_label(sprintf("Predicted %s flux", product_label))
      ) +
      theme_nature() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(
        colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
        shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
      )
    labelled <- feasible[is.finite(feasible$display_rank), , drop = FALSE]
    if (nrow(labelled) > 0L) {
      labelled$display_label <- sprintf(
        "D%d %s",
        as.integer(labelled$display_rank),
        as.character(labelled$target_id)
      )
      plot <- plot + ggrepel::geom_text_repel(
        data = labelled,
        ggplot2::aes(label = display_label),
        size = 1.8,
        family = font_family,
        colour = "black",
        box.padding = 0.3,
        point.padding = 0.2,
        min.segment.length = 0,
        segment.size = 0.2,
        max.overlaps = Inf,
        max.iter = 10000,
        force = 1.2,
        show.legend = FALSE,
        seed = 1
      )
    }
    plot
  }
)

render_panel(
  "fig03_strain_design",
  "Figure 3",
  "OptKnock and RobustKnock reaction-knockout designs. Each horizontal interval runs from the guaranteed (left) to maximum (right) product flux at the reported growth optimum; the guaranteed endpoint is drawn on top and its orange fill marks growth coupling.",
  "Point-range chart comparing guaranteed and maximum product for OptKnock and RobustKnock designs.",
  c("optknock", "robustknock"),
  180,
  90,
  TRUE,
  function() {
    opt <- read_artifact_csv("optknock")
    robust <- read_artifact_csv("robustknock")
    opt$method <- rep("OptKnock", nrow(opt))
    robust$method <- rep("RobustKnock", nrow(robust))
    designs <- rbind(opt, robust)
    designs$max_product <- suppressWarnings(as.numeric(designs$max_product))
    designs$guaranteed_product <- suppressWarnings(as.numeric(designs$guaranteed_product))
    designs$growth_coupled <- as_flag(designs$growth_coupled)
    designs <- designs[
      is.finite(designs$max_product) & is.finite(designs$guaranteed_product),
      ,
      drop = FALSE
    ]
    if (nrow(designs) == 0L) {
      return(
        ggplot2::ggplot() +
          ggplot2::annotate(
            "text", x = 0, y = 0, label = "No feasible strain design returned", size = 2.1,
            family = font_family
          ) +
          ggplot2::theme_void(base_family = font_family, base_size = 6)
      )
    }
    designs$display <- factor(designs$knockouts, levels = rev(unique(designs$knockouts)))
    ggplot2::ggplot(designs, ggplot2::aes(y = display)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = guaranteed_product, xend = max_product, yend = display),
        linewidth = 0.55,
        colour = light_grey
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = max_product, shape = "Maximum product"),
        size = 1.8,
        stroke = 0.5,
        colour = blue
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = guaranteed_product,
          fill = growth_coupled,
          shape = "Guaranteed product"
        ),
        size = 1.8,
        stroke = 0.35,
        colour = "black"
      ) +
      ggplot2::facet_wrap(~method, scales = "free_y", nrow = 1) +
      ggplot2::scale_shape_manual(
        values = c("Maximum product" = 4, "Guaranteed product" = 21),
        name = "Endpoint"
      ) +
      ggplot2::scale_fill_manual(
        values = c(`FALSE` = "white", `TRUE` = orange),
        labels = c(`FALSE` = "Not coupled", `TRUE` = "Growth coupled"),
        name = "Guaranteed endpoint"
      ) +
      ggplot2::labs(
        x = flux_axis_label(sprintf("%s flux", product_label)),
        y = "Knockout set"
      ) +
      theme_nature() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(
        shape = ggplot2::guide_legend(order = 1),
        fill = ggplot2::guide_legend(order = 2, override.aes = list(shape = 21))
      )
  }
)

render_panel(
  "fig04_amplification",
  "Figure 4",
  "Independent method-specific amplification trajectories: panel a shows FSEOF ranks 1-10 and panel b FVSEOF ranks 1-10; target intersection is not required. D labels give within-method rank, and a target shared by both methods retains one colour. Loop-flagged targets are marked [loop] with black crosses and placed on a separate free-y diagnostic scale so extreme cycles do not compress eligible trajectories; they are retained in flux-response validation but excluded from support and recommendation eligibility. FVSEOF separates mean flux from forced-minimum magnitude; axes report mmol gDW⁻¹ h⁻¹.",
  "Independent top-ten FSEOF and FVSEOF amplification-target trajectories with direct rank labels.",
  amplification_roles,
  180,
  145,
  TRUE,
  function() {
    fseof <- read_artifact_csv("fseof_tidy")
    fvseof <- read_artifact_csv("fvseof_tidy")
    fseof$enforced_product_flux <- as.numeric(fseof$enforced_product_flux)
    fseof$reaction_flux <- as.numeric(fseof$reaction_flux)
    fvseof$enforced_product_flux <- as.numeric(fvseof$enforced_product_flux)
    fvseof$mean_flux <- as.numeric(fvseof$mean_flux)
    fvseof$forced_min_flux <- as.numeric(fvseof$forced_min_flux)
    fseof$loop_flagged <- if ("loop_artifact_flag" %in% names(fseof)) {
      as_flag(fseof$loop_artifact_flag)
    } else {
      rep(FALSE, nrow(fseof))
    }
    fvseof$loop_flagged <- if ("loop_artifact_flag" %in% names(fvseof)) {
      as_flag(fvseof$loop_artifact_flag)
    } else {
      rep(FALSE, nrow(fvseof))
    }
    if (artifact_available("amplification_loop_diagnostic")) {
      diagnostic <- read_artifact_csv("amplification_loop_diagnostic")
      loop_flagged <- diagnostic$target[as_flag(diagnostic$loop_artifact_flag)]
      fseof$loop_flagged <- fseof$loop_flagged | fseof$target %in% loop_flagged
      fvseof$loop_flagged <- fvseof$loop_flagged | fvseof$target %in% loop_flagged
    }
    select_method_top <- function(data) {
      targets <- unique(as.character(data$target))
      rank_column <- if ("method_rank" %in% names(data)) {
        suppressWarnings(as.numeric(data$method_rank))
      } else {
        rep(NA_real_, nrow(data))
      }
      ranked <- data.frame(
        target = targets,
        rank = vapply(
          targets,
          function(target) {
            values <- rank_column[as.character(data$target) == target]
            values <- values[is.finite(values) & values >= 1]
            if (length(values) > 0L) min(values) else NA_real_
          },
          numeric(1)
        ),
        first = match(targets, as.character(data$target)),
        stringsAsFactors = FALSE
      )
      if (all(!is.finite(ranked$rank))) {
        ranked$rank <- seq_len(nrow(ranked))
      } else {
        unranked <- which(!is.finite(ranked$rank))
        if (length(unranked) > 0L) {
          ranked$rank[unranked] <- max(ranked$rank[is.finite(ranked$rank)]) +
            seq_along(unranked)
        }
      }
      ranked <- ranked[order(ranked$rank, ranked$first, ranked$target), , drop = FALSE]
      ranked <- utils::head(ranked, 10L)
      selected <- data[as.character(data$target) %in% ranked$target, , drop = FALSE]
      selected$method_rank <- ranked$rank[match(as.character(selected$target), ranked$target)]
      selected
    }
    fseof <- select_method_top(fseof)
    fvseof <- select_method_top(fvseof)
    if (nrow(fseof) == 0L || nrow(fvseof) == 0L) {
      stop("both FSEOF and FVSEOF require at least one eligible trajectory")
    }
    union_targets <- sort(unique(c(fseof$target, fvseof$target)))
    if (length(union_targets) > 20L) {
      stop("independent FSEOF/FVSEOF display is limited to 20 unique targets")
    }
    target_colours <- stats::setNames(
      grDevices::hcl.colors(length(union_targets), palette = "viridis"),
      union_targets
    )
    scale_levels <- c("Eligible top-ranked targets", "Loop diagnostic only")
    fseof$scale_group <- factor(
      ifelse(fseof$loop_flagged, scale_levels[[2]], scale_levels[[1]]),
      levels = scale_levels
    )
    fvseof$scale_group <- factor(
      ifelse(fvseof$loop_flagged, scale_levels[[2]], scale_levels[[1]]),
      levels = scale_levels
    )
    endpoint_rows <- function(data, value_column) {
      indices <- vapply(
        split(seq_len(nrow(data)), as.character(data$target)),
        function(index) index[which.max(data$enforced_product_flux[index])],
        integer(1)
      )
      endpoints <- data[unname(indices), , drop = FALSE]
      endpoints$label_y <- endpoints[[value_column]]
      endpoints$display_label <- sprintf(
        "D%d %s%s",
        as.integer(endpoints$method_rank),
        as.character(endpoints$target),
        ifelse(endpoints$loop_flagged, " [loop]", "")
      )
      endpoints
    }
    fseof_endpoints <- endpoint_rows(fseof, "reaction_flux")
    panel_a <- ggplot2::ggplot(
      fseof,
      ggplot2::aes(
        x = enforced_product_flux,
        y = reaction_flux,
        colour = target,
        group = target
      )
    ) +
      ggplot2::geom_line(ggplot2::aes(alpha = loop_flagged), linewidth = 0.5) +
      ggplot2::geom_point(ggplot2::aes(alpha = loop_flagged), size = 0.8) +
      ggplot2::geom_point(
        data = fseof_endpoints[fseof_endpoints$loop_flagged, , drop = FALSE],
        ggplot2::aes(x = enforced_product_flux, y = label_y),
        inherit.aes = FALSE,
        shape = 4,
        size = 1.5,
        stroke = 0.45,
        colour = "black"
      ) +
      ggrepel::geom_text_repel(
        data = fseof_endpoints,
        ggplot2::aes(
          x = enforced_product_flux,
          y = label_y,
          label = display_label
        ),
        inherit.aes = FALSE,
        size = 1.85,
        colour = "#222222",
        family = font_family,
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0,
        segment.size = 0.18,
        max.overlaps = Inf,
        max.iter = 10000,
        seed = 2,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = target_colours,
        limits = union_targets,
        drop = FALSE,
        guide = "none"
      ) +
      ggplot2::scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.38), guide = "none") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.2))) +
      ggplot2::facet_wrap(~scale_group, ncol = 1, scales = "free_y", drop = TRUE) +
      ggplot2::labs(
        title = "FSEOF: independent top ten",
        x = flux_axis_label(sprintf("Enforced %s flux", product_label)),
        y = flux_axis_label("Target-reaction flux")
      ) +
      theme_nature() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 7, face = "bold", hjust = 0)
      )
    fv_mean <- data.frame(
      target = fvseof$target,
      enforced_product_flux = fvseof$enforced_product_flux,
      flux = fvseof$mean_flux,
      metric = "Mean flux",
      loop_flagged = fvseof$loop_flagged,
      scale_group = fvseof$scale_group,
      stringsAsFactors = FALSE
    )
    fv_forced <- data.frame(
      target = fvseof$target,
      enforced_product_flux = fvseof$enforced_product_flux,
      flux = fvseof$forced_min_flux,
      metric = "Forced minimum",
      method_rank = fvseof$method_rank,
      loop_flagged = fvseof$loop_flagged,
      scale_group = fvseof$scale_group,
      stringsAsFactors = FALSE
    )
    fv_mean$method_rank <- fvseof$method_rank
    fv_long <- rbind(fv_mean, fv_forced)
    fvseof_endpoints <- endpoint_rows(fvseof, "mean_flux")
    panel_b <- ggplot2::ggplot(
      fv_long,
      ggplot2::aes(
        x = enforced_product_flux,
        y = flux,
        colour = target,
        linetype = metric,
        group = interaction(target, metric)
      )
    ) +
      ggplot2::geom_line(
        ggplot2::aes(alpha = loop_flagged),
        linewidth = 0.5,
        show.legend = c(colour = FALSE, linetype = TRUE)
      ) +
      ggplot2::geom_point(
        data = fvseof_endpoints[fvseof_endpoints$loop_flagged, , drop = FALSE],
        ggplot2::aes(x = enforced_product_flux, y = label_y),
        inherit.aes = FALSE,
        shape = 4,
        size = 1.5,
        stroke = 0.45,
        colour = "black"
      ) +
      ggrepel::geom_text_repel(
        data = fvseof_endpoints,
        ggplot2::aes(
          x = enforced_product_flux,
          y = label_y,
          label = display_label
        ),
        inherit.aes = FALSE,
        size = 1.85,
        colour = "#222222",
        family = font_family,
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0,
        segment.size = 0.18,
        max.overlaps = Inf,
        max.iter = 10000,
        seed = 3,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = target_colours,
        limits = union_targets,
        drop = FALSE,
        guide = "none"
      ) +
      ggplot2::scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.38), guide = "none") +
      ggplot2::scale_linetype_manual(values = c("Mean flux" = "solid", "Forced minimum" = "22")) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.2))) +
      ggplot2::facet_wrap(~scale_group, ncol = 1, scales = "free_y", drop = TRUE) +
      ggplot2::labs(
        title = "FVSEOF: independent top ten",
        x = flux_axis_label(sprintf("Enforced %s flux", product_label)),
        y = flux_axis_label("Target-reaction flux"),
        linetype = "FVSEOF metric"
      ) +
      theme_nature() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 7, face = "bold", hjust = 0),
        legend.position = "bottom"
      )
    combined <- panel_a + panel_b +
      patchwork::plot_annotation(tag_levels = "a")
    combined & ggplot2::theme(
      plot.tag = ggplot2::element_text(family = font_family, size = 7, face = "bold")
    )
  }
)

render_panel(
  "fig05_flux_response",
  "Figure 5",
  "Standard target-to-product flux-response curves for every completed compatible candidate analysis. Every facet uses enforced candidate-reaction flux (target_flux) on the x-axis and optimized target-product flux (response_flux) on the y-axis. Biomass flux is a recorded secondary value under the configured minimum-growth constraint, not a plot axis. The index candidate scope separates amplification hypotheses from wild-type pre-deletion titrations of single-reaction knockout candidates; multi-reaction signatures remain explicit skipped or unavailable index rows. Legacy schema-v2 product-to-growth scans are retained in their CSV but are not relabelled as product-response panels.",
  "Multi-row enforced-candidate-reaction-flux versus target-product-flux panels for all completed compatible amplification and knockout-derived candidate scans.",
  c(
    "flux_response_tidy",
    if (artifact_available("flux_response_validation_index")) {
      "flux_response_validation_index"
    }
  ),
  180,
  150,
  FALSE,
  function() {
    response <- read_artifact_csv("flux_response_tidy")
    response$target_flux <- as.numeric(response$target_flux)
    response$response_flux <- as.numeric(response$response_flux)
    response$biomass_flux <- as.numeric(response$biomass_flux)
    if (any(!(response$background %in% c("wild_type", "gene_knockout")))) {
      stop("flux-response table contains an unknown model background")
    }
    response <- response[
      response$status %in% "optimal" &
        is.finite(response$target_flux) &
        is.finite(response$response_flux),
      ,
      drop = FALSE
    ]
    if (nrow(response) == 0L) stop("no optimal flux-response point was available")
    response$candidate_scope <- if ("candidate_scope" %in% names(response)) {
      as.character(response$candidate_scope)
    } else {
      rep(NA_character_, nrow(response))
    }
    response_index <- data.frame()
    if (artifact_available("flux_response_validation_index")) {
      response_index <- read_artifact_csv("flux_response_validation_index")
      if ("candidate_scope" %in% names(response_index)) {
        indexed <- response_index[
          response_index$status == "complete" &
            response_index$candidate_scope %in% c(
              "all_report_selected_candidates",
              "all_display_ranked_candidates"
            ),
          c("target", "candidate_scope"),
          drop = FALSE
        ]
        ambiguous <- names(which(vapply(
          split(indexed$candidate_scope, as.character(indexed$target)),
          function(values) length(unique(values)) > 1L,
          logical(1)
        )))
        if (length(ambiguous) > 0L) {
          stop(sprintf(
            "flux-response index assigns multiple candidate scopes to: %s",
            paste(ambiguous, collapse = ", ")
          ))
        }
        scope_map <- stats::setNames(
          as.character(indexed$candidate_scope),
          as.character(indexed$target)
        )
        indexed_scope <- unname(scope_map[as.character(response$target)])
        conflicts <- !is.na(response$candidate_scope) &
          nzchar(response$candidate_scope) &
          !is.na(indexed_scope) &
          response$candidate_scope != indexed_scope
        if (any(conflicts)) {
          stop(sprintf(
            "flux-response tidy/index candidate scopes disagree for: %s",
            paste(unique(as.character(response$target[conflicts])), collapse = ", ")
          ))
        }
        missing_scope <- is.na(response$candidate_scope) | !nzchar(response$candidate_scope)
        response$candidate_scope[missing_scope] <- indexed_scope[missing_scope]
      }
    }
    exploratory_targets <- character()
    noncomplete_knockout_targets <- character()
    if (nrow(response_index) > 0L && "candidate_scope" %in% names(response_index)) {
      knockout_index <- response_index[
        response_index$candidate_scope == "all_display_ranked_candidates",
        ,
        drop = FALSE
      ]
      noncomplete_knockout_targets <- unique(as.character(
        knockout_index$target[knockout_index$status != "complete"]
      ))
      if ("scan_reference_flux" %in% names(knockout_index)) {
        reference_flux <- suppressWarnings(as.numeric(knockout_index$scan_reference_flux))
        exploratory_targets <- unique(as.character(
          knockout_index$target[
            knockout_index$status == "complete" &
              is.finite(reference_flux) &
              abs(reference_flux) <= zero_reference_flux_tolerance
          ]
        ))
      }
    }
    response$panel_class <- ifelse(
      response$candidate_scope == "all_report_selected_candidates",
      "amplification",
      ifelse(
        response$candidate_scope == "all_display_ranked_candidates",
        "knockout",
        NA_character_
      )
    )
    summary <- if (artifact_available("summary")) read_artifact_json("summary") else list()
    product_reaction <- as.character(summary$product %||% "")
    unclassified <- is.na(response$panel_class)
    product_response <- if (nzchar(product_reaction)) {
      as.character(response$response_reaction) == product_reaction
    } else {
      response$background == "wild_type"
    }
    response$panel_class[
      unclassified & product_response & response$background == "wild_type"
    ] <- "amplification"
    response$panel_class[
      unclassified & product_response & response$background == "gene_knockout"
    ] <- "knockout"
    standard_rows <- !is.na(response$panel_class)
    incompatible_rows <- standard_rows & !product_response
    response$panel_class[incompatible_rows] <- NA_character_
    omitted_targets <- unique(as.character(response$target[is.na(response$panel_class)]))
    plotted <- response[!is.na(response$panel_class), , drop = FALSE]
    if (nrow(plotted) == 0L) {
      stop("no standard target-to-product flux-response point was available")
    }
    amplification <- plotted[plotted$panel_class == "amplification", , drop = FALSE]
    knockout <- plotted[plotted$panel_class == "knockout", , drop = FALSE]
    knockout_is_current_scope <- nrow(knockout) > 0L && all(
      !is.na(knockout$candidate_scope) &
        knockout$candidate_scope == "all_display_ranked_candidates"
    )
    knockout_panel_title <- if (knockout_is_current_scope) {
      "Knockout-derived candidate reactions (wild type, pre-deletion)"
    } else {
      "Knockout-derived candidate reaction scans (recorded backgrounds)"
    }
    knockout_context <- if (knockout_is_current_scope) {
      paste0(
        "Knockout-derived panels are wild-type pre-deletion single-reaction titrations; ",
        "their facet titles pair the representative gene with its scanned reaction."
      )
    } else {
      paste0(
        "Current scoped knockout-derived panels are wild-type pre-deletion single-reaction ",
        "titrations. Legacy unscoped panels retain their recorded model background and are ",
        "not reinterpreted as current scoped knockout evidence."
      )
    }
    amplification_targets <- unique(as.character(amplification$target))
    knockout_targets <- unique(as.character(knockout$target))
    amplification$facet_label <- as.character(amplification$target)
    knockout$facet_label <- ifelse(
      as.character(knockout$target) == as.character(knockout$scan_reaction),
      as.character(knockout$target),
      sprintf(
        "%s (%s)",
        as.character(knockout$target),
        as.character(knockout$scan_reaction)
      )
    )
    amplification$target <- factor(amplification$target, levels = amplification_targets)
    knockout$target <- factor(knockout$target, levels = knockout_targets)
    amplification$facet_label <- factor(
      amplification$facet_label,
      levels = unique(as.character(amplification$facet_label))
    )
    knockout$facet_label <- factor(
      knockout$facet_label,
      levels = unique(as.character(knockout$facet_label))
    )
    amplification <- amplification[
      order(amplification$target, amplification$target_flux),
      ,
      drop = FALSE
    ]
    knockout <- knockout[order(knockout$target, knockout$target_flux), , drop = FALSE]
    finite_limits <- function(values) {
      limits <- range(values[is.finite(values)])
      if (diff(limits) == 0) {
        padding <- max(abs(limits[[1]]) * 0.05, 0.05)
        limits <- limits + c(-padding, padding)
      }
      limits
    }
    facet_flux_breaks <- function(limits) {
      if (length(limits) != 2L || any(!is.finite(limits))) return(numeric())
      pretty(limits, n = 3)
    }
    shared_response_limits <- finite_limits(plotted$response_flux)
    panels <- list()
    if (nrow(amplification) > 0L) {
      amplification_columns <- min(4L, length(amplification_targets))
      panels[[length(panels) + 1L]] <- ggplot2::ggplot(
        amplification,
        ggplot2::aes(x = target_flux, y = response_flux, group = target)
      ) +
        ggplot2::geom_path(linewidth = 0.55, colour = blue) +
        ggplot2::geom_point(size = 0.85, colour = blue) +
        ggplot2::facet_wrap(
          ~facet_label,
          scales = "free_x",
          ncol = amplification_columns
        ) +
        ggplot2::scale_x_continuous(breaks = facet_flux_breaks) +
        ggplot2::scale_y_continuous(
          limits = shared_response_limits,
          n.breaks = 4,
          expand = ggplot2::expansion(mult = 0.04)
        ) +
        ggplot2::labs(
          title = "Amplification candidates (wild type)",
          x = flux_axis_label("Enforced candidate-reaction flux"),
          y = flux_axis_label("Target-product flux")
        ) +
        theme_nature() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 7, face = "bold", hjust = 0),
          panel.spacing.x = grid::unit(4, "mm")
        )
    }
    if (nrow(knockout) > 0L) {
      knockout_columns <- min(4L, length(knockout_targets))
      panels[[length(panels) + 1L]] <- ggplot2::ggplot(
        knockout,
        ggplot2::aes(x = target_flux, y = response_flux, group = target)
      ) +
        ggplot2::geom_path(linewidth = 0.55, colour = orange) +
        ggplot2::geom_point(size = 0.85, colour = orange) +
        ggplot2::facet_wrap(
          ~facet_label,
          scales = "free_x",
          ncol = knockout_columns
        ) +
        ggplot2::scale_x_continuous(breaks = facet_flux_breaks) +
        ggplot2::scale_y_continuous(
          limits = shared_response_limits,
          n.breaks = 4,
          expand = ggplot2::expansion(mult = 0.04)
        ) +
        ggplot2::labs(
          title = knockout_panel_title,
          x = flux_axis_label("Enforced candidate-reaction flux"),
          y = flux_axis_label("Target-product flux")
        ) +
        theme_nature() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 7, face = "bold", hjust = 0),
          panel.spacing.x = grid::unit(4, "mm")
        )
    }
    if (length(panels) == 0L) stop("no supported flux-response background was available")
    combined <- patchwork::wrap_plots(panels, ncol = 1) +
      patchwork::plot_annotation(tag_levels = "a") &
      ggplot2::theme(
        plot.tag = ggplot2::element_text(family = font_family, size = 7, face = "bold")
      )
    total_rows <-
      if (length(amplification_targets) > 0L) {
        ceiling(length(amplification_targets) / 4)
      } else {
        0L
      }
    total_rows <- total_rows +
      if (length(knockout_targets) > 0L) ceiling(length(knockout_targets) / 4) else 0L
    attr(combined, "cmm_height_mm") <- min(230, max(150, 27 * total_rows + 35))
    amplification_count <- length(amplification_targets)
    knockout_count <- length(knockout_targets)
    omitted_note <- if (length(omitted_targets) > 0L) {
      sprintf(
        paste0(
          " %d legacy or incompatible %s (%s) were retained in the source CSV but not ",
          "relabeled as target-to-product response panels."
        ),
        length(omitted_targets),
        plural_noun(length(omitted_targets), "target"),
        paste(omitted_targets, collapse = ", ")
      )
    } else {
      ""
    }
    exploratory_note <- if (length(exploratory_targets) > 0L) {
      sprintf(
        paste0(
          " %d zero-reference knockout-derived %s (%s) use the full feasible reaction ",
          "domain and are exploratory; they cannot causally support deletion."
        ),
        length(exploratory_targets),
        plural_noun(length(exploratory_targets), "scan"),
        paste(exploratory_targets, collapse = ", ")
      )
    } else {
      ""
    }
    noncomplete_note <- if (length(noncomplete_knockout_targets) > 0L) {
      sprintf(
        " %d knockout-derived %s had no panel and remain explicit non-complete index rows.",
        length(noncomplete_knockout_targets),
        plural_noun(length(noncomplete_knockout_targets), "candidate")
      )
    } else {
      ""
    }
    attr(combined, "cmm_caption") <- sprintf(
      paste0(
        "Standard target-to-product flux-response curves for every completed compatible ",
        "candidate analysis: %d total, comprising %d amplification %s and %d knockout-derived ",
        "%s, arranged in multi-row grids of at most four facets per row. Every facet maps ",
        "enforced candidate-reaction flux (target_flux) on the x-axis to optimized target-product ",
        "flux (response_flux) on the y-axis. Product-flux limits are shared; candidate-",
        "reaction x ranges vary by facet. Biomass flux records the secondary state under the ",
        "configured minimum-growth constraint and is not a plot axis. %s Multi-reaction signatures remain ",
        "explicit skipped or unavailable index rows. Fluxes are in mmol gDW⁻¹ h⁻¹. Candidate ",
        "inclusion is independent of recommendation status. Optimal points are connected in ",
        "target_flux scan order; all rows and failed execution reasons remain in ",
        "flux_response_tidy.csv and flux_response_index.csv.%s%s%s ",
        "Gene IDs sharing one blocked-reaction signature share a simulation facet; the index ",
        "and report table list every candidate alias."
      ),
      amplification_count + knockout_count,
      amplification_count,
      plural_noun(amplification_count, "target"),
      knockout_count,
      plural_noun(knockout_count, "candidate-reaction titration"),
      knockout_context,
      omitted_note,
      exploratory_note,
      noncomplete_note
    )
    combined
  }
)

render_panel(
  "fig06_sampling_shift",
  "Figure 6",
  "Paired feasible-flux distributions for the target-product exchange (top; mmol gDW⁻¹ h⁻¹) and biomass reaction (bottom; h⁻¹) before and after every completed display-ranked single-gene-knockout sampling analysis. Candidate inclusion is independent of recommendation status. The main figure intentionally omits other reactions, which remain available in sampling_tidy.csv. Samples are feasible states, not biological replicates.",
  "Multi-row product-flux and biomass-growth distributions for all completed display-ranked single-gene-knockout sampling analyses.",
  c(
    "sampling_tidy",
    "summary",
    if (artifact_available("single_knockout_sampling_validation_index")) {
      "single_knockout_sampling_validation_index"
    }
  ),
  180,
  118,
  FALSE,
  function() {
    samples <- read_artifact_csv("sampling_tidy")
    summary <- read_artifact_json("summary")
    product_reaction <- as.character(summary$product %||% "")
    biomass_reaction <- as.character(summary$biomass %||% "")
    if (!nzchar(product_reaction) || !nzchar(biomass_reaction)) {
      stop("summary must identify product and biomass reactions for the sampling figure")
    }
    samples$flux <- as.numeric(samples$flux)
    samples <- samples[
      is.finite(samples$flux) &
        samples$reaction_id %in% c(product_reaction, biomass_reaction),
      ,
      drop = FALSE
    ]
    if (nrow(samples) == 0L) stop("sampling table contained no finite flux value")
    samples$condition <- factor(
      samples$condition,
      levels = c("wild_type", "knockout"),
      labels = c("Wild type", "Knockout")
    )
    target_levels <- unique(as.character(samples$target))
    samples$target <- factor(samples$target, levels = target_levels)
    facet_columns <- min(5L, length(target_levels))
    distribution_panel <- function(data, title, y_label) {
      ggplot2::ggplot(
        data,
        ggplot2::aes(x = condition, y = flux, fill = condition)
      ) +
        ggplot2::geom_violin(
          scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.72
        ) +
        ggplot2::geom_boxplot(
          width = 0.16, outlier.shape = NA, linewidth = 0.3, fill = "white"
        ) +
        ggplot2::facet_wrap(~target, scales = "free_y", ncol = facet_columns) +
        ggplot2::scale_fill_manual(
          values = c("Wild type" = grey, "Knockout" = orange),
          drop = FALSE
        ) +
        ggplot2::labs(title = title, x = NULL, y = y_label, fill = "Condition") +
        theme_nature() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 7, face = "bold", hjust = 0),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
    product_samples <- samples[samples$reaction_id == product_reaction, , drop = FALSE]
    biomass_samples <- samples[samples$reaction_id == biomass_reaction, , drop = FALSE]
    if (nrow(product_samples) == 0L || nrow(biomass_samples) == 0L) {
      stop("sampling figure requires both product and biomass distributions")
    }
    product_panel <- distribution_panel(
      product_samples,
      "Product exchange",
      flux_axis_label("Product flux")
    )
    biomass_panel <- distribution_panel(
      biomass_samples,
      "Biomass reaction",
      growth_axis_label()
    )
    combined <- product_panel + biomass_panel +
      patchwork::plot_layout(ncol = 1, guides = "collect") +
      patchwork::plot_annotation(tag_levels = "a")
    combined <- combined & ggplot2::theme(
      legend.position = "bottom",
      plot.tag = ggplot2::element_text(family = font_family, size = 7, face = "bold")
    )
    target_rows <- ceiling(length(target_levels) / 5)
    attr(combined, "cmm_height_mm") <- min(230, max(118, 58 * target_rows + 58))
    attr(combined, "cmm_caption") <- sprintf(
      paste0(
        "Paired feasible-flux distributions for the product exchange (top; mmol gDW⁻¹ h⁻¹) ",
        "and biomass reaction (bottom; h⁻¹) before and after all %d completed display-ranked ",
        "single-gene-knockout sampling analyses, arranged in grids of at most five targets per ",
        "row. Candidate inclusion is independent of recommendation status. The main figure ",
        "intentionally omits other reactions, which remain in sampling_tidy.csv; samples are ",
        "feasible states, not biological replicates. ",
        "Gene IDs sharing one blocked-reaction signature share a simulation facet; the index ",
        "and report table list every candidate alias."
      ),
      length(target_levels)
    )
    combined
  }
)

renderer_versions <- as.list(vapply(
  required_packages,
  function(package) as.character(utils::packageVersion(package)),
  character(1)
))

figure_manifest <- list(
  schema_version = 2L,
  renderer = list(
    engine = "R/ggplot2",
    r = R.version.string,
    minimum_r = minimum_r_version,
    packages = renderer_versions,
    minimum_packages = as.list(minimum_package_versions),
    script_sha256 = Sys.getenv("CMM_RENDERER_SHA256", unset = "unknown"),
    font_requested = font_requested,
    font_family = font_family,
    font_resolved = font_resolved,
    specification = list(
      widths_mm = c(89L, 180L),
      raster_dpi = 300L,
      vector_formats = c("pdf", "svg"),
      colour_space = "RGB"
    )
  ),
  figures = figures
)

jsonlite::write_json(
  figure_manifest,
  path = figure_manifest_path,
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  digits = NA
)

if (required_failed) {
  stop("one or more required publication figures failed; see figures/figure_manifest.json")
}
