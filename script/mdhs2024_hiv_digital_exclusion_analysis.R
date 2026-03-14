#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 140,
  survey.lonely.psu = "adjust",
  dplyr.summarise.inform = FALSE,
  repos = c(CRAN = "https://cloud.r-project.org")
)

required_pkgs <- c(
  "haven", "labelled", "dplyr", "tidyr", "stringr", "purrr", "forcats",
  "readr", "tibble", "ggplot2", "scales", "glue", "survey", "srvyr",
  "broom", "gt", "htmltools", "fs", "janitor", "rlang", "cli", "ggrepel"
)

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs, dependencies = TRUE)
}
install_if_missing(required_pkgs)
invisible(lapply(required_pkgs, library, character.only = TRUE))

ROOT_DIR <- normalizePath(".", winslash = "/", mustWork = TRUE)
DIR_OUT  <- file.path(ROOT_DIR, "outputs_h7_digital_exclusion_hiv_FINAL")
DIR_FIG  <- file.path(DIR_OUT, "figures")
DIR_TAB  <- file.path(DIR_OUT, "tables")
DIR_LOG  <- file.path(DIR_OUT, "logs")
DIR_RDS  <- file.path(DIR_OUT, "rds")
fs::dir_create(c(DIR_OUT, DIR_FIG, DIR_TAB, DIR_LOG, DIR_RDS), recurse = TRUE)

log_file <- file.path(DIR_LOG, paste0("analysis_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
log_con  <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

assert_that <- function(condition, msg) {
  if (!isTRUE(condition)) stop(msg, call. = FALSE)
}

fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ formatC(p, format = "e", digits = 2),
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}
fmt_num <- function(x, digits = 2) ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
fmt_pct <- function(x, digits = 1) ifelse(is.na(x), NA_character_, paste0(formatC(100 * x, format = "f", digits = digits), "%"))
or_ci_string <- function(est, lo, hi, digits = 2) paste0(fmt_num(est, digits), " (", fmt_num(lo, digits), ", ", fmt_num(hi, digits), ")")
pr_ci_string <- or_ci_string

pick_strata_var <- function(dat) {
  cand <- c("v022", "v023")
  existing <- cand[cand %in% names(dat)]
  assert_that(length(existing) >= 1, "Neither V022 nor V023 was found in the IR file.")
  existing[[1]]
}

relevel_if_present <- function(x, preferred) {
  if (!is.factor(x)) return(x)
  existing <- levels(x)
  ord <- c(preferred[preferred %in% existing], setdiff(existing, preferred))
  factor(x, levels = ord)
}

weighted_mean_wide <- function(dat, var, by_var) {
  dat |>
    dplyr::group_by(.data[[by_var]]) |>
    dplyr::summarise(value = stats::weighted.mean(.data[[var]], wt, na.rm = TRUE), .groups = "drop") |>
    dplyr::rename(group = !!rlang::sym(by_var))
}

weighted_dist_wide <- function(dat, var, by_var) {
  overall <- dat |>
    dplyr::count(level = .data[[var]], wt = wt, name = "weighted_n") |>
    dplyr::mutate(overall = weighted_n / sum(weighted_n, na.rm = TRUE)) |>
    dplyr::select(level, overall)
  
  by_group <- dat |>
    dplyr::count(group = .data[[by_var]], level = .data[[var]], wt = wt, name = "weighted_n") |>
    dplyr::group_by(group) |>
    dplyr::mutate(prop = weighted_n / sum(weighted_n, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::select(level, group, prop) |>
    tidyr::pivot_wider(names_from = group, values_from = prop)
  
  overall |>
    dplyr::left_join(by_group, by = "level")
}

make_design <- function(dat, strata_var) {
  survey::svydesign(
    ids = ~v021,
    strata = stats::as.formula(paste0("~", strata_var)),
    weights = ~wt,
    nest = TRUE,
    data = dat
  )
}

weighted_binary_prop <- function(design_obj, var_name) {
  out <- survey::svymean(stats::as.formula(paste0("~", var_name)), design = design_obj, na.rm = TRUE)
  ci <- suppressWarnings(stats::confint(out))
  tibble::tibble(estimate = as.numeric(out)[1], ci_low = ci[1, 1], ci_high = ci[1, 2])
}

coef_table <- function(model, exponentiate = FALSE) {
  broom::tidy(model, conf.int = TRUE, exponentiate = exponentiate) |>
    dplyr::mutate(
      estimate_fmt = fmt_num(estimate, 3),
      conf_low_fmt = fmt_num(conf.low, 3),
      conf_high_fmt = fmt_num(conf.high, 3),
      p_value_fmt = fmt_p(p.value)
    )
}

compute_evalue_rr <- function(rr, lo) {
  rr <- as.numeric(rr); lo <- as.numeric(lo)
  if (is.na(rr) || rr <= 0) return(tibble::tibble(evalue_point = NA_real_, evalue_ci = NA_real_))
  e_point <- if (rr >= 1) rr + sqrt(rr * (rr - 1)) else (1 / rr) + sqrt((1 / rr) * ((1 / rr) - 1))
  e_ci <- if (is.na(lo)) {
    NA_real_
  } else if (lo <= 1 && rr >= 1) {
    1
  } else if (lo >= 1) {
    lo + sqrt(lo * (lo - 1))
  } else {
    invlo <- 1 / lo
    invlo + sqrt(invlo * (invlo - 1))
  }
  tibble::tibble(evalue_point = e_point, evalue_ci = e_ci)
}

extract_effect <- function(model, term = "internet_any", scale = c("OR", "PR")) {
  scale <- match.arg(scale)
  td <- broom::tidy(model, conf.int = TRUE, exponentiate = TRUE)
  row <- td |>
    dplyr::filter(.data$term == !!term) |>
    dplyr::slice(1)
  assert_that(nrow(row) == 1, paste0("Term not found in model: ", term))
  tibble::tibble(
    estimate = row$estimate,
    conf_low = row$conf.low,
    conf_high = row$conf.high,
    p_value = row$p.value,
    effect_scale = scale
  )
}

predict_svyglm_safe <- function(model, newdata, type = c("link", "response")) {
  type <- match.arg(type)
  pred <- stats::predict(model, newdata = newdata, type = type, se.fit = TRUE)
  
  if (is.list(pred) && !is.null(pred$fit)) {
    fit <- as.numeric(pred$fit)
    se  <- as.numeric(pred$se.fit)
  } else {
    fit <- as.numeric(pred)
    pred_var <- attr(pred, "var")
    if (is.null(pred_var)) {
      se <- rep(NA_real_, length(fit))
    } else if (is.matrix(pred_var)) {
      se <- sqrt(pmax(diag(pred_var), 0))
    } else {
      se <- sqrt(pmax(as.numeric(pred_var), 0))
    }
  }
  
  out <- dplyr::bind_cols(newdata, tibble::tibble(fit = fit, se = se))
  if (type == "link") {
    out |>
      dplyr::mutate(
        prob = plogis(fit),
        prob_lo = plogis(fit - 1.96 * se),
        prob_hi = plogis(fit + 1.96 * se)
      )
  } else {
    out |>
      dplyr::mutate(
        prob = fit,
        prob_lo = pmax(0, fit - 1.96 * se),
        prob_hi = pmin(1, fit + 1.96 * se)
      )
  }
}

collapse_term_labels <- function(term) {
  dplyr::case_when(
    term == "(Intercept)" ~ "Intercept",
    term == "internet_any" ~ "Any internet use vs none",
    term == "age_years" ~ "Age (per 1-year increase)",
    stringr::str_starts(term, "education") ~ stringr::str_replace(term, "education", "Educational attainment: "),
    stringr::str_starts(term, "wealth_quintile") ~ stringr::str_replace(term, "wealth_quintile", "Household wealth quintile: "),
    stringr::str_starts(term, "residence") ~ stringr::str_replace(term, "residence", "Residence: "),
    stringr::str_starts(term, "marital_status") ~ stringr::str_replace(term, "marital_status", "Marital status: "),
    stringr::str_starts(term, "region") ~ stringr::str_replace(term, "region", "Region: "),
    TRUE ~ term
  )
}

safe_regterm_p <- function(model, term_formula) {
  out <- tryCatch(survey::regTermTest(model, term_formula), error = function(e) NULL)
  if (is.null(out)) return(NA_real_)
  if (!is.null(out$p)) return(out$p)
  NA_real_
}

save_table_dual <- function(df, stem, title, subtitle = NULL, note = NULL) {
  csv_path  <- file.path(DIR_TAB, paste0(stem, ".csv"))
  html_path <- file.path(DIR_TAB, paste0(stem, ".html"))
  readr::write_csv(df, csv_path, na = "")
  gt_tbl <- gt::gt(df) |>
    gt::tab_header(
      title = gt::md(title),
      subtitle = if (!is.null(subtitle)) gt::md(subtitle) else NULL
    ) |>
    gt::opt_row_striping() |>
    gt::fmt_missing(everything(), missing_text = "") |>
    gt::tab_options(
      table.width = gt::pct(100),
      table.font.names = c("Arial", "Helvetica", "sans-serif"),
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(16),
      heading.subtitle.font.size = gt::px(12),
      data_row.padding = gt::px(4),
      column_labels.font.weight = "bold",
      table_body.hlines.color = "#E5E7EB",
      table.border.top.width = gt::px(2),
      table.border.bottom.width = gt::px(2),
      source_notes.font.size = gt::px(10)
    )
  if (!is.null(note)) gt_tbl <- gt_tbl |> gt::tab_source_note(gt::md(note))
  gt::gtsave(gt_tbl, html_path)
  invisible(list(csv = csv_path, html = html_path))
}

save_plot650 <- function(plot_obj, filename, width = 10.5, height = 7.2, bg = "white") {
  ggplot2::ggsave(
    filename = file.path(DIR_FIG, filename),
    plot = plot_obj,
    width = width, height = height, units = "in", dpi = 650, bg = bg, limitsize = FALSE
  )
}

path_ir <- file.path(ROOT_DIR, "MWIR81DT", "MWIR81FL.dta")
assert_that(file.exists(path_ir), paste0("IR file not found at: ", path_ir))
raw_ir <- haven::read_dta(path_ir) |> janitor::clean_names()

strata_var <- pick_strata_var(raw_ir)
needed_vars <- c("v005", "v021", strata_var, "v012", "v024", "v025", "v106", "v171a", "v190", "v501", "v781")
missing_vars <- setdiff(needed_vars, names(raw_ir))
assert_that(length(missing_vars) == 0, paste0("Required variables missing: ", paste(missing_vars, collapse = ", ")))

analytic0 <- raw_ir |>
  dplyr::transmute(
    wt = v005 / 1000000,
    v021 = as.numeric(v021),
    strata = as.numeric(.data[[strata_var]]),
    age_years = as.numeric(v012),
    region = haven::as_factor(v024, levels = "labels") |> forcats::fct_drop(),
    residence = haven::as_factor(v025, levels = "labels") |> forcats::fct_drop(),
    education = haven::as_factor(v106, levels = "labels") |> forcats::fct_drop(),
    wealth_quintile = haven::as_factor(v190, levels = "labels") |> forcats::fct_drop(),
    marital_status = haven::as_factor(v501, levels = "labels") |> forcats::fct_drop(),
    v171a_num = suppressWarnings(as.numeric(v171a)),
    v781_num = suppressWarnings(as.numeric(v781))
  ) |>
  dplyr::mutate(
    internet_any = dplyr::case_when(
      is.na(v171a_num) ~ NA_real_,
      v171a_num == 0 ~ 0,
      v171a_num > 0 ~ 1,
      TRUE ~ NA_real_
    ),
    ever_tested_hiv = dplyr::case_when(
      is.na(v781_num) ~ NA_real_,
      v781_num == 0 ~ 0,
      v781_num == 1 ~ 1,
      TRUE ~ NA_real_
    ),
    internet_label = factor(
      dplyr::case_when(
        internet_any == 0 ~ "No internet use",
        internet_any == 1 ~ "Any internet use",
        TRUE ~ NA_character_
      ),
      levels = c("No internet use", "Any internet use")
    ),
    age_group = cut(
      age_years,
      breaks = c(15, 20, 25, 30, 35, 40, 45, 50),
      right = FALSE, include.lowest = TRUE,
      labels = c("15–19", "20–24", "25–29", "30–34", "35–39", "40–44", "45–49")
    )
  ) |>
  dplyr::mutate(
    education = relevel_if_present(education, c("No education", "Primary", "Secondary", "Higher", "More than secondary")),
    wealth_quintile = relevel_if_present(wealth_quintile, c("Lowest", "Second", "Middle", "Fourth", "Highest")),
    residence = relevel_if_present(residence, c("Urban", "Rural")),
    marital_status = relevel_if_present(marital_status, c("Never married", "Married", "Living together", "Married/living together", "Divorced", "Separated", "Widowed", "Divorced/separated/widowed/nullified")),
    region = relevel_if_present(region, c("Northern", "Central", "Southern"))
  )

cc_vars <- c("wt", "v021", "strata", "age_years", "region", "residence", "education",
             "wealth_quintile", "marital_status", "internet_any", "ever_tested_hiv",
             "internet_label", "age_group")

analytic <- analytic0 |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(cc_vars), ~ !is.na(.x))) |>
  droplevels()

design <- make_design(analytic, strata_var = "strata")
design_srvyr <- srvyr::as_survey_design(design)

qc_table <- tibble::tibble(
  analytic_step = c(
    "Rows in raw IR file",
    "Rows with non-missing HIV-testing outcome (V781)",
    "Rows with non-missing internet-use exposure (V171A)",
    "Rows with non-missing prespecified covariates",
    "Rows in final complete-case analytic sample",
    "Sum of normalized weights in analytic sample"
  ),
  value = c(
    nrow(raw_ir),
    sum(!is.na(analytic0$ever_tested_hiv)),
    sum(!is.na(analytic0$internet_any)),
    sum(stats::complete.cases(analytic0[, c("age_years", "education", "wealth_quintile", "residence", "marital_status", "region")])),
    nrow(analytic),
    round(sum(analytic$wt), 2)
  )
)

prev_internet <- weighted_binary_prop(design, "internet_any")
prev_tested <- weighted_binary_prop(design, "ever_tested_hiv")

qc_anchor <- tibble::tibble(
  metric = c("Weighted prevalence of any internet use", "Weighted prevalence of ever having tested for HIV"),
  estimate = c(fmt_pct(prev_internet$estimate, 1), fmt_pct(prev_tested$estimate, 1)),
  `95% CI` = c(
    paste0(fmt_pct(prev_internet$ci_low, 1), " to ", fmt_pct(prev_internet$ci_high, 1)),
    paste0(fmt_pct(prev_tested$ci_low, 1), " to ", fmt_pct(prev_tested$ci_high, 1))
  )
)

save_table_dual(qc_table, "table1_sample_flow",
                "**Table 1. Analytic sample derivation and run-level quality control**",
                "Malawi DHS 2024 women’s IR file; complete-case primary analysis",
                "The DHS sample weights are normalized, so their sum approximates the interviewed analytic sample size rather than a national population count.")

save_table_dual(qc_anchor, "table2_weighted_prevalence_anchors",
                "**Table 2. Weighted design-based prevalence anchors for the H7 analysis**",
                "Primary exposure and outcome prevalences in the final analytic sample",
                "Confidence intervals are based on the complex DHS design using Taylor-series linearization.")

age_row <- weighted_mean_wide(analytic, "age_years", "internet_label") |>
  tidyr::pivot_wider(names_from = group, values_from = value) |>
  dplyr::mutate(
    Characteristic = "Age (years), weighted mean",
    Level = "Mean",
    Overall = stats::weighted.mean(analytic$age_years, analytic$wt, na.rm = TRUE)
  ) |>
  dplyr::select(Characteristic, Level, Overall, `No internet use`, `Any internet use`)

make_baseline_block <- function(var, label) {
  weighted_dist_wide(analytic, var, "internet_label") |>
    dplyr::mutate(Characteristic = label, Level = as.character(level)) |>
    dplyr::select(Characteristic, Level, overall, `No internet use`, `Any internet use`) |>
    dplyr::rename(Overall = overall)
}

baseline_tbl <- dplyr::bind_rows(
  age_row,
  make_baseline_block("education", "Educational attainment"),
  make_baseline_block("wealth_quintile", "Household wealth quintile"),
  make_baseline_block("residence", "Residence"),
  make_baseline_block("marital_status", "Marital status"),
  make_baseline_block("region", "Region"),
  make_baseline_block("age_group", "Age group")
) |>
  dplyr::mutate(
    Overall = ifelse(Characteristic == "Age (years), weighted mean", fmt_num(Overall, 2), fmt_pct(Overall, 1)),
    `No internet use` = ifelse(Characteristic == "Age (years), weighted mean", fmt_num(`No internet use`, 2), fmt_pct(`No internet use`, 1)),
    `Any internet use` = ifelse(Characteristic == "Age (years), weighted mean", fmt_num(`Any internet use`, 2), fmt_pct(`Any internet use`, 1))
  )

save_table_dual(baseline_tbl, "table3_compact_weighted_baseline_characteristics",
                "**Table 3. Compact weighted baseline characteristics by internet-use status**",
                "All proportions are weighted column percentages; age is a weighted mean",
                "This table is intentionally wide and compact to avoid bulky long-format descriptive outputs.")

formula_primary <- ever_tested_hiv ~ internet_any + age_years + education + wealth_quintile + residence + marital_status + region
formula_crude <- ever_tested_hiv ~ internet_any
formula_pr <- ever_tested_hiv ~ internet_any + age_years + education + wealth_quintile + residence + marital_status + region
formula_age_int <- ever_tested_hiv ~ internet_any * age_group + education + wealth_quintile + residence + marital_status + region
formula_wealth_int <- ever_tested_hiv ~ internet_any * wealth_quintile + age_years + education + residence + marital_status + region
formula_res_int <- ever_tested_hiv ~ internet_any * residence + age_years + education + wealth_quintile + marital_status + region
formula_no_marital <- ever_tested_hiv ~ internet_any + age_years + education + wealth_quintile + residence + region
formula_no_residence <- ever_tested_hiv ~ internet_any + age_years + education + wealth_quintile + marital_status + region
formula_no_wealth <- ever_tested_hiv ~ internet_any + age_years + education + residence + marital_status + region

fit_crude <- survey::svyglm(formula_crude, design = design, family = quasibinomial())
fit_adj <- survey::svyglm(formula_primary, design = design, family = quasibinomial())
fit_pr <- survey::svyglm(formula_pr, design = design, family = quasipoisson(link = "log"))
fit_age_int <- survey::svyglm(formula_age_int, design = design, family = quasibinomial())
fit_wealth_int <- survey::svyglm(formula_wealth_int, design = design, family = quasibinomial())
fit_res_int <- survey::svyglm(formula_res_int, design = design, family = quasibinomial())

analytic_married <- analytic |>
  dplyr::filter(stringr::str_detect(tolower(as.character(marital_status)), "married|living together|union")) |>
  droplevels()
analytic_not_married <- analytic |>
  dplyr::filter(!stringr::str_detect(tolower(as.character(marital_status)), "married|living together|union")) |>
  droplevels()
analytic_urban <- analytic |>
  dplyr::filter(stringr::str_detect(tolower(as.character(residence)), "urban")) |>
  droplevels()
analytic_rural <- analytic |>
  dplyr::filter(stringr::str_detect(tolower(as.character(residence)), "rural")) |>
  droplevels()

fit_married <- survey::svyglm(formula_no_marital, design = make_design(analytic_married, "strata"), family = quasibinomial())
fit_not_married <- survey::svyglm(formula_no_marital, design = make_design(analytic_not_married, "strata"), family = quasibinomial())
fit_urban <- survey::svyglm(formula_no_residence, design = make_design(analytic_urban, "strata"), family = quasibinomial())
fit_rural <- survey::svyglm(formula_no_residence, design = make_design(analytic_rural, "strata"), family = quasibinomial())

main_effects <- dplyr::bind_rows(
  extract_effect(fit_crude, scale = "OR") |> dplyr::mutate(model_name = "Crude survey-weighted logistic"),
  extract_effect(fit_adj, scale = "OR") |> dplyr::mutate(model_name = "Adjusted survey-weighted logistic"),
  extract_effect(fit_pr, scale = "PR") |> dplyr::mutate(model_name = "Adjusted modified Poisson"),
  extract_effect(fit_married, scale = "OR") |> dplyr::mutate(model_name = "Adjusted logistic: married/living together subset"),
  extract_effect(fit_not_married, scale = "OR") |> dplyr::mutate(model_name = "Adjusted logistic: not currently in union subset"),
  extract_effect(fit_urban, scale = "OR") |> dplyr::mutate(model_name = "Adjusted logistic: urban subset"),
  extract_effect(fit_rural, scale = "OR") |> dplyr::mutate(model_name = "Adjusted logistic: rural subset")
) |>
  dplyr::select(model_name, effect_scale, estimate, conf_low, conf_high, p_value)

pr_row <- main_effects |> dplyr::filter(model_name == "Adjusted modified Poisson") |> dplyr::slice(1)
ev_pr <- compute_evalue_rr(pr_row$estimate, pr_row$conf_low)

main_effects_pub <- main_effects |>
  dplyr::mutate(
    effect_estimate = ifelse(effect_scale == "OR", or_ci_string(estimate, conf_low, conf_high), pr_ci_string(estimate, conf_low, conf_high)),
    p_value = fmt_p(p_value)
  ) |>
  dplyr::select(`Model specification` = model_name, `Effect scale` = effect_scale, `Estimate (95% CI)` = effect_estimate, `P value` = p_value)

save_table_dual(main_effects_pub, "table4_primary_and_sensitivity_models",
                "**Table 4. Primary and sensitivity effect estimates for internet use and ever HIV testing**",
                "Main contrasts from design-corrected regression models",
                paste0("For the adjusted prevalence ratio model, E-value = ", fmt_num(ev_pr$evalue_point, 2),
                       "; E-value for the confidence-limit closest to the null = ", fmt_num(ev_pr$evalue_ci, 2), "."))

interaction_tests <- tibble::tibble(
  domain = c("Age group", "Household wealth quintile", "Residence"),
  raw_p = c(
    safe_regterm_p(fit_age_int, ~internet_any:age_group),
    safe_regterm_p(fit_wealth_int, ~internet_any:wealth_quintile),
    safe_regterm_p(fit_res_int, ~internet_any:residence)
  )
) |>
  dplyr::mutate(fdr_q = stats::p.adjust(raw_p, method = "BH"))

fit_by_group <- function(group_var) {
  levs <- levels(analytic[[group_var]])
  subgroup_formula <- switch(
    group_var,
    "wealth_quintile" = formula_no_wealth,
    "residence" = formula_no_residence,
    "age_group" = formula_primary,
    formula_primary
  )
  
  purrr::map_dfr(levs, function(g) {
    sub_dat <- analytic |>
      dplyr::filter(.data[[group_var]] == g) |>
      droplevels()
    
    if (length(unique(sub_dat$internet_any)) < 2) {
      return(tibble::tibble(domain = group_var, stratum = g, estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_))
    }
    
    fit <- tryCatch(
      survey::svyglm(subgroup_formula, design = make_design(sub_dat, "strata"), family = quasibinomial()),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      tibble::tibble(domain = group_var, stratum = g, estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_)
    } else {
      eff <- tryCatch(extract_effect(fit, scale = "OR"), error = function(e) NULL)
      if (is.null(eff)) {
        tibble::tibble(domain = group_var, stratum = g, estimate = NA_real_, conf_low = NA_real_, conf_high = NA_real_, p_value = NA_real_)
      } else {
        eff |>
          dplyr::transmute(domain = group_var, stratum = g, estimate, conf_low, conf_high, p_value)
      }
    }
  })
}

stratum_effects <- dplyr::bind_rows(
  fit_by_group("age_group"),
  fit_by_group("wealth_quintile"),
  fit_by_group("residence")
) |>
  dplyr::mutate(
    domain = dplyr::recode(domain, age_group = "Age group", wealth_quintile = "Household wealth quintile", residence = "Residence"),
    effect = or_ci_string(estimate, conf_low, conf_high),
    p_value = fmt_p(p_value)
  )

interaction_pub <- stratum_effects |>
  dplyr::left_join(
    interaction_tests |>
      dplyr::mutate(`Interaction P value` = fmt_p(raw_p), `BH q value` = fmt_p(fdr_q)) |>
      dplyr::select(domain, `Interaction P value`, `BH q value`),
    by = "domain"
  ) |>
  dplyr::select(Domain = domain, Stratum = stratum, `Adjusted OR (95% CI)` = effect, `Stratum-specific P value` = p_value, `Interaction P value`, `BH q value`)

save_table_dual(interaction_pub, "table5_effect_modification_and_inequality_gradients",
                "**Table 5. Effect modification and inequality gradients**",
                "Stratum-specific adjusted odds ratios plus interaction tests",
                "This table addresses whether the association between digital access and prior HIV testing differs materially across age, wealth, and rural-urban residence.")

main_effects_for_plot <- main_effects |>
  dplyr::mutate(
    model_name = factor(model_name, levels = rev(model_name)),
    label = dplyr::if_else(
      effect_scale == "OR",
      paste0("OR ", or_ci_string(estimate, conf_low, conf_high)),
      paste0("PR ", pr_ci_string(estimate, conf_low, conf_high))
    )
  )

modal_value <- function(x) {
  if (is.factor(x)) names(sort(table(x), decreasing = TRUE))[1] else stats::median(x, na.rm = TRUE)
}
modal_edu <- modal_value(analytic$education)
modal_wealth <- modal_value(analytic$wealth_quintile)
modal_res <- modal_value(analytic$residence)
modal_marital <- modal_value(analytic$marital_status)
modal_region <- modal_value(analytic$region)

new_overall <- tibble::tibble(
  internet_any = c(0, 1),
  age_years = mean(analytic$age_years, na.rm = TRUE),
  education = factor(modal_edu, levels = levels(analytic$education)),
  wealth_quintile = factor(modal_wealth, levels = levels(analytic$wealth_quintile)),
  residence = factor(modal_res, levels = levels(analytic$residence)),
  marital_status = factor(modal_marital, levels = levels(analytic$marital_status)),
  region = factor(modal_region, levels = levels(analytic$region)),
  internet_label = factor(c("No internet use", "Any internet use"), levels = levels(analytic$internet_label)),
  age_group = cut(
    mean(analytic$age_years, na.rm = TRUE),
    breaks = c(15, 20, 25, 30, 35, 40, 45, 50),
    right = FALSE, include.lowest = TRUE, labels = levels(analytic$age_group)
  )
)
pred_overall <- predict_svyglm_safe(fit_adj, new_overall, type = "link")

age_grid <- expand.grid(age_years = seq(15, 49, by = 1), internet_any = c(0, 1)) |>
  dplyr::mutate(
    education = factor(modal_edu, levels = levels(analytic$education)),
    wealth_quintile = factor(modal_wealth, levels = levels(analytic$wealth_quintile)),
    residence = factor(modal_res, levels = levels(analytic$residence)),
    marital_status = factor(modal_marital, levels = levels(analytic$marital_status)),
    region = factor(modal_region, levels = levels(analytic$region)),
    age_group = cut(age_years, breaks = c(15, 20, 25, 30, 35, 40, 45, 50), right = FALSE, include.lowest = TRUE, labels = levels(analytic$age_group)),
    internet_label = factor(ifelse(internet_any == 1, "Any internet use", "No internet use"), levels = levels(analytic$internet_label))
  )
pred_age <- predict_svyglm_safe(fit_adj, age_grid, type = "link")

wealth_res_grid <- expand.grid(
  wealth_quintile = levels(analytic$wealth_quintile),
  residence = levels(analytic$residence),
  internet_any = c(0, 1),
  stringsAsFactors = FALSE
) |>
  dplyr::mutate(
    age_years = mean(analytic$age_years, na.rm = TRUE),
    education = factor(modal_edu, levels = levels(analytic$education)),
    marital_status = factor(modal_marital, levels = levels(analytic$marital_status)),
    region = factor(modal_region, levels = levels(analytic$region)),
    wealth_quintile = factor(wealth_quintile, levels = levels(analytic$wealth_quintile)),
    residence = factor(residence, levels = levels(analytic$residence)),
    internet_label = factor(ifelse(internet_any == 1, "Any internet use", "No internet use"), levels = levels(analytic$internet_label)),
    age_group = cut(age_years, breaks = c(15, 20, 25, 30, 35, 40, 45, 50), right = FALSE, include.lowest = TRUE, labels = levels(analytic$age_group))
  )
pred_wealth_res <- predict_svyglm_safe(fit_adj, wealth_res_grid, type = "link")

region_prev <- design_srvyr |>
  dplyr::group_by(region) |>
  dplyr::summarise(
    internet_prev = srvyr::survey_mean(internet_any, vartype = "ci", proportion = TRUE),
    hiv_prev = srvyr::survey_mean(ever_tested_hiv, vartype = "ci", proportion = TRUE)
  ) |>
  dplyr::mutate(region = forcats::fct_reorder(region, hiv_prev))

observed_prev <- design_srvyr |>
  dplyr::group_by(internet_label) |>
  dplyr::summarise(hiv_prev = srvyr::survey_mean(ever_tested_hiv, vartype = "ci", proportion = TRUE))

theme_pub <- function() {
  ggplot2::theme_minimal(base_size = 13, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 10)),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9, colour = "gray30", lineheight = 1.1),
      plot.margin = ggplot2::margin(18, 28, 18, 18)
    )
}

fig1 <- ggplot2::ggplot(observed_prev, ggplot2::aes(x = internet_label, y = hiv_prev, fill = internet_label)) +
  ggplot2::geom_col(width = 0.60, colour = "white") +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = hiv_prev_low, ymax = hiv_prev_upp), width = 0.08, linewidth = 0.7, colour = "gray20") +
  ggplot2::geom_text(ggplot2::aes(label = fmt_pct(hiv_prev, 1)), vjust = -0.5, size = 4.0, fontface = "bold") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::scale_fill_manual(values = c("No internet use" = "#8F2D56", "Any internet use" = "#20639B")) +
  ggplot2::labs(
    title = "Observed HIV-testing coverage is higher among women with internet access",
    subtitle = "Complex-survey weighted prevalence of ever having tested for HIV",
    x = NULL, y = "Weighted prevalence of ever HIV testing",
    caption = "Bars show weighted survey estimates; whiskers show 95% confidence intervals."
  ) + theme_pub() + ggplot2::theme(legend.position = "none")
save_plot650(fig1, "figure1_observed_weighted_hiv_testing_by_internet.png", 9.2, 7.0)

fig2 <- ggplot2::ggplot(main_effects_for_plot, ggplot2::aes(x = estimate, y = model_name, xmin = conf_low, xmax = conf_high)) +
  ggplot2::geom_vline(xintercept = 1, linetype = 2, linewidth = 0.5, colour = "gray45") +
  ggplot2::geom_errorbarh(height = 0.16, linewidth = 0.8, colour = "gray30") +
  ggplot2::geom_point(size = 3.0, shape = 21, stroke = 0.7, fill = "#4E79A7", colour = "black") +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = label),
    direction = "y", nudge_x = 0.08, hjust = 0, size = 3.2,
    box.padding = 0.20, point.padding = 0.10, segment.alpha = 0.5, min.segment.length = 0
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "The internet-use signal is stable across core sensitivity analyses",
    subtitle = "Effect estimates for the contrast any internet use vs no internet use",
    x = "Effect estimate on logarithmic scale", y = NULL,
    caption = "Points denote point estimates and horizontal lines denote 95% confidence intervals."
  ) + theme_pub()
save_plot650(fig2, "figure2_forest_primary_and_sensitivity_models.png", 12.0, 7.8)

fig3 <- ggplot2::ggplot(pred_age, ggplot2::aes(x = age_years, y = prob, colour = internet_label, fill = internet_label)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = prob_lo, ymax = prob_hi), alpha = 0.14, linewidth = 0) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::scale_x_continuous(breaks = seq(15, 49, by = 5)) +
  ggplot2::scale_colour_manual(values = c("No internet use" = "#8F2D56", "Any internet use" = "#1B9E77")) +
  ggplot2::scale_fill_manual(values = c("No internet use" = "#8F2D56", "Any internet use" = "#1B9E77")) +
  ggplot2::labs(
    title = "Digital exclusion remains consequential across the reproductive-age life course",
    subtitle = "Adjusted predicted probability of ever HIV testing by age and internet-use status",
    x = "Age (years)", y = "Predicted probability of ever HIV testing",
    caption = "Predictions come from the adjusted survey-weighted logistic model with covariates fixed to a standardized reference profile."
  ) + theme_pub()
save_plot650(fig3, "figure3_adjusted_probability_by_age_and_internet.png", 11.0, 7.2)

fig4 <- ggplot2::ggplot(pred_wealth_res, ggplot2::aes(x = wealth_quintile, y = prob, colour = internet_label, group = internet_label)) +
  ggplot2::geom_line(linewidth = 1.1, position = ggplot2::position_dodge(width = 0.18)) +
  ggplot2::geom_point(size = 2.8, position = ggplot2::position_dodge(width = 0.18)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = prob_lo, ymax = prob_hi), width = 0.08, linewidth = 0.6, position = ggplot2::position_dodge(width = 0.18)) +
  ggplot2::facet_wrap(~ residence, nrow = 1) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::scale_colour_manual(values = c("No internet use" = "#B23A48", "Any internet use" = "#2A9D8F")) +
  ggplot2::labs(
    title = "Digital access intersects with wealth and place to shape testing inequality",
    subtitle = "Adjusted predicted probabilities across wealth quintiles, stratified by urban-rural residence",
    x = "Household wealth quintile", y = "Predicted probability of ever HIV testing",
    caption = "Non-parallel patterns across wealth and residence provide a direct visual test of inequality amplification."
  ) + theme_pub() + ggplot2::theme(panel.grid.major.x = ggplot2::element_line(colour = "gray90"))
save_plot650(fig4, "figure4_wealth_residence_digital_inequality.png", 12.0, 7.5)

fig5 <- ggplot2::ggplot(region_prev, ggplot2::aes(x = internet_prev, y = hiv_prev, label = region)) +
  ggplot2::geom_point(size = 3.2, colour = "#1D3557", alpha = 0.95) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, colour = "#E76F51", linetype = 2) +
  ggrepel::geom_text_repel(size = 3.2, min.segment.length = 0, box.padding = 0.25, point.padding = 0.20) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, max(region_prev$internet_prev, na.rm = TRUE) * 1.15)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(min(region_prev$hiv_prev, na.rm = TRUE) * 0.95, 1)) +
  ggplot2::labs(
    title = "Regions with higher digital penetration also show higher HIV-testing coverage",
    subtitle = "Ecologic alignment of regional internet access and prior HIV testing",
    x = "Regional prevalence of any internet use", y = "Regional prevalence of ever HIV testing",
    caption = "This ecologic figure is descriptive rather than causal; it complements the individual-level design-based models."
  ) + theme_pub()
save_plot650(fig5, "figure5_region_level_internet_vs_hiv_testing.png", 10.8, 7.4)

pred_export <- pred_wealth_res |>
  dplyr::transmute(
    wealth_quintile, residence, internet_label,
    predicted_probability = fmt_pct(prob, 1),
    ci_low = fmt_pct(prob_lo, 1),
    ci_high = fmt_pct(prob_hi, 1)
  )
save_table_dual(pred_export, "supplement_adjusted_predicted_probabilities_wealth_residence_internet",
                "**Supplementary table. Adjusted predicted probabilities by wealth, residence, and internet-use status**",
                "Derived from the adjusted survey-weighted logistic model")

coef_primary <- coef_table(fit_adj, exponentiate = TRUE) |>
  dplyr::mutate(term_label = collapse_term_labels(term)) |>
  dplyr::select(Term = term_label, `Adjusted OR` = estimate_fmt, `95% CI low` = conf_low_fmt, `95% CI high` = conf_high_fmt, `P value` = p_value_fmt)
save_table_dual(coef_primary, "supplement_full_primary_model_coefficients",
                "**Supplementary table. Full adjusted model coefficients**",
                "Exponentiated coefficients from the primary survey-weighted logistic model")

saveRDS(list(
  fit_crude = fit_crude, fit_adj = fit_adj, fit_pr = fit_pr,
  fit_age_int = fit_age_int, fit_wealth_int = fit_wealth_int, fit_res_int = fit_res_int,
  fit_married = fit_married, fit_not_married = fit_not_married, fit_urban = fit_urban, fit_rural = fit_rural
), file = file.path(DIR_RDS, "model_objects_h7_final.rds"))

saveRDS(analytic, file = file.path(DIR_RDS, "analytic_dataset_complete_case_h7_final.rds"))
saveRDS(list(
  main_effects = main_effects, interaction_tests = interaction_tests,
  pred_overall = pred_overall, pred_age = pred_age, pred_wealth_res = pred_wealth_res,
  region_prev = region_prev, observed_prev = observed_prev
), file = file.path(DIR_RDS, "derived_results_h7_final.rds"))

manifest <- tibble::tibble(
  item = c("Analysis root directory", "Primary data file", "Analytic sample size", "Detected strata variable", "Output directory", "R version"),
  value = c(ROOT_DIR, path_ir, nrow(analytic), "strata", DIR_OUT, R.version.string)
)
save_table_dual(manifest, "manifest_reproducibility", "**Reproducibility manifest**", "Essential run-level metadata")
writeLines(capture.output(sessionInfo()), con = file.path(DIR_LOG, "sessionInfo.txt"))

key_adj <- extract_effect(fit_adj, scale = "OR")
key_pr <- extract_effect(fit_pr, scale = "PR")
cat("\nAdjusted OR for any internet use vs none: ", or_ci_string(key_adj$estimate, key_adj$conf_low, key_adj$conf_high), " | p = ", fmt_p(key_adj$p_value), "\n", sep = "")
cat("Adjusted PR for any internet use vs none: ", pr_ci_string(key_pr$estimate, key_pr$conf_low, key_pr$conf_high), " | p = ", fmt_p(key_pr$p_value), "\n", sep = "")
cat("Outputs written to: ", DIR_OUT, "\n", sep = "")
