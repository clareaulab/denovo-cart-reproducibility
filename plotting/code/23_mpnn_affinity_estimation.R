library(tidyr)
library(dplyr)
library(stringr)
library(purrr)

## Load file
experimental_readout = read.csv("../data/2025_10_21_bcma_cd19_cd22_readout_combined.csv",check.names = FALSE)

mapped_id = read.csv("../data/2025_10_31_denovo_CAR_T_ID_conversions.csv")
mapped_id$order_binders = str_replace_all(mapped_id$order_binders, pattern = "_CD22", replacement = "")
experimental_readout_mapped = merge(experimental_readout,mapped_id,by.x="binder_name",by.y="order_binders",all.x=TRUE)
experimental_readout_mapped = experimental_readout_mapped %>% mutate(
  new_name_final = case_when(
    is.na(new_name_final) ~ str_split_i(binder_name, pattern = "_", i = -1),
    TRUE ~ new_name_final
  )
)
# Based on co-culture experiment results, compute summarizing metrics

## BCMA positive lines
experimental_readout$mean_BCMA_positive_activity = ((
  experimental_readout$`coculture DG-75`+
    experimental_readout$`coculture K562 BCMA OE` +
    experimental_readout$`coculture Raji` +
    experimental_readout$`coculture Ramos` +
    experimental_readout$`coculture RPMI 8226`
)/5)
experimental_readout$max_BCMA_negative_activity = pmax(
  experimental_readout$`coculture CAR Only`,
  experimental_readout$`coculture SET-2`,
  experimental_readout$`coculture U937`,
  experimental_readout$`coculture A549`,
  experimental_readout$`coculture Hela`,
  experimental_readout$`coculture K562`,
  experimental_readout$`coculture MDA-MB-468`
)

## CD19 positive lines
experimental_readout$mean_CD19_positive_activity = ((
  experimental_readout$`coculture K562 CD19 OE`+
    experimental_readout$`coculture Raji` #+
  #experimental_readout$`coculture Raji CD19 OE`
)/2)
experimental_readout$max_CD19_negative_activity = pmax(
  experimental_readout$`coculture CAR Only`,
  experimental_readout$`coculture Raji CD19 KO`,
  experimental_readout$`coculture K562`
)

## CD22 positive lines
experimental_readout$mean_CD22_positive_activity = ((
  experimental_readout$`coculture K562 CD22 OE`+
    experimental_readout$`coculture Raji` +
    experimental_readout$`coculture Nalm6`+
    experimental_readout$`coculture Nalm6 CD22 OE`
)/4)

experimental_readout$max_CD22_negative_activity = pmax(
  experimental_readout$`coculture CAR Only`,
  experimental_readout$`coculture RPMI 8226`,
  experimental_readout$`coculture K562`
)

## Compute ratio-based metrics
## BCMA
experimental_readout$BCMA_activity_gain = experimental_readout$mean_BCMA_positive_activity - experimental_readout$`coculture CAR Only`
experimental_readout$BCMA_activity_ratio = experimental_readout$mean_BCMA_positive_activity / experimental_readout$`coculture CAR Only`

## CD19
experimental_readout$CD19_activity_gain = experimental_readout$mean_CD19_positive_activity - experimental_readout$`coculture CAR Only`
experimental_readout$CD19_activity_ratio = experimental_readout$mean_CD19_positive_activity / experimental_readout$`coculture CAR Only`

## CD22
experimental_readout$CD22_activity_gain = experimental_readout$mean_CD22_positive_activity - experimental_readout$`coculture CAR Only`
experimental_readout$CD22_activity_ratio = experimental_readout$mean_CD22_positive_activity / experimental_readout$`coculture CAR Only`

## Overall
experimental_readout = experimental_readout %>%
  mutate(
    overall_activity_gain = case_when(
      binder_target == "BCMA" ~ BCMA_activity_gain,
      binder_target == "CD19"  ~ CD19_activity_gain,
      binder_target == "CD22" ~ CD22_activity_gain
    ),
    overall_activity_ratio = case_when(
      binder_target == "BCMA" ~ BCMA_activity_ratio,
      binder_target == "CD19"  ~ CD19_activity_ratio,
      binder_target == "CD22" ~ CD22_activity_ratio
    ),
    overall_mean_positive_activity = case_when(
      binder_target == "BCMA" ~ mean_BCMA_positive_activity,
      binder_target == "CD19"  ~ mean_CD19_positive_activity,
      binder_target == "CD22" ~ mean_CD22_positive_activity
    ),
    overall_max_negative_activity = case_when(
      binder_target == "BCMA" ~ max_BCMA_negative_activity,
      binder_target == "CD19"  ~ max_CD19_negative_activity,
      binder_target == "CD22" ~ max_CD22_negative_activity
    )
  )

#' Estimates binding affinity (Kd) using a hybrid approach.
#'
#' This function encapsulates all logic in a single function, first attempting
#' a nonlinear least squares (NLS) fit for each binder. If the fit is reliable,
#' it uses the NLS-estimated Kd. Otherwise, it falls back to a single-point
#' estimation using a provided reference binder to normalize the data.
#'
#' @param df A data frame with at least "binder_name", "binder_experiment",
#'   and columns for staining values at different concentrations (e.g., "staining_0_nM").
#' @param ref_mapping A named list mapping each unique "binder_experiment"
#'   to a list containing its specific reference binder details and single-point
#'   concentration. The structure should be:
#'   list("experiment_name" = list(
#'     reference_name = "ref_binder",
#'     reference_kd = numeric_value,
#'     ref_conc = numeric_value
#'   ), ...)
#' @return A data frame containing the binder name, the final estimated Kd,
#'   the NLS-derived Kd, the single-point derived Kd, a quality flag for the fit,
#'   and the R-squared value from the curve fit.
#' @examples
#' # Sample data frame with a reference binder and multiple experiments
#' binder_data <- data.frame(
#'   binder_name = c("BCMA_UTD", "BCMA_NB", "BCMA_Abecma", "BCMA_l59_int_wildtype", "BCMA_exp2_good", "BCMA_exp2_ref"),
#'   binder_experiment = c("BCMA", "BCMA", "BCMA", "BCMA", "BCMA_Exp2", "BCMA_Exp2"),
#'   staining_0_nM = c(1.09, 0.91, 0.94, 1.53, 0.5, 0.6),
#'   staining_0.1_nM = c(1.15, 1.24, 3.24, 2.08, 1.5, 1.8),
#'   staining_1_nM = c(1.09, 0.93, 80.50, 53.60, 50.0, 55.0),
#'   staining_10_nM = c(1.41, 0.95, 98.70, 99.60, 95.0, 98.0),
#'   staining_100_nM = c(1.97, 1.07, 98.80, 100.00, 98.0, 99.5),
#'   check.names = FALSE
#' )
#'
#' # Define the reference mapping for each experiment
#' # Note: BCMA_Abecma is the reference for BCMA experiment,
#' # and BCMA_exp2_ref is the reference for BCMA_Exp2 experiment.
#' ref_mapping <- list(
#'   "BCMA" = list(
#'     reference_name = "BCMA_Abecma",
#'     reference_kd = 2,
#'     ref_conc = 10
#'   ),
#'   "BCMA_Exp2" = list(
#'     reference_name = "BCMA_exp2_ref",
#'     reference_kd = 5,
#'     ref_conc = 100
#'   )
#' )
#'
#' # Run the function
#' estimate_hybrid_kd(df = binder_data, ref_mapping = ref_mapping)
estimate_hybrid_kd = function(df, ref_mapping) {
  
  # Helper function: Defines the binding model for nonlinear curve fitting.
  binding_model <- function(conc, Kd, max_signal) {
    (max_signal * conc) / (Kd + conc)
  }
  
  # Helper function: Inlined NLS Fit Calculation
  calculate_nls_fit <- function(fit_data) {
    nls_kd <- NA
    nls_r2 <- NA
    if (nrow(fit_data) >= 2) {
      tryCatch({
        fit_result <- nls(staining ~ binding_model(concentrations, Kd, max_signal),
                          data = fit_data,
                          start = list(Kd = 1, max_signal = max(fit_data$staining, na.rm = TRUE)))
        
        nls_kd <- coef(fit_result)["Kd"]
        ss_res <- sum(residuals(fit_result)^2)
        ss_tot <- sum((fit_data$staining - mean(fit_data$staining))^2)
        nls_r2 <- 1 - (ss_res / ss_tot)
      }, error = function(e) {
        message(paste("NLS fit failed:", e$message))
      })
    }
    return(list(kd = nls_kd, r2 = nls_r2))
  }
  
  # Helper function: Inlined Single-Point Estimate Calculation
  calculate_single_point_kd <- function(ref_data, target_data, ref_kd, ref_conc) {
    sp_kd <- NA
    ref_staining_col <- paste0("staining_", ref_conc, "_nM")
    ref_staining <- as.numeric(ref_data[[ref_staining_col]])
    
    if (!is.na(ref_staining) && ref_staining > 0) {
      max_signal_ref_est <- (ref_staining * (ref_kd + ref_conc)) / ref_conc
      binder_staining_col <- paste0("staining_", ref_conc, "_nM")
      binder_staining <- as.numeric(target_data[[binder_staining_col]])
      
      if (!is.na(binder_staining) && binder_staining > 0 && binder_staining <= max_signal_ref_est) {
        sp_kd <- (ref_conc * (max_signal_ref_est - binder_staining)) / binder_staining
      }
    }
    return(sp_kd)
  }
  
  # Find the staining columns and concentrations
  staining_cols <- grep("staining_.*_nM", names(df), value = TRUE)
  concentrations <- as.numeric(gsub("staining_|_nM", "", staining_cols))
  
  # Initialize the final results data frame
  final_results_df <- data.frame(
    binder_name = character(),
    estimated_Kd = numeric(),
    nls_kd = numeric(),
    nls_kd_scaled = numeric(), # New column for scaled NLS KD
    sp_kd = numeric(),
    sp_kd_scaled = numeric(),  # New column for scaled SP KD
    kd_fit_quality = character(),
    kd_R2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Split the data by experiment
  df_by_exp <- split(df, df$binder_experiment)
  
  for (exp_name in names(df_by_exp)) {
    exp_df <- df_by_exp[[exp_name]]
    
    # Get the specific reference information for this experiment
    current_ref_info <- ref_mapping[[exp_name]]
    if (is.null(current_ref_info)) {
      warning(paste("No reference mapping found for experiment:", exp_name, ". Skipping..."))
      next
    }
    
    current_ref_name <- current_ref_info$reference_name
    current_ref_kd <- current_ref_info$reference_kd
    current_ref_conc <- current_ref_info$ref_conc
    
    ref_data <- exp_df[exp_df$binder_name == current_ref_name, , drop = FALSE]
    
    if (nrow(ref_data) == 0) {
      warning(paste("Reference binder '", current_ref_name, "' not found in experiment:", exp_name, ". Skipping..."))
      next
    }
    
    # Calculate the scaling factor for the current experiment
    ref_staining_values <- as.numeric(ref_data[1, staining_cols])
    ref_fit_data <- data.frame(concentrations = concentrations, staining = ref_staining_values)
    ref_fit_data <- na.omit(ref_fit_data)
    
    ref_nls_fit <- calculate_nls_fit(ref_fit_data)
    scaling_factor <- NA
    if (!is.na(ref_nls_fit$kd) && ref_nls_fit$kd > 0) {
      scaling_factor <- current_ref_kd / ref_nls_fit$kd
    } else {
      warning(paste("NLS fit for reference binder '", current_ref_name, "' failed or returned a non-positive value. Skipping scaling for this experiment."))
    }
    
    # Process each binder in the current experiment
    for (i in 1:nrow(exp_df)) {
      binder_name <- exp_df$binder_name[i]
      staining_values <- as.numeric(exp_df[i, staining_cols])
      fit_data <- data.frame(concentrations = concentrations, staining = staining_values)
      fit_data <- na.omit(fit_data)
      
      # Step 1: NLS Fit
      nls_fit_result <- calculate_nls_fit(fit_data)
      
      # Step 2: Single-Point Estimate
      sp_kd <- calculate_single_point_kd(ref_data, exp_df[i, ], current_ref_kd, current_ref_conc)
      
      # Step 3: Scale the NLS and SP KDs if a valid scaling factor exists
      nls_kd_scaled <- NA
      sp_kd_scaled <- NA
      if (!is.na(scaling_factor)) {
        if (!is.na(nls_fit_result$kd)) {
          nls_kd_scaled <- nls_fit_result$kd * scaling_factor
        }
        if (!is.na(sp_kd)) {
          sp_kd_scaled <- sp_kd * scaling_factor
        }
      }
      
      # Step 4: Select Best Estimate (using the scaled values)
      if (binder_name == current_ref_name) {
        estimated_kd <- current_ref_kd
        fit_quality <- "Known"
      } else if (!is.na(nls_fit_result$kd) && !is.na(nls_fit_result$r2) && nls_fit_result$r2 > 0.95 && nls_fit_result$kd > 0) {
        estimated_kd <- nls_kd_scaled
        fit_quality <- "Good (Curve Fit)"
      } else {
        estimated_kd <- sp_kd_scaled
        fit_quality <- "Compromise (Single Point Estimate)"
        if (is.na(sp_kd_scaled)) {
          fit_quality <- "Failure"
        }
      }
      
      # Add to final results
      final_results_df <- rbind(final_results_df, data.frame(
        binder_name = binder_name,
        estimated_Kd = estimated_kd,
        nls_kd = nls_fit_result$kd,
        nls_kd_scaled = nls_kd_scaled, # New column for scaled NLS KD
        sp_kd = sp_kd,
        sp_kd_scaled = sp_kd_scaled,  # New column for scaled SP KD
        kd_fit_quality = fit_quality,
        kd_R2 = nls_fit_result$r2,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(final_results_df)
}

reference_mapping = list(
  "BCMA" = list(
    reference_name = "BCMA_Abecma",
    reference_kd = 10,
    ref_conc = 10
  ),
  "CD19_l82" = list(
    reference_name = "CD19_l82_FMC63",
    reference_kd = 5, # Value from: Seigner J et al Sci Rep 2023
    ref_conc = 100
  ),
  "CD19_l112" = list(
    reference_name = "CD19_l112_FMC63",
    reference_kd = 5, # Value from: Seigner J et al Sci Rep 2023
    ref_conc = 100
  ),
  "CD19_l115" = list(
    reference_name = "CD19_l115_FMC63",
    reference_kd = 5, # Value from: Seigner J et al Sci Rep 2023
    ref_conc = 100
  ),
  "CD22" = list(
    reference_name = "CD22_m971",
    reference_kd = 24.8, # Value from: June Ereño-Orbea et al JBC 2021, Figure 4
    ref_conc = 0.1
  )
)


hybrid_kd_estimation_df = estimate_hybrid_kd(experimental_readout,reference_mapping)
hybrid_kd_estimation_df$Kd_estimate_nM_clipped = pmax(hybrid_kd_estimation_df$estimated_Kd,1000)
hybrid_kd_estimation_df$Kd_estimate_molar = hybrid_kd_estimation_df$estimated_Kd*1e-9

hybrid_kd_estimation_df$neg_log10_Kd_M = -log10(hybrid_kd_estimation_df$Kd_estimate_molar) 
hybrid_kd_estimation_df$clipped_neg_log10_Kd_M = pmax(hybrid_kd_estimation_df$neg_log10_Kd_M,6)
hybrid_kd_estimation_df

experimental_readout_merged = merge(experimental_readout_mapped,hybrid_kd_estimation_df,by="binder_name",all.x=TRUE)

if (FALSE) {
  ## These are the affinity estimations used in the main ananlysis file
  write.csv(experimental_readout_merged,"../data/mpnn_affinity_estimations.csv",quote=TRUE,row.names = FALSE)
}

experimental_palette = c(
  "interface" = "dodgerblue3",
  "non-interface" = "firebrick",
  "parental" = "black",
  "pos_control" = "orange",
  "neg_control" = "gray"
)

## This is the part that visualizes this
plot_staining_lineplot <- function(data, experiment_name, y_value = "staining_value", title = NULL) {
  # Filter the data for the specified experiment
  plot_data <- data %>%
    filter(binder_experiment == experiment_name)
  
  # Extract unique numeric concentrations to set levels and breaks
  concentration_levels <- sort(unique(plot_data$concentration_nM_numeric))
  
  # Create the ggplot object
  p <- ggplot(plot_data, aes(
    x = factor(concentration_nM_numeric, levels = concentration_levels),
    y = .data[[y_value]],
    group = binder_name,
    color = binder_evolution_type,
    #color = experimental_setting
  )) +
    geom_line() +
    geom_point() +
    scale_x_discrete(
      breaks = concentration_levels,
      labels = as.character(concentration_levels)
    ) +
    labs(
      x = "Concentration (nM)",
      y = "Staining Value",
      color = "Evolution Strategy",
      title = experiment_name
      #title = paste0("Staining Outcomes for ", experiment_name)
    ) +
    scale_color_manual(values=experimental_palette) +
    theme_minimal()
  
  return(p)
}

staining_cols = c("staining_0_nM", "staining_0.1_nM", "staining_1_nM", "staining_10_nM", "staining_100_nM")

# Prepare the data for plotting
annotated_df_long = experimental_readout_merged %>%
  pivot_longer(
    cols = all_of(staining_cols),
    names_to = "concentration_nM",
    values_to = "staining_value"
  ) %>%
  mutate(concentration_nM_numeric = as.numeric(gsub("staining_|_nM", "", concentration_nM)))


bcma_staining_lineplot = plot_staining_lineplot(annotated_df_long, "BCMA") + pretty_plot()+ L_border() + theme(legend.position = "none")
bcma_staining_lineplot
cd19_l82_staining_lineplot = plot_staining_lineplot(annotated_df_long, "CD19_l82")+ pretty_plot()+ L_border() + theme(legend.position = "none")
cd19_l82_staining_lineplot
cd19_l112_staining_lineplot = plot_staining_lineplot(annotated_df_long, "CD19_l112")+ pretty_plot()+ L_border() + theme(legend.position = "none")
cd19_l112_staining_lineplot
cd19_l115_staining_lineplot = plot_staining_lineplot(annotated_df_long, "CD19_l115")+ pretty_plot()+ L_border() + theme(legend.position = "none")
cd19_l115_staining_lineplot
cd22_staining_lineplot = plot_staining_lineplot(annotated_df_long, "CD22")+ pretty_plot()+ L_border() + theme(legend.position = "none")

lineplot_legend = get_legend(bcma_staining_lineplot)
# lineplots_without_legends = list(
#   bcma_staining_lineplot + theme(legend.position = "none"),
#   #cd19_staining_lineplot + theme(legend.position = "none"),
#   cd22_staining_lineplot + theme(legend.position = "none")
# )

lineplots_without_legends = list(
  bcma_staining_lineplot + theme(legend.position = "none",plot.title = element_text(hjust = 0.5)),# + labs(title=element_blank()),
  #cd19_l82_staining_lineplot + theme(legend.position = "none"),
  cd19_l112_staining_lineplot + theme(legend.position = "none",plot.title = element_text(hjust = 0.5)),
  cd19_l115_staining_lineplot + theme(legend.position = "none",plot.title = element_text(hjust = 0.5)),
  cd22_staining_lineplot + theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
)

# staining_lineplots_combined = cowplot::plot_grid(
#   cowplot::plot_grid(plotlist = lineplots_without_legends, ncol = 5, align = "h"),
#   #lineplot_legend,
#   ncol = 2,
#   rel_widths = c(1, 0.1) # Adjust the width ratio as needed
# )
# staining_lineplots_combined

staining_lineplots_combined = cowplot::plot_grid(plotlist = lineplots_without_legends, ncol = 4, align = "h")
staining_lineplots_combined
cowplot::ggsave2("../plots/staining_lineplots.pdf",staining_lineplots_combined,dpi=300,height=1.8,width=1.8*4+0.3)

## Make affinity estimation barplots
experimental_readout_merged_filtered = experimental_readout_merged %>% filter(binder_experiment != "CD19_l82")

library(tidytext)
estimated_affinity_plot = ggplot(
  experimental_readout_merged_filtered, # %>% filter(binder_experiment!="CD19"),
  aes(
    x = neg_log10_Kd_M,
    y = tidytext::reorder_within(new_name_final, neg_log10_Kd_M, binder_experiment),
    #y = reorder(binder_name, neg_log10_Kd_M),
    fill = binder_evolution_type
  )
) +
  geom_bar(stat = "identity") +
  tidytext::scale_y_reordered() +
  # Use geom_text_repel for intelligent, non-overlapping annotations.
  # The `label` is created with a conditional statement to show "Kd >= 1000 nM"
  # when the affinity is very low.
  geom_text_repel(
    aes(
      label = ifelse(
        estimated_Kd >= 1000,
        paste0("≥ 1000 nM"),
        paste0(sprintf("%.2f", estimated_Kd), " nM")
      )
    ),
    # Position the text slightly to the right of the bar.
    hjust = 0,
    nudge_x = 0.5, # Nudge the text right by 0.5 units on the x-axis.
    direction = "y", # Only allow movement along the y-axis.
    size = 3,
    force = 1, # Increase the repulsion force to ensure separation.
    min.segment.length = Inf # Set this to a very large value to hide the segment.
  ) +
  
  # Separate plots for each experiment.
  facet_wrap(~ binder_experiment, scales = "free_y", ncol = 4) +
  
  # Manually set the x-axis limits.
  # Note: `coord_cartesian` is used to zoom in without discarding data,
  # so annotations outside the viewable area will still be plotted.
  coord_cartesian(xlim = c(5, 10.5)) +
  
  # Customize the x-axis ticks and labels.
  scale_x_continuous(
    breaks = c(5, 6, 7, 8, 9, 10),
    labels = function(pKd) {
      kd_nM <- 10^(-pKd) * 1e9
      paste0(scales::comma(kd_nM), "")
    }
  ) +
  
  # Set the color palette for the bars.
  #scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values=experimental_palette) +
  
  # Set the labels for the axes and legend.
  labs(
    x = "Kd (nM, lower = tighter binding)",
    y = "Binder Name",
    fill = "Experimental Setting"
  ) +
  
  # Use a clean theme.
  theme_classic() +
  theme(legend.position = "none")

# Print the plot
estimated_affinity_plot
cowplot::ggsave2("../plots/estimated_affinity_plot.pdf",estimated_affinity_plot,dpi=300,height=1.8*4,width=1.8*10)


