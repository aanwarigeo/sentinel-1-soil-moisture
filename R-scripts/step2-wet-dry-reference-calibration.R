################################################################################
# Step 2: Wet & Dry Reference Calibration
################################################################################
# 
# This script outlines a step-by-step workflow for step 2 to estimate field-scale
# soil moisture using sentinel-1 GRD SAR data.
#
# Key steps in this workflow:
# 1. Calculated dry reference backscatter intensity (σ°dry)
# 2. Calculated backscatter change (Δσ) relative to dry reference (σ°dry)
# 3. Determine relationship between maximum backscatter change and 
# Dual-pol Radar Vegetation Index (DpRVIc) 

# For a detailed explanation of each step of the methodology, please refer to the 
# research paper by Bhogapurapu et al. (2022).

# Author: aanwari
# Date: Feb 28, 2025

################################################################################


# ==============================================
# 1. Load Required Libraries
# ==============================================
load_required_libraries <- function() {
  # List of required packages
  packages <- c(
    "terra",      # For raster processing
    "sf",         # For spatial data handling
    "openxlsx",   # For Excel file operations
    "ggplot2",    # For plotting
    "dplyr",      # For data manipulation
    "tidyr",      # For data tidying
    "lubridate"   # For date handling
  )
  
  # Install missing packages
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      message("Installing package: ", pkg)
      install.packages(pkg)
      library(pkg, character.only = TRUE)
    }
  }
}

# Run the function to load required libraries
load_required_libraries()



# ==============================================
# 2. Extract Dates from Sentinel-1 GeoTIFF files downloaded in step 1
# ==============================================

extract_dates_from_filenames <- function(directory) {
  # Check if directory exists
  if (!dir.exists(directory)) {
    stop("Input directory does not exist: ", directory)
  }
  
  # Get list of GeoTIFF files
  tiff_files <- list.files(directory, pattern = "\\.tif$", full.names = FALSE)
  if (length(tiff_files) == 0) {
    stop("No GeoTIFF files found in: ", directory)
  }
  
  # Function to extract datetime from a single filename
  extract_datetime <- function(filename) {
    # Extract datetime string using regex pattern
    datetime_str <- stringr::str_extract(filename, "\\d{8}T\\d{6}(?=_)")
    
    # Check if datetime was found
    if (is.na(datetime_str)) {
      warning("Datetime not found in filename: ", filename)
      return(NA)
    }
    
    # Convert to POSIXct datetime
    tryCatch({
      as.POSIXct(datetime_str, format = "%Y%m%dT%H%M%S", tz = "UTC")
    }, error = function(e) {
      warning("Error parsing datetime from filename: ", filename, "\n", e$message)
      NA
    })
  }
  
  # Extract datetimes for all files
  datetimes <- sapply(tiff_files, extract_datetime)
  
  # Create dataframe with valid datetimes only
  raster_info_df <- tibble(
    filename = tiff_files,
    full_path = file.path(directory, tiff_files),
    start_datetime = datetimes
  ) %>% 
    filter(!is.na(start_datetime)) %>%
    arrange(start_datetime)
  
  # Check if any valid datetimes were found
  if (nrow(raster_info_df) == 0) {
    stop("No valid datetimes found in any of the files")
  }
  
  # Add formatted date columns
  raster_info_df <- raster_info_df %>%
    mutate(
      # Format as DD/MM/YYYY
      date_formatted = tryCatch({
        format(as.Date(start_datetime), "%d/%m/%Y")
      }, error = function(e) {
        warning("Error formatting date: ", e$message)
        NA_character_
      }),
      
      # Format as YYYY-MM-DD_HH:MM:SS
      datetime_formatted = tryCatch({
        paste0(
          strftime(start_datetime, "%Y-%m-%d", tz = "UTC"),
          "_",
          strftime(start_datetime, "%H:%M:%S", tz = "UTC")
        )
      }, error = function(e) {
        warning("Error formatting datetime: ", e$message)
        NA_character_
      })
    )
  
  message("Extracted dates from ", nrow(raster_info_df), " files.")
  return(raster_info_df)
}

# ==============================================
# Use extract_dates_from_filenames function to check the results
# ==============================================

# check the result:
data_dir = "Replace your input directory path here, where you save geotiff stack from step 1"
raster_info <- extract_dates_from_filenames(data_dir)
print(raster_info) 




# ==============================================
# 3. Load Raster Bands and Create Stacks
# ==============================================

load_raster_bands <- function(file_paths, datetime_formatted) {
  # Input validation
  if (length(file_paths) == 0) {
    stop("No file paths provided")
  }
  
  if (length(file_paths) != length(datetime_formatted)) {
    stop("Number of file paths (", length(file_paths), 
         ") does not match number of datetime strings (", 
         length(datetime_formatted), ")")
  }
  
  # Function to load individual file using terra
  load_single_file <- function(file_path) {
    # Check if file exists
    if (!file.exists(file_path)) {
      stop("File does not exist: ", file_path)
    }
    
    # Load raster data with error handling
    raster_data <- tryCatch({
      terra::rast(file_path)
    }, error = function(e) {
      stop("Error loading raster file: ", e$message)
    })
    
    # Check if file has at least 2 layers
    if (terra::nlyr(raster_data) < 2) {
      stop("File must have at least 2 layers: ", file_path)
    }
    
    # Return list with VV_dB and DpRVIc bands
    list(VV_dB = raster_data[[1]], DpRVIc = raster_data[[2]])
  }
  
  # Load all files
  message("Loading raster files...")
  all_bands <- lapply(file_paths, load_single_file)
  
  # Extract VV_dB and DpRVIc bands
  VV_rasters <- lapply(all_bands, function(x) x$VV_dB)
  DpRVIc_rasters <- lapply(all_bands, function(x) x$DpRVIc)
  
  # Check if all rasters have the same extent and resolution
  message("Checking raster extents and resolutions...")
  
  # Function to get extent and resolution info
  get_raster_info <- function(r) {
    list(
      extent = terra::ext(r),
      resolution = terra::res(r),
      crs = terra::crs(r)
    )
  }
  
  # Get info for all rasters
  VV_info <- lapply(VV_rasters, get_raster_info)
  DpRVIc_info <- lapply(DpRVIc_rasters, get_raster_info)
  
  # Function to check if all extents are the same
  all_same_extent <- function(extents) {
    if (length(extents) <= 1) return(TRUE)
    for (i in 2:length(extents)) {
      if (!all(extents[[1]] == extents[[i]])) return(FALSE)
    }
    return(TRUE)
  }
  
  # Function to check if all resolutions are the same
  all_same_resolution <- function(resolutions) {
    if (length(resolutions) <= 1) return(TRUE)
    for (i in 2:length(resolutions)) {
      if (!all(resolutions[[1]] == resolutions[[i]])) return(FALSE)
    }
    return(TRUE)
  }
  
  # Function to check if all CRS are the same
  all_same_crs <- function(crs_list) {
    if (length(crs_list) <= 1) return(TRUE)
    for (i in 2:length(crs_list)) {
      if (crs_list[[1]] != crs_list[[i]]) return(FALSE)
    }
    return(TRUE)
  }
  
  # Check if we need to resample
  need_resample <- !all_same_extent(lapply(VV_info, function(x) x$extent)) || 
    !all_same_resolution(lapply(VV_info, function(x) x$resolution)) ||
    !all_same_crs(lapply(VV_info, function(x) x$crs))
  
  if (need_resample) {
    message("Rasters have different extents, resolutions, or CRS. Resampling to a common reference...")
    
    # Find the common extent (intersection of all extents)
    common_extent <- VV_info[[1]]$extent
    for (i in 2:length(VV_info)) {
      common_extent <- terra::intersect(common_extent, VV_info[[i]]$extent)
    }
    
    # Use the resolution of the first raster as reference
    reference_resolution <- VV_info[[1]]$resolution
    reference_crs <- VV_info[[1]]$crs
    
    # Create a template raster with the common extent and reference resolution
    template <- terra::rast(common_extent, resolution = reference_resolution, crs = reference_crs)
    
    # Resample all rasters to the template
    message("Resampling VV backscatter rasters...")
    VV_rasters <- lapply(VV_rasters, function(x) {
      terra::resample(x, template, method = "bilinear")
    })
    
    message("Resampling DpRVIc rasters...")
    DpRVIc_rasters <- lapply(DpRVIc_rasters, function(x) {
      terra::resample(x, template, method = "bilinear")
    })
  } else {
    message("All rasters have the same extent, resolution, and CRS.")
  }
  
  # Create VV_dB stack using terra
  message("Creating VV backscatter stack...")
  VV_dB_stack <- tryCatch({
    terra::rast(VV_rasters)
  }, error = function(e) {
    message("Error creating VV_dB stack: ", e$message)
    message("Trying alternative stacking method...")
    
    # Alternative method: create a stack manually
    stack <- VV_rasters[[1]]
    for (i in 2:length(VV_rasters)) {
      stack <- terra::c(stack, VV_rasters[[i]])
    }
    return(stack)
  })
  
  # Set names for the layers
  names(VV_dB_stack) <- paste0("VV_", datetime_formatted)
  message("VV backscatter stack created with ", terra::nlyr(VV_dB_stack), " layers")
  
  # Create DpRVIc stack using terra
  message("Creating DpRVIc stack...")
  DpRVIc_stack <- tryCatch({
    terra::rast(DpRVIc_rasters)
  }, error = function(e) {
    message("Error creating DpRVIc stack: ", e$message)
    message("Trying alternative stacking method...")
    
    # Alternative method: create a stack manually
    stack <- DpRVIc_rasters[[1]]
    for (i in 2:length(DpRVIc_rasters)) {
      stack <- terra::c(stack, DpRVIc_rasters[[i]])
    }
    return(stack)
  })
  
  # Set names for the layers
  names(DpRVIc_stack) <- paste0("DpRVIc_", datetime_formatted)
  message("DpRVIc stack created with ", terra::nlyr(DpRVIc_stack), " layers")
  
  # Return both stacks
  list(VV_dB_stack = VV_dB_stack, DpRVIc_stack = DpRVIc_stack)
}


# ==============================================
# Use load_raster_bands function to check the results
# ==============================================

# Load raster bands and create stacks
message("Loading raster bands and creating stacks...")
stacks <- load_raster_bands(raster_info$full_path, raster_info$datetime_formatted)

# Access the stacks
VV_dB_stack <- stacks$VV_dB_stack
DpRVIc_stack <- stacks$DpRVIc_stack

# Print summary of the stacks
message("Summary of VV_dB_stack:")
print(summary(VV_dB_stack))

message("Summary of DpRVIc_stack:")
print(summary(DpRVIc_stack))

message("Processing completed successfully!") 




# ==============================================
# 4. Calculate Dry Reference Backscatter
# ==============================================

calculate_dry_reference <- function(VV_dB_stack, lower_quantile = 0.02, upper_quantile = 0.98) {
  message("Calculating dry reference backscatter...")
  
  # Input validation
  if (!inherits(VV_dB_stack, "SpatRaster")) {
    stop("Input must be a terra SpatRaster")
  }
  if (lower_quantile < 0 || lower_quantile > 1 || upper_quantile < 0 || upper_quantile > 1) {
    stop("Quantiles must be between 0 and 1")
  }
  if (lower_quantile >= upper_quantile) {
    stop("Lower quantile must be less than upper quantile")
  }
  
  # Function to calculate dry reference for a single pixel
  calc_dry_reference <- function(x) {
    # Remove NA values
    x_clean <- x[!is.na(x)]
    
    if (length(x_clean) >= 3) {
      # Remove top and bottom percentiles to exclude outliers
      q <- quantile(x_clean, probs = c(lower_quantile, upper_quantile), na.rm = TRUE)
      x_filtered <- x_clean[x_clean >= q[1] & x_clean <= q[2]]
      
      # Return minimum of filtered values, or minimum of all values if filtering removed everything
      return(ifelse(length(x_filtered) >= 1, min(x_filtered), min(x_clean)))
    } else if (length(x_clean) > 0) {
      # If we have fewer than 3 values, just take the minimum
      return(min(x_clean))
    } else {
      # No valid data
      return(NA)
    }
  }
  
  # Apply the function to each pixel across all layers
  message("Computing dry reference backscatter values...")
  sigma_dry <- terra::app(VV_dB_stack, fun = calc_dry_reference)
  
  # Set a descriptive name for the output raster
  names(sigma_dry) <- "dry_reference_backscatter"
  
  message("Dry reference backscatter calculation completed.")
  return(sigma_dry)
}

# ==============================================
# Use calculate_dry_reference function to check the results
# ==============================================

# Calculate dry reference backscatter
message("Calculating dry reference backscatter...")
sigma_dry <- calculate_dry_reference(VV_dB_stack)

# Print summary of the dry reference backscatter
message("Summary of dry reference backscatter values:")
print(summary(sigma_dry))

message("Processing completed successfully!") 




# ==============================================
# 5. Calculate backscatter change (Δσ)
# ==============================================
#' Calculate backscatter change (Δσ) relative to dry reference

calculate_delta_sigma <- function(VV_dB_stack, sigma_dry, raster_info_df) {
  message("Calculating backscatter change (Δσ)...")
  
  # Input validation
  if (!inherits(VV_dB_stack, "SpatRaster")) {
    stop("VV_dB_stack must be a terra SpatRaster")
  }
  if (!inherits(sigma_dry, "SpatRaster")) {
    stop("sigma_dry must be a terra SpatRaster")
  }
  if (!is.data.frame(raster_info_df) || !"date_formatted" %in% names(raster_info_df)) {
    stop("raster_info_df must be a data frame with a 'date_formatted' column")
  }
  if (terra::nlyr(VV_dB_stack) != nrow(raster_info_df)) {
    stop("Number of layers in VV_dB_stack (", terra::nlyr(VV_dB_stack), 
         ") does not match number of rows in raster_info_df (", 
         nrow(raster_info_df), ")")
  }
  
  # Function to calculate delta sigma for a single layer
  calc_delta_sigma_layer <- function(i) {
    # Get the current layer
    current_layer <- VV_dB_stack[[i]]
    
    # Calculate delta sigma (backscatter change relative to dry reference)
    delta_sigma <- current_layer - sigma_dry
    
    # Set negative values to 0 (physically implausible)
    # Using terra::app for pixel-wise operations
    delta_sigma <- terra::app(delta_sigma, fun = function(x) {
      ifelse(is.na(x), NA, pmax(0, x))
    })
    
    return(delta_sigma)
  }
  
  # Create a list to store individual delta sigma layers
  delta_sigma_layers <- list()
  
  # Calculate delta sigma for each layer
  message("Computing delta sigma for each layer...")
  for (i in 1:terra::nlyr(VV_dB_stack)) {
    message("Processing layer ", i, " of ", terra::nlyr(VV_dB_stack))
    delta_sigma_layers[[i]] <- calc_delta_sigma_layer(i)
  }
  
  # Combine all layers into a single SpatRaster
  message("Combining delta sigma layers into a stack...")
  delta_sigma_stack <- tryCatch({
    terra::rast(delta_sigma_layers)
  }, error = function(e) {
    message("Error creating delta sigma stack: ", e$message)
    message("Trying alternative stacking method...")
    
    # Alternative method: create a stack manually
    stack <- delta_sigma_layers[[1]]
    for (i in 2:length(delta_sigma_layers)) {
      stack <- terra::c(stack, delta_sigma_layers[[i]])
    }
    return(stack)
  })
  
  # Assign names based on dates
  names(delta_sigma_stack) <- paste0("Delta_sigma_", raster_info_df$date_formatted)
  
  message("Backscatter change calculation completed with ", terra::nlyr(delta_sigma_stack), " layers.")
  return(delta_sigma_stack)
}


# ==============================================
# Use calculate_delta_sigma function to check the results
# ==============================================

# Calculate delta sigma
message("Calculating delta sigma...")
delta_sigma_stack <- calculate_delta_sigma(VV_dB_stack, sigma_dry, raster_info)

# Print summary of the delta sigma stack
message("Summary of delta sigma values:")
print(summary(delta_sigma_stack))

message("Processing completed successfully!") 





# ==============================================
# 6. Calculate wet reference
# ==============================================
#' Calculate wet reference of backscatter change vs. DpRVIc

calculate_upper_envelope <- function(delta_sigma_stack, DpRVIc_stack, num_bins = 100, percentile = 0.98) { 
  # adjust number of bins according to the fitting of regression line
  message("Calculating upper envelope using binning approach...")
  
  # Input validation
  if (!inherits(delta_sigma_stack, "SpatRaster") || !inherits(DpRVIc_stack, "SpatRaster")) {
    stop("Inputs must be terra SpatRasters")
  }
  if (terra::nlyr(delta_sigma_stack) != terra::nlyr(DpRVIc_stack)) {
    stop("Stacks must have the same number of layers")
  }
  if (num_bins < 2) stop("Number of bins must be at least 2")
  if (percentile <= 0 || percentile >= 1) stop("Percentile must be between 0 and 1")
  
  # Combine data from all time steps for regression
  message("Extracting values from all layers...")
  combined_df <- do.call(rbind, lapply(1:terra::nlyr(delta_sigma_stack), function(i) {
    # Extract values from both rasters
    delta_sigma_values <- as.vector(terra::values(delta_sigma_stack[[i]]))
    DpRVIc_values <- as.vector(terra::values(DpRVIc_stack[[i]]))
    
    # Create a data frame with the values
    df <- data.frame(
      DpRVIc = DpRVIc_values,
      backscatterChange = delta_sigma_values
    )
    
    # Filter for DpRVIc >= 0.1 and remove NA values
    df <- df[!is.na(df$DpRVIc) & !is.na(df$backscatterChange) & df$DpRVIc >= 0.1,] 
    
    return(df)
  }))
  
  if (nrow(combined_df) == 0) {
    stop("No valid data points found after filtering")
  }
  
  message("Combined ", nrow(combined_df), " data points for envelope calculation.")
  
  # Create equal-width bins based on DpRVIc range
  DpRVIc_range <- range(combined_df$DpRVIc, na.rm = TRUE)
  bin_width <- diff(DpRVIc_range) / num_bins
  combined_df$DpRVIc_bin <- floor((combined_df$DpRVIc - DpRVIc_range[1]) / bin_width) + 1
  
  # Calculate percentile for each bin
  bin_thresholds <- aggregate(
    backscatterChange ~ DpRVIc_bin,
    data = combined_df,
    FUN = function(x) quantile(x, probs = percentile, na.rm = TRUE)
  )
  names(bin_thresholds)[2] <- "threshold"
  
  # Identify all points that are above or equal to their bin's percentile
  upper_envelope_df <- merge(combined_df, bin_thresholds, by = "DpRVIc_bin")
  upper_envelope_df <- upper_envelope_df[upper_envelope_df$backscatterChange >= upper_envelope_df$threshold, 
                                         c("DpRVIc", "backscatterChange")]
  
  # Sort by DpRVIc for visualization
  upper_envelope_df <- upper_envelope_df[order(upper_envelope_df$DpRVIc), ]
  
  message("Number of Upper Envelope Points: ", nrow(upper_envelope_df))
  
  return(list(
    upper_envelope_df = upper_envelope_df,
    combined_df = combined_df
  ))
}

# ==============================================
# Use calculate_upper_envelope function to check the results
# ==============================================

# Calculate upper envelope
message("Calculating upper envelope...")
envelope_results <- calculate_upper_envelope(delta_sigma_stack, DpRVIc_stack)

# Access the results
upper_envelope_df <- envelope_results$upper_envelope_df
combined_df <- envelope_results$combined_df

# Print summary of the upper envelope
message("Summary of upper envelope points:")
print(summary(upper_envelope_df))


message("Processing completed successfully!") 




# ==============================================
# 7. Fit Regression Model
# ==============================================
#' Fit regression model to upper envelope data

fit_regression_model <- function(upper_envelope_df, force_quadratic = FALSE) {
  message("Fitting regression models to upper envelope...")
  
  # Fit both linear and quadratic models
  linear_model <- lm(backscatterChange ~ DpRVIc, data = upper_envelope_df)
  quadratic_model <- lm(backscatterChange ~ poly(DpRVIc, 2, raw = TRUE), data = upper_envelope_df)
  
  # Compare model performance using multiple criteria
  linear_r2 <- summary(linear_model)$r.squared
  quadratic_r2 <- summary(quadratic_model)$r.squared
  
  # Calculate AIC and BIC for model comparison
  linear_aic <- AIC(linear_model)
  quadratic_aic <- AIC(quadratic_model)
  linear_bic <- BIC(linear_model)
  quadratic_bic <- BIC(quadratic_model)
  
  # Calculate RMSE
  linear_rmse <- sqrt(mean(residuals(linear_model)^2))
  quadratic_rmse <- sqrt(mean(residuals(quadratic_model)^2))
  
  # Print model comparison statistics
  message("\nModel Comparison Statistics:")
  message("Linear Model:")
  message("  R-squared: ", round(linear_r2, 4))
  message("  AIC: ", round(linear_aic, 2))
  message("  BIC: ", round(linear_bic, 2))
  message("  RMSE: ", round(linear_rmse, 4))
  
  message("\nQuadratic Model:")
  message("  R-squared: ", round(quadratic_r2, 4))
  message("  AIC: ", round(quadratic_aic, 2))
  message("  BIC: ", round(quadratic_bic, 2))
  message("  RMSE: ", round(quadratic_rmse, 4))
  
  # Select the better model based on multiple criteria
  if (force_quadratic) {
    message("\nForcing quadratic model as specified.")
    selected_model <- quadratic_model
    model_type <- "quadratic"
  } else {
    # Compare models using multiple criteria
    linear_score <- (linear_r2 * 0.4) + 
      (1 - (linear_aic - min(linear_aic, quadratic_aic)) / max(linear_aic, quadratic_aic)) * 0.3 +
      (1 - linear_rmse / max(linear_rmse, quadratic_rmse)) * 0.3
    
    quadratic_score <- (quadratic_r2 * 0.4) + 
      (1 - (quadratic_aic - min(linear_aic, quadratic_aic)) / max(linear_aic, quadratic_aic)) * 0.3 +
      (1 - quadratic_rmse / max(linear_rmse, quadratic_rmse)) * 0.3
    
    if (quadratic_score > linear_score) {
      message("\nSelecting quadratic model based on combined criteria.")
      selected_model <- quadratic_model
      model_type <- "quadratic"
    } else {
      message("\nSelecting linear model based on combined criteria.")
      selected_model <- linear_model
      model_type <- "linear"
    }
  }
  
  # Add predictions to the data frame
  upper_envelope_df$predicted_backscatterChange <- predict(selected_model, newdata = upper_envelope_df)
  
  # Generate equation text for plotting
  coefficients <- coef(selected_model)
  if (length(coefficients) == 2) {  # Linear
    intercept <- coefficients[1]
    slope <- coefficients[2]
    equation_text <- bquote(Y == .(round(slope, 2)) * x ~ 
                              .(ifelse(intercept < 0, "-", "+")) ~ 
                              .(abs(round(intercept, 2))))
  } else if (length(coefficients) == 3) {  # Quadratic
    a <- coefficients[3]
    b <- coefficients[2]
    c <- coefficients[1]
    equation_text <- bquote(Y == .(round(a, 2)) * x^2 ~ 
                              .(ifelse(b < 0, "-", "+")) ~ 
                              .(abs(round(b, 2))) * x ~ 
                              .(ifelse(c < 0, "-", "+")) ~ 
                              .(abs(round(c, 2))))
  }
  
  return(list(
    model = selected_model,
    model_type = model_type,
    upper_envelope_df = upper_envelope_df,
    equation_text = equation_text,
    coefficients = coefficients,
    model_stats = list(
      r2 = if(model_type == "linear") linear_r2 else quadratic_r2,
      aic = if(model_type == "linear") linear_aic else quadratic_aic,
      bic = if(model_type == "linear") linear_bic else quadratic_bic,
      rmse = if(model_type == "linear") linear_rmse else quadratic_rmse
    )
  ))
}


# ==============================================
# Use fit_regression_model function to check the results
# ==============================================

# Fit regression model
message("Fitting regression model to upper envelope...")
regression_results <- fit_regression_model(upper_envelope_df)

# Access the results
model <- regression_results$model
model_type <- regression_results$model_type
upper_envelope_df <- regression_results$upper_envelope_df
equation_text <- regression_results$equation_text
coefficients <- regression_results$coefficients
model_stats <- regression_results$model_stats

# Print summary of the regression model
message("Regression Model:")
message("Model Type: ", model_type)
message("Equation: ", equation_text)
message("Coefficients:")
print(coefficients)
message("Model Statistics:")
print(model_stats)

message("Processing completed successfully!") 





# ==============================================
# 8. Create Scatter Plot with Selected Model regression equation
# ==============================================

create_scatter_plot_with_model <- function(combined_df, upper_envelope_df, selected_model, model_type) {
  message("Creating scatter plot with selected model...")
  
  # Add predictions from the selected model
  upper_envelope_df$predicted_backscatterChange <- predict(selected_model, newdata = upper_envelope_df)
  
  # Generate equation text based on model type
  coefficients <- coef(selected_model)
  if (model_type == "linear") {
    intercept <- coefficients[1]
    slope <- coefficients[2]
    equation_text <- bquote(Y == .(round(slope, 2)) * x ~ 
                              .(ifelse(intercept < 0, "-", "+")) ~ 
                              .(abs(round(intercept, 2))))
  } else {
    a <- coefficients[3]
    b <- coefficients[2]
    c <- coefficients[1]
    equation_text <- bquote(Y == .(round(a, 2)) * x^2 ~ 
                              .(ifelse(b < 0, "-", "+")) ~ 
                              .(abs(round(b, 2))) * x ~ 
                              .(ifelse(c < 0, "-", "+")) ~ 
                              .(abs(round(c, 2))))
  }
  
  # Create the plot
  p <- ggplot() +
    # Plot all points in gray
    geom_point(data = combined_df, 
               aes(x = DpRVIc, y = backscatterChange), 
               color = "gray", size = 1, alpha = 0.6) +
    
    # Overlay the upper envelope points in red
    geom_point(data = upper_envelope_df, 
               aes(x = DpRVIc, y = backscatterChange), 
               color = "red", shape = 21, size = 2, stroke = 0.8) +
    
    # Add the regression line with confidence interval
    geom_smooth(data = upper_envelope_df, 
                aes(x = DpRVIc, y = backscatterChange),
                method = "lm", 
                formula = if(model_type == "quadratic") {y ~ poly(x, 2)} else {y ~ x},
                color = "white", 
                linewidth = 1,
                se = TRUE,
                alpha = 0.2) +
    
    # Add model equation
    annotate("text", 
             x = max(combined_df$DpRVIc, na.rm = TRUE) * 0.92, 
             y = max(combined_df$backscatterChange, na.rm = TRUE) * 0.92, 
             label = equation_text, 
             hjust = 1, 
             size = 5, 
             color = "white", 
             fontface = "bold") +
    
    # Set the axes limits and breaks
    scale_x_continuous(limits = c(0, 1), 
                       breaks = seq(0, 1, by = 0.2), 
                       expand = c(0, 0.01)) + 
    scale_y_continuous(limits = c(0, 16), 
                       breaks = seq(0, 16, by = 2), 
                       expand = c(0, 0)) +
    
    # Add axis labels
    labs(x = "DpRVIc", 
         y = expression(Delta * sigma ~ "(dB)")) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA),
      axis.line = element_line(color = "white", linewidth = 0.5),
      axis.ticks = element_line(color = "white", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title = element_text(color = "white", face = "bold", size = 12),
      axis.text = element_text(color = "white", face = "bold", size = 10)
    )
  
  return(p)
}



# ==============================================
# Use create_scatter_plot_with_model function to plot the scatter plot
# ==============================================

# Create scatter plot with selected model
message("Creating scatter plot with selected model...")
scatter_plot <- create_scatter_plot_with_model(combined_df, upper_envelope_df, model, model_type)
print(scatter_plot)

message("Processing completed successfully!") 


# Save the plot with high resolution, make sure to set your output directory
message("\nSaving plot...")
ggsave("FinalModelDeveloped.png", 
       plot = scatter_plot,
       width = 8, height = 5, 
       units = "in", 
       dpi = 600, 
       bg = "black")

message("\nAnalysis completed successfully!")
message("Results have been saved as 'FinalModelDeveloped.png'")



# ==============================================
# End of step 2, from here you can start your step 3
# ==============================================
