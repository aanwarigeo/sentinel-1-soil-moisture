################################################################################
# Step 3: Soil Moisture Estimation
################################################################################
# 
# This script outlines a step-by-step workflow for step 3 to estimate field-scale
# soil moisture using sentinel-1 GRD SAR data.
#
# Key steps in this workflow:
# 1. Calculate relative soil moisture as ratio of backscatter change to maximum
# 2. Convert relative soil moisture to volumeteric soil moisture using ground 
# measurements

# For a detailed explanation of each step of the methodology, please refer to the 
# research paper by Bhogapurapu et al. (2022).

# Author: aanwari
# Date: Feb 28, 2025

################################################################################


# ==============================================
# 10. Load Required Libraries
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
# 11. Calculate Relative Soil Moisture (Θ)
# ==============================================
#' Calculate relative soil moisture using the change detection approach 
calculate_relative_soil_moisture <- function(delta_sigma_stack, DpRVIc_stack, regression_model, model_type) {
  message("Calculating relative soil moisture (Θ)...")
  
  # Input validation
  if (!inherits(delta_sigma_stack, "SpatRaster") || !inherits(DpRVIc_stack, "SpatRaster")) {
    stop("Inputs must be terra SpatRasters")
  }
  if (terra::nlyr(delta_sigma_stack) != terra::nlyr(DpRVIc_stack)) {
    stop("Stacks must have the same number of layers")
  }
  
  # Create a function to predict maximum backscatter change for each DpRVIc value
  predict_max_backscatter <- function(x) {
    # Create a data frame with DpRVIc values
    df <- data.frame(DpRVIc = x)
    
    # Predict using the appropriate model formula
    if (model_type == "quadratic") {
      df$DpRVIc_squared <- df$DpRVIc^2
      predicted <- predict(regression_model, newdata = df)
    } else {
      predicted <- predict(regression_model, newdata = df)
    }
    
    return(predicted)
  }
  
  # Calculate predicted maximum backscatter change for each layer
  message("Calculating predicted maximum backscatter change...")
  predicted_sigma_max_stack <- terra::app(DpRVIc_stack, fun = predict_max_backscatter)
  
  # Calculate relative soil moisture (Θ) for each layer
  message("Calculating relative soil moisture...")
  Theta_layers <- list()
  
  for (i in 1:terra::nlyr(delta_sigma_stack)) {
    message("Processing layer ", i, " of ", terra::nlyr(delta_sigma_stack))
    
    # Get current layers
    current_delta_sigma <- delta_sigma_stack[[i]]
    current_predicted_max <- predicted_sigma_max_stack[[i]]
    
    # Extract values as vectors
    delta_sigma_values <- as.vector(terra::values(current_delta_sigma))
    predicted_max_values <- as.vector(terra::values(current_predicted_max))
    
    # Ensure both vectors have the same length
    if (length(delta_sigma_values) != length(predicted_max_values)) {
      stop("Non-conformable arrays: delta_sigma and predicted_max have different lengths")
    }
    
    # Calculate theta values
    theta_values <- rep(NA, length(delta_sigma_values))
    
    # Process in chunks to avoid memory issues
    chunk_size <- 1000000
    num_chunks <- ceiling(length(delta_sigma_values) / chunk_size)
    
    for (j in 1:num_chunks) {
      start_idx <- (j-1) * chunk_size + 1
      end_idx <- min(j * chunk_size, length(delta_sigma_values))
      
      # Get current chunk
      delta_chunk <- delta_sigma_values[start_idx:end_idx]
      pred_chunk <- predicted_max_values[start_idx:end_idx]
      
      # Calculate ratio (equation 10)
      theta_chunk <- delta_chunk / pred_chunk
      
      # Handle NA values
      theta_chunk[is.na(delta_chunk) | is.na(pred_chunk)] <- NA
      
      # Cap at 1 (fully saturated)
      theta_chunk <- pmin(theta_chunk, 1)
      
      # Set negative values to 0 (physically implausible)
      theta_chunk <- pmax(theta_chunk, 0)
      
      # Store in result vector
      theta_values[start_idx:end_idx] <- theta_chunk
    }
    
    # Create a new SpatRaster with the same properties as the input
    theta_raster <- terra::rast(
      nrows = terra::nrow(current_delta_sigma),
      ncols = terra::ncol(current_delta_sigma),
      xmin = terra::xmin(current_delta_sigma),
      xmax = terra::xmax(current_delta_sigma),
      ymin = terra::ymin(current_delta_sigma),
      ymax = terra::ymax(current_delta_sigma),
      crs = terra::crs(current_delta_sigma)
    )
    
    # Set values
    terra::values(theta_raster) <- theta_values
    
    # Add to list
    Theta_layers[[i]] <- theta_raster
  }
  
  # Combine all layers into a single SpatRaster
  message("Combining layers into a single SpatRaster...")
  Theta_stack <- terra::rast(Theta_layers)
  
  # Assign names to the layers
  names(Theta_stack) <- paste0("Theta_", names(delta_sigma_stack))
  
  message("Relative soil moisture calculation completed.")
  return(list(
    Theta_stack = Theta_stack,
    predicted_sigma_max_stack = predicted_sigma_max_stack
  ))
}


# ==============================================
# Run calculate_relative_soil_moisture function 
# ==============================================

# Calculate relative soil moisture
message("Calculating relative soil moisture...")
soil_moisture_results <- calculate_relative_soil_moisture(delta_sigma_stack, DpRVIc_stack, model, model_type)

# Access the results
Theta_stack <- soil_moisture_results$Theta_stack
predicted_sigma_max_stack <- soil_moisture_results$predicted_sigma_max_stack

# Print summary of the relative soil moisture stack
message("Summary of relative soil moisture values:")
print(summary(Theta_stack))


message("Processing completed successfully!") 



# ==============================================
# 12. Calculate Volumeteric Soil Moisture (Θv)
# ==============================================
#' Calculate volumeteric soil moisture using ground-based wilting point/field capacity values

calculate_absolute_soil_moisture <- function(Theta_stack, min_ground_sm = 0.176, max_ground_sm = 0.383, raster_info_df) {
  message("Calculating absolute soil moisture (θ)...")
  
  # Input validation
  if (!inherits(Theta_stack, "SpatRaster")) {
    stop("Input must be a terra SpatRaster")
  }
  
  if (min_ground_sm >= max_ground_sm) {
    stop("Minimum ground soil moisture must be less than maximum")
  }
  
  # Calculate absolute soil moisture using equation: θv = θ_min + Θ * (θ_max - θ_min)
  message("Converting relative to absolute soil moisture...")
  
  # Use vectorized operations instead of function for better performance
  range_diff <- max_ground_sm - min_ground_sm
  
  # Apply calculation to the entire stack at once
  absolute_sm_stack <- min_ground_sm + Theta_stack * range_diff
  
  # Ensure values are within physical bounds 
  absolute_sm_stack <- terra::clamp(absolute_sm_stack, min_ground_sm, max_ground_sm)
  
  # Assign appropriate names to the layers
  n_layers <- terra::nlyr(absolute_sm_stack)
  
  if (!is.null(raster_info_df) && !is.data.frame(raster_info_df)) {
    warning("raster_info_df must be a data frame. Using default layer names.")
    names(absolute_sm_stack) <- paste0("absolute_soil_moisture_", 1:n_layers)
  } else if (is.null(raster_info_df) || n_layers != nrow(raster_info_df)) {
    warning("Number of layers (", n_layers, 
            ") does not match number of rows in raster_info_df (", 
            ifelse(is.null(raster_info_df), "NULL", nrow(raster_info_df)), ")")
    names(absolute_sm_stack) <- paste0("absolute_soil_moisture_", 1:n_layers)
  } else {
    # Determine which date column to use
    date_col <- if ("date_formatted" %in% names(raster_info_df)) "date_formatted" else "datetime_formatted"
    
    if (date_col %in% names(raster_info_df)) {
      names(absolute_sm_stack) <- paste0("absolute_soil_moisture_", raster_info_df[[date_col]])
    } else {
      names(absolute_sm_stack) <- paste0("absolute_soil_moisture_", 1:n_layers)
      warning("Neither 'date_formatted' nor 'datetime_formatted' found in raster_info_df")
    }
  }
  
  message("Absolute soil moisture calculation completed.")
  return(absolute_sm_stack)
}

# ==============================================
# 13. Save Individual Soil Moisture Layers
# ==============================================

save_individual_layers <- function(raster_stack, output_dir, raster_info_df = NULL, 
                                   prefix = "absolute_soil_moisture") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    message("Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get the number of layers
  num_layers <- terra::nlyr(raster_stack)
  message("Saving ", num_layers, " individual layers...")
  
  # Prepare vector to store filenames
  saved_files <- character(num_layers)
  
  # Determine date column if raster_info_df is provided
  date_col <- NULL
  if (!is.null(raster_info_df) && is.data.frame(raster_info_df)) {
    if ("datetime_formatted" %in% names(raster_info_df)) {
      date_col <- "datetime_formatted"
    } else if ("date_formatted" %in% names(raster_info_df)) {
      date_col <- "date_formatted"
    }
  }
  
  # Save each layer
  for (i in 1:num_layers) {
    # Create filename based on datetime
    if (!is.null(date_col) && i <= nrow(raster_info_df)) {
      datetime_str <- raster_info_df[[date_col]][i]
      # Replace problematic characters for filesystem compatibility
      safe_datetime_str <- gsub("[:/]", "-", datetime_str)
      filename <- file.path(output_dir, paste0(prefix, "_", safe_datetime_str, ".tif"))
    } else {
      # Fallback if no date info available
      filename <- file.path(output_dir, paste0(prefix, "_layer_", i, ".tif"))
    }
    
    # Save the layer
    message("Saving layer ", i, " to '", filename, "'")
    terra::writeRaster(raster_stack[[i]], filename, overwrite = TRUE)
    
    saved_files[i] <- filename
  }
  
  message("All individual layers saved to: ", output_dir)
  return(invisible(saved_files))
}

# ==============================================
# Process and Save Soil Moisture Data
# ==============================================

process_soil_moisture <- function(Theta_stack, raster_info, output_dir = "output") {
  # Calculate absolute soil moisture
  message("Calculating absolute soil moisture...")
  absolute_soil_moisture_stack <- calculate_absolute_soil_moisture(
    Theta_stack, 
    raster_info_df = raster_info
  )
  
  # Print summary of the absolute soil moisture stack
  message("Summary of absolute soil moisture values:")
  print(summary(absolute_soil_moisture_stack))
  
  # Save individual layers
  saved_files <- save_individual_layers(
    absolute_soil_moisture_stack, 
    output_dir, 
    raster_info_df = raster_info
  )
  
  message("Processing completed successfully!")
  return(invisible(absolute_soil_moisture_stack))
}


# use process_soil_moisture to process and save absolute soil moisture to output directory
message("Processing and saving absolute soil moisture...")
absolute_soil_moisture_stack <- process_soil_moisture(
  Theta_stack, 
  raster_info, 
  # Replace the path below with your own output folder location
  "/Users/replace your own output folder"
)




# ==============================================
# 14. Extract Soil Moisture at Specific Coordinates
# ==============================================
# Replace the path below with your input folder where you save raster stack from previous step
extract_soil_moisture_at_coords <- function(output_dir = "replace your own output folder",
                                            target_coords,
                                            file_pattern = "*.tif",
                                            save_excel = TRUE) {
  message("Extracting soil moisture values at coordinates: (", 
          target_coords[1], ", ", target_coords[2], ")")
  
  # Input validation
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist: ", output_dir)
  }
  
  if (length(target_coords) != 2) {
    stop("target_coords must be a vector of length 2 (longitude, latitude)")
  }
  
  # Create a point from coordinates (using sf)
  point <- sf::st_point(target_coords)
  point_sf <- sf::st_sfc(point, crs = "EPSG:4326")
  
  # Get list of GeoTIFF files
  tif_files <- list.files(output_dir, pattern = file_pattern, full.names = TRUE)
  
  # Filter out any non-GeoTIFF files
  tif_files <- tif_files[grepl("\\.tif$", tif_files, ignore.case = TRUE)]
  
  if (length(tif_files) == 0) {
    stop("No GeoTIFF files found in directory: ", output_dir)
  }
  
  message("Found ", length(tif_files), " GeoTIFF files")
  
  # Process each file and collect results using purrr-like approach with dplyr
  results <- dplyr::tibble(
    file_path = tif_files
  ) %>%
    dplyr::mutate(
      # Extract filename
      filename = basename(file_path),
      # Extract datetime from filename (remove prefix and extension)
      datetime_str = gsub("^.*_(\\d{4}-\\d{2}-\\d{2}.*)\\.tif$", "\\1", filename),
      # Format datetime string for proper parsing
      datetime_str_cleaned = gsub("-(?=\\d{2}-\\d{2}$)", " ", datetime_str, perl = TRUE),
      datetime_str_formatted = gsub("-", ":", datetime_str_cleaned)
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # Read raster and extract value
      raster_value = {
        message("Processing file: ", filename)
        # Use terra instead of raster package
        r <- terra::rast(file_path)
        # Convert sf point to terra's expected format if needed
        if (terra::crs(r) != sf::st_crs(point_sf)$wkt) {
          point_transformed <- sf::st_transform(point_sf, terra::crs(r))
          v <- terra::extract(r, terra::vect(point_transformed))[[2]]
        } else {
          v <- terra::extract(r, terra::vect(point_sf))[[2]]
        }
        v
      }
    ) %>%
    dplyr::ungroup() %>%
    # Select and rename columns
    dplyr::transmute(
      absoluteSoilMoistureDateTime = datetime_str_formatted,
      absoluteSoilMoisture = as.numeric(raster_value)
    )
  
  # Use lubridate for proper datetime handling
  results <- results %>%
    dplyr::mutate(
      parsed_datetime = lubridate::parse_date_time(
        absoluteSoilMoistureDateTime, 
        orders = c("Y-m-d H:M:S", "Y-m-d", "Y-m-d_H:M:S")
      ),
      # Format datetime consistently
      absoluteSoilMoistureDateTime = format(parsed_datetime, "%Y-%m-%d %H:%M:%S")
    ) %>%
    # Sort by datetime
    dplyr::arrange(parsed_datetime) %>%
    # Remove the helper column
    dplyr::select(-parsed_datetime)
  
  # Save to Excel if requested
  if (save_excel) {
    output_file <- file.path(output_dir, "estimatedSoilMoisture.xlsx")
    message("\nSaving extracted values to: ", output_file)
    
    # Create workbook and add worksheet
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Soil Moisture")
    openxlsx::writeData(wb, "Soil Moisture", results)
    
    # Add formatting
    header_style <- openxlsx::createStyle(
      textDecoration = "bold",
      halign = "center",
      border = "bottom"
    )
    openxlsx::addStyle(wb, "Soil Moisture", style = header_style, rows = 1, cols = 1:2)
    
    # Save workbook
    openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  }
  
  # Print summary
  message("\nSummary of extracted values:")
  message("Total files processed: ", nrow(results))
  message("NA values: ", sum(is.na(results$absoluteSoilMoisture)))
  
  # Handle case when all values are NA
  if (all(is.na(results$absoluteSoilMoisture))) {
    message("All extracted values are NA. Check coordinate system compatibility.")
  } else {
    message("Range of soil moisture values: ", 
            round(min(results$absoluteSoilMoisture, na.rm = TRUE), 4), " to ", 
            round(max(results$absoluteSoilMoisture, na.rm = TRUE), 4))
  }
  
  # Print the first few rows
  message("\nFirst few rows of the extracted data:")
  print(head(results))
  
  return(results)
}



# ==============================================
# Run extract_soil_moisture_at_coords function to extract soil moisture values at specific coordinates
# ==============================================

# Extract soil moisture values at specific coordinates
message("\nExtracting soil moisture values at specific coordinates...")
# Define your own target coordinates here
target_coords <- c(-7.0211013, 37.3819269)
extracted_values <- extract_soil_moisture_at_coords(
  output_dir = "replace your own output folder",
  target_coords = target_coords,
  save_excel = TRUE
)
message("Extraction completed successfully!") 
