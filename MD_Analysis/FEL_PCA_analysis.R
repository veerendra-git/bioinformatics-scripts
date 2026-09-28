# FEL and PCA analysis of MD trajectory
# Input: PDB structure and DCD trajectory
# Output: RMSD, RMSF, PCA, clustering, DCCM and free energy plots

library(bio3d)
library(MASS)
library(fields)
library(plotly)
library(htmlwidgets)
library(webshot2)

# -----------------------------
# Input files and output folder
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

pdbfile <- if (length(args) >= 1) args[1] else "Ref.pdb"
dcdfile <- if (length(args) >= 2) args[2] else "Output.dcd"
outdir <- if (length(args) >= 3) args[3] else "MD_Results"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Constants
# -----------------------------

k_B_kj <- 0.0083144621
T <- 300
epsilon <- 1e-10

# -----------------------------
# Read PDB and DCD
# -----------------------------

pdb <- tryCatch(
  read.pdb(pdbfile),
  error = function(e) stop("Error reading PDB: ", e$message)
)

trj <- tryCatch(
  read.dcd(dcdfile),
  error = function(e) stop("Error reading DCD: ", e$message)
)

if (ncol(trj) != length(pdb$xyz)) {
  stop("Trajectory atom count does not match PDB atoms.")
}

# -----------------------------
# Select alpha carbon atoms
# -----------------------------

ca.inds <- atom.select(pdb, elety = "CA")

if (length(ca.inds$xyz) == 0) {
  stop("No alpha carbon atoms found.")
}

# -----------------------------
# Fit trajectory
# -----------------------------

xyz <- tryCatch(
  fit.xyz(
    fixed = pdb$xyz,
    mobile = trj,
    fixed.inds = ca.inds$xyz,
    mobile.inds = ca.inds$xyz
  ),
  error = function(e) stop("Trajectory fitting failed: ", e$message)
)

# -----------------------------
# RMSD
# -----------------------------

rd <- rmsd(
  xyz[1, ca.inds$xyz],
  xyz[, ca.inds$xyz]
)

png(
  file.path(outdir, "RMSD.png"),
  width = 2000,
  height = 1500,
  res = 300
)

plot(
  rd,
  type = "l",
  ylab = "RMSD (Å)",
  xlab = "Frame Number",
  main = "RMSD Over Time"
)

lines(
  lowess(rd),
  lty = 2,
  lwd = 2
)

dev.off()

# -----------------------------
# RMSF
# -----------------------------

rf <- rmsf(xyz[, ca.inds$xyz])

png(
  file.path(outdir, "RMSF.png"),
  width = 2000,
  height = 1500,
  res = 300
)

plot(
  rf,
  ylab = "RMSF (Å)",
  xlab = "Residue Index",
  type = "l",
  main = "RMSF per Residue"
)

dev.off()

# -----------------------------
# PCA
# -----------------------------

if (nrow(xyz) < 2) {
  stop("Not enough frames for PCA.")
}

pc <- tryCatch(
  pca.xyz(xyz[, ca.inds$xyz]),
  error = function(e) stop("PCA failed: ", e$message)
)

png(
  file.path(outdir, "PCA.png"),
  width = 2000,
  height = 1500,
  res = 300
)

plot(
  pc,
  col = rainbow(nrow(xyz)),
  main = "Principal Component Analysis"
)

dev.off()

pc_coords <- pc$z

# -----------------------------
# PCA clustering
# -----------------------------

if (ncol(pc_coords) >= 2) {

  k <- min(2, nrow(pc_coords))

  hc <- hclust(
    dist(pc_coords[, 1:2])
  )

  grps <- cutree(
    hc,
    k = k
  )

  png(
    file.path(outdir, "Clustered_PCA.png"),
    width = 2000,
    height = 1500,
    res = 300
  )

  plot(
    pc,
    col = grps,
    main = paste0("PCA Clustering (k=", k, ")")
  )

  dev.off()

} else {

  grps <- rep(
    1,
    nrow(pc_coords)
  )

  warning(
    "Less than two principal components available."
  )
}

# -----------------------------
# PC1 and PC2 displacement
# -----------------------------

if (ncol(pc$au) >= 2) {

  pc_disp_df <- data.frame(
    Residue_Index = seq_len(nrow(pc$au)),
    PC1_Displacement = pc$au[, 1],
    PC2_Displacement = pc$au[, 2]
  )

  write.csv(
    pc_disp_df,
    file.path(
      outdir,
      "Eigenvector_Displacement_PC1_PC2.csv"
    ),
    row.names = FALSE
  )

  png(
    file.path(outdir, "PC1_PC2_Displacement.png"),
    width = 2000,
    height = 1500,
    res = 300
  )

  plot(
    pc_disp_df$Residue_Index,
    pc_disp_df$PC1_Displacement,
    type = "l",
    lwd = 2,
    xlab = "Residue Index",
    ylab = "Displacement (Å)",
    main = "PC1 and PC2 Displacement"
  )

  lines(
    pc_disp_df$Residue_Index,
    pc_disp_df$PC2_Displacement,
    lwd = 2
  )

  legend(
    "topright",
    legend = c("PC1", "PC2"),
    lwd = 2
  )

  grid()

  dev.off()
}

# -----------------------------
# Dynamic cross-correlation matrix
# -----------------------------

cij <- tryCatch(
  dccm(xyz[, ca.inds$xyz]),
  error = function(e) {
    warning(
      "DCCM calculation failed: ",
      e$message
    )
    NULL
  }
)

if (!is.null(cij)) {

  png(
    file.path(outdir, "DCCM.png"),
    width = 2000,
    height = 1500,
    res = 300
  )

  plot(
    cij,
    main = "Dynamic Cross-Correlation Matrix"
  )

  dev.off()
}

# -----------------------------
# Free energy from RMSD
# -----------------------------

rmsd_nm <- rd / 10

dens_rmsd <- density(
  rmsd_nm,
  n = 200
)

free_energy <- -k_B_kj * T *
  log(dens_rmsd$y + epsilon)

free_energy <- free_energy -
  min(
    free_energy,
    na.rm = TRUE
  )

png(
  file.path(outdir, "FEL_RMSD.png"),
  width = 2000,
  height = 1500,
  res = 300
)

plot(
  dens_rmsd$x,
  free_energy,
  type = "l",
  lwd = 2,
  xlab = "RMSD (nm)",
  ylab = "Free Energy (kJ/mol)",
  main = "Free Energy Profile from RMSD"
)

dev.off()

# -----------------------------
# 2D PCA Free Energy Landscape
# -----------------------------

dens2d <- NULL
free_energy_2d <- NULL

if (ncol(pc_coords) >= 2) {

  dens2d <- tryCatch(
    kde2d(
      pc_coords[, 1],
      pc_coords[, 2],
      n = 200
    ),
    error = function(e) {
      warning(
        "2D KDE failed: ",
        e$message
      )
      NULL
    }
  )

  if (!is.null(dens2d)) {

    free_energy_2d <- -k_B_kj * T *
      log(dens2d$z + epsilon)

    free_energy_2d <- free_energy_2d -
      min(
        free_energy_2d,
        na.rm = TRUE
      )

    colorscale <- colorRampPalette(
      c(
        "blue",
        "cyan",
        "green",
        "yellow",
        "red"
      )
    )

    # -----------------------------
    # 2D FEL
    # -----------------------------

    png(
      file.path(outdir, "FEL_2D.png"),
      width = 2500,
      height = 2000,
      res = 300
    )

    par(
      mar = c(5, 5, 4, 6)
    )

    image.plot(
      x = dens2d$x,
      y = dens2d$y,
      z = free_energy_2d,
      col = colorscale(200),
      zlim = range(
        free_energy_2d,
        na.rm = TRUE
      ),
      xlab = "PC1",
      ylab = "PC2",
      main = "2D Free Energy Landscape",
      legend.lab = "Free Energy (kJ/mol)",
      legend.line = 4,
      legend.mar = 5
    )

    points(
      pc_coords[, 1],
      pc_coords[, 2],
      pch = 21,
      bg = ifelse(
        grps == 1,
        "black",
        "red"
      ),
      col = "black",
      cex = 1.2
    )

    dev.off()

    # -----------------------------
    # Zoomed 2D FEL
    # -----------------------------

    padding <- 0.5

    xlim_zoom <- range(
      pc_coords[, 1]
    ) + c(
      -padding,
      padding
    )

    ylim_zoom <- range(
      pc_coords[, 2]
    ) + c(
      -padding,
      padding
    )

    png(
      file.path(
        outdir,
        "FEL_2D_ZoomedIn.png"
      ),
      width = 2500,
      height = 2000,
      res = 300
    )

    par(
      mar = c(5, 5, 4, 6)
    )

    image.plot(
      x = dens2d$x,
      y = dens2d$y,
      z = free_energy_2d,
      col = colorscale(200),
      zlim = range(
        free_energy_2d,
        na.rm = TRUE
      ),
      xlim = xlim_zoom,
      ylim = ylim_zoom,
      xlab = "PC1",
      ylab = "PC2",
      main = "Zoomed-In 2D Free Energy Landscape",
      legend.lab = "Free Energy (kJ/mol)",
      legend.line = 4,
      legend.mar = 5
    )

    points(
      pc_coords[, 1],
      pc_coords[, 2],
      pch = 21,
      bg = ifelse(
        grps == 1,
        "black",
        "red"
      ),
      col = "black",
      cex = 1.2
    )

    dev.off()
  }
}

# -----------------------------
# 3D Free Energy Landscape
# -----------------------------

if (
  ncol(pc_coords) >= 2 &&
  !is.null(dens2d) &&
  !is.null(free_energy_2d)
) {

  energy_points <- fields::interp.surface(
    list(
      x = dens2d$x,
      y = dens2d$y,
      z = free_energy_2d
    ),
    loc = pc_coords[, 1:2]
  )

  clines <- contourLines(
    x = dens2d$x,
    y = dens2d$y,
    z = free_energy_2d,
    nlevels = 5
  )

  z_level <- min(
    free_energy_2d,
    na.rm = TRUE
  ) - 0.5

  fig <- plot_ly(
    x = dens2d$x,
    y = dens2d$y,
    z = free_energy_2d,
    type = "surface",
    colorscale = list(
      list(0, "blue"),
      list(0.25, "cyan"),
      list(0.5, "green"),
      list(0.75, "yellow"),
      list(1, "red")
    ),
    showscale = TRUE
  )

  fig <- fig %>%
    add_markers(
      x = pc_coords[, 2],
      y = pc_coords[, 1],
      z = energy_points,
      marker = list(
        color = ifelse(
          grps == 1,
          "black",
          "red"
        ),
        size = 4
      ),
      name = "Structures"
    )

  for (cline in clines) {

    fig <- fig %>%
      add_trace(
        x = cline$y,
        y = cline$x,
        z = rep(
          z_level,
          length(cline$x)
        ),
        type = "scatter3d",
        mode = "lines",
        showlegend = FALSE
      )
  }

  fig <- fig %>%
    layout(
      title = "3D Free Energy Landscape",
      scene = list(
        xaxis = list(
          title = "PC1"
        ),
        yaxis = list(
          title = "PC2"
        ),
        zaxis = list(
          title = "Free Energy (kJ/mol)"
        )
      )
    )

  htmlwidgets::saveWidget(
    fig,
    file.path(
      outdir,
      "3D_FEL_Funnel.html"
    ),
    selfcontained = TRUE
  )

  print(fig)
}

# -----------------------------
# Finished
# -----------------------------

cat(
  "\nMD analysis completed.\n"
)

cat(
  "Results saved in:",
  outdir,
  "\n"
)
