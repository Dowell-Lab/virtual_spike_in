{ pkgs ? import <unstable> { overlays = [ (import ./overlay.nix) ]; } }:
let
  # Allow broken packages
  # allowBroken = true;
  # Use darwin frameworks
  frameworks = pkgs.darwin.apple_sdk.frameworks;
  # Fix hdf5 output
  hdf5 = pkgs.symlinkJoin {
    name = "hdf5";
    paths = [
      pkgs.hdf5
      pkgs.hdf5.dev
    ];
  };
  geneplotter = pkgs.rPackages.geneplotter.overrideAttrs (old: {
    buildInputs = old.buildInputs ++ [ pkgs.libiconv ];
  });
  genefilter = pkgs.rPackages.genefilter.overrideAttrs (old: {
    buildInputs = old.buildInputs ++ [ pkgs.libiconv ];
  });
  # Specify R packages
  my-r-pkgs = pkgs.rWrapper.override {
    packages = [
      pkgs.rPackages.tidyverse
      pkgs.rPackages.extrafont
      pkgs.rPackages.plyr
      pkgs.rPackages.ggthemes
      pkgs.rPackages.KernSmooth
      pkgs.rPackages.BiocManager
      pkgs.rPackages.Biobase
      pkgs.rPackages.BiocGenerics
      pkgs.rPackages.BiocParallel
      pkgs.rPackages.DESeq2
      pkgs.rPackages.scico
      genefilter
      geneplotter
      pkgs.rPackages.GenomicRanges
      pkgs.rPackages.ggplot2
      pkgs.rPackages.IRanges
      pkgs.rPackages.locfit
      pkgs.rPackages.Rcpp
      pkgs.rPackages.RcppArmadillo
      pkgs.rPackages.S4Vectors
      pkgs.rPackages.SummarizedExperiment
      pkgs.rPackages.styler
      pkgs.rPackages.patchwork
      pkgs.rPackages.languageserver
      # DESeq2
    ];
  };
in pkgs.mkShell {
  name = "virtualSpikeInEnv";
  buildInputs = [
    # Generic dependencies
    pkgs.git
    # C/C++ Dependencies
    pkgs.gcc
    pkgs.gfortran
    pkgs.bedtools
    pkgs.blas
    pkgs.openblas
    pkgs.gettext
    pkgs.libiconv
    pkgs.libpng
    # Library Dependencies
    pkgs.netcdf # override
    pkgs.graphviz
    hdf5
    # Python Dependencies
    pkgs.python39
    pkgs.poetry
    # pkgs.python39Packages.python-lsp-server
    # R dependencies
    my-r-pkgs
    # System dependencies
    frameworks.Accelerate
    frameworks.CoreFoundation
  ];
  shellHook = ''
    # Update library paths
    export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:${pkgs.graphviz}/lib/pkgconfig
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.stdenv.cc.cc.lib}/lib
    export HDF5_DIR=${hdf5}
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.libiconv}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.libpng}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.hdf5}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.netcdf}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.graphviz}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.gettext}/lib
    mkdir -p "$(pwd)/_libs"
    export R_LIBS_USER="$(pwd)/_libs"
    Rscript -e "BiocManager::install('DESeq2', force=TRUE)"
  '';

}
