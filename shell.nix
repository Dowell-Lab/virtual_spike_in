{ pkgs ? import <unstable> { overlays = [ (import ./overlay.nix) ]; } }:
let
	# Use clang
	stdenv = pkgs.gccStdenv;
  # Allow broken packages
  allowBroken = true;
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
in pkgs.mkShell {
  name = "virtualSpikeInEnv";
  buildInputs = [
		# Generic dependencies
		pkgs.git
		# C/C++ Dependencies
    pkgs.gcc
    pkgs.bedtools
    pkgs.blas
    # Library Dependencies
    pkgs.netcdf
    pkgs.graphviz
    hdf5
	  # Python Dependencies
    pkgs.python39
    pkgs.poetry
    pkgs.python39Packages.python-lsp-server
		# R Dependencies
		pkgs.R
    # System dependencies
    frameworks.Accelerate
  ];
  shellHook = ''
    # Update library paths
    export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:${pkgs.graphviz}/lib/pkgconfig
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.stdenv.cc.cc.lib}/lib
    export HDF5_DIR=${hdf5}
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.hdf5}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.netcdf}/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${pkgs.graphviz}/lib
  '';

}
