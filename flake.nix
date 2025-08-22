{
  description = "Flake for a R environment";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs { inherit system; };
        breedSimulatR_buildInputs = with pkgs.rPackages; [
          data_table
          R6
          vcfR
          covr
          gaston
          knitr
          pkgdown
          plotly
          rmarkdown
          roxygen2
          spelling
          testthat
        ];
        R-packages =
          with pkgs.rPackages;
          [
            # developement packages:
            devtools
            roxygen2
            roxygen2md
            usethis
            languageserver
            styler
            tidyverse
            ggplot2
          ]
          ++ breedSimulatR_buildInputs;
        breedSimulatR = (
          pkgs.rPackages.buildRPackage {
            name = "breedSimulatR";
            src = ./.;
            propagatedBuildInputs = breedSimulatR_buildInputs;
            doCheck = false;
          }
        );
        R-with-packages = pkgs.rWrapper.override { packages = R-packages; };
        Rstudio-with-packages = pkgs.rstudioWrapper.override { packages = R-packages; };
      in
      rec {
        devShells.default = pkgs.mkShell {
          LOCALE_ARCHIVE =
            if "${system}" == "x86_64-linux" then "${pkgs.glibcLocalesUtf8}/lib/locale/locale-archive" else "";
          R_LIBS_USER = "''"; # to no use users' installed R packages
          R_PROFILE_USER = "''"; # to disable`.Rprofile` files (eg. when the project already use `renv`)
          nativeBuildInputs = [ pkgs.bashInteractive ];
          buildInputs = [
            R-with-packages
          ]
          ++ pkgs.lib.optionals (pkgs.stdenv.isLinux) [
            Rstudio-with-packages
          ];
        };

        apps = {
        };
      }
    );
}
