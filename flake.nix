{
  description = "rhosa";

  inputs.nixpkgs.url = github:NixOS/nixpkgs/nixos-unstable;

  outputs = { self, nixpkgs }: let

    allSystems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];

    forAllSystems = f: nixpkgs.lib.genAttrs allSystems (system: f (import nixpkgs { inherit system; }));

    rhosa-dev = pkgs: with pkgs; let

      custom-R = rWrapper.override {
        packages = with rPackages; [
          devtools
          ggplot2
          knitr
          multitaper
          rmarkdown
          testthat
          tidyverse
        ];
      };

      rhosa-dev-shell = mkShell {
        packages = [
          custom-R
          gettext
          qpdf
        ];
      };

    in {

      default = rhosa-dev-shell;

    };

  in {

    packages = forAllSystems rhosa-dev;

    devShells = forAllSystems rhosa-dev;

  };
}
