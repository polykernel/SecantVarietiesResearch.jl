{ inputs, config, pkgs, ... }:

let
  pkgs-unstable = import inputs.nixpkgs-unstable {
    inherit (pkgs.stdenv) system;
    config = {};
    overlays = [];
  };
in
{
  languages.julia = {
    enable = true;
    package = pkgs-unstable.julia-bin;
  };

  pre-commit.hooks = {
    treefmt = {
      enable = true;
    };
  };
}
