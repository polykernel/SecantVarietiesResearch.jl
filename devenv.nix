{ config, ... }:

{
  languages.julia.enable = true;

  pre-commit.hooks = {
    treefmt = {
      enable = true;
    };
  };
}