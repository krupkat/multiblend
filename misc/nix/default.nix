{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell rec {

  simde = pkgs.callPackage ./simde.nix {};
  thread-pool = pkgs.callPackage ./bshoshany-thread-pool.nix {};

  buildInputs = with pkgs; [
    ninja
    cmake
    libjpeg
    libpng
    libtiff
    catch2_3
    cereal
    spdlog
    thread-pool
    simde
    (python3.withPackages (pkgs: with pkgs; [ pillow ]))
  ];
}
