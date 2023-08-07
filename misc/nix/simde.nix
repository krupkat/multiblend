{ stdenv
, fetchFromGitHub
}:

stdenv.mkDerivation rec {
  pname = "simde";
  version = "0.7.6";

  src = fetchFromGitHub {
    owner = "simd-everywhere";
    repo = "simde";
    rev = "fefc7857ff3e785b988a61a8f5f3c5bd5eb24342";
    sha256 = "pj+zaD5o9XYkTavezcQFzM6ao0IdQP1zjP9L4vcCyEY=";
  };

  dontBuild = true;

  installPhase = ''
    mkdir -p $out/include
    cp -r simde $out/include/
  '';
}
