{ stdenv
, fetchFromGitHub
}:

stdenv.mkDerivation rec {
  pname = "bshoshany-thread-pool";
  version = "3.5.0";

  src = fetchFromGitHub {
    owner = "bshoshany";
    repo = "thread-pool";
    rev = "cabb3df5876c9a6824b07fcb0ff73d4a0e506ca0";
    sha256 = "7esTJEDh37kiuFUDO0ouyMQO0LJ7H43a0NtFB5KkbsY=";
  };

  dontBuild = true;

  installPhase = ''
    mkdir -p $out/include
    cp include/* $out/include
  '';
}
