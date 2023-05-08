$tidy_runner = (Get-Command run-clang-tidy).Path

python $tidy_runner -p .\build\Release -export-fixes fixes.yaml
python .\misc\python\clang_tidy_deduplicate.py fixes.yaml
clang-apply-replacements .

.\misc\format.ps1
