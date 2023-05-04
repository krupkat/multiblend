$multiblend_sources = @(Get-ChildItem -Recurse -Path src/ -Include *.cpp,*.h).fullname

clang-format -i @($multiblend_sources)
