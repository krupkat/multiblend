Param
(
    [Parameter(Position=0)]
    $BuildDirParam = '.\build\Release',
    [Parameter(Position=1)]
    $ExportFixesParam = 'fixes.yaml'
)

$tidy_runner = (Get-Command run-clang-tidy).Path
python $tidy_runner src -p $BuildDirParam -quiet -export-fixes $ExportFixesParam
