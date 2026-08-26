$out = 'rd_extract.txt'
Remove-Item $out -ErrorAction SilentlyContinue
Get-ChildItem ..\man\*.Rd | Sort-Object Name | ForEach-Object {
  $c = Get-Content $_.FullName -Raw
  $t = ''
  if ($c -match '\\title\{([^}]*)\}') { $t = $Matches[1] }
  $u = ''
  if ($c -match '(?s)\\usage\{(.*?)\n\}') { $u = $Matches[1] }
  $u = $u -replace '\\dontrun\{', '' -replace '\\donttest\{', ''
  $u = $u -replace '\\link\{[^}]*\}', ''
  $lines = ($u.Trim() -split "`n") | Where-Object { $_.Trim() -ne '' } | Select-Object -First 8
  Add-Content $out ("=== " + $_.BaseName + " :: " + $t)
  Add-Content $out ($lines -join "`n")
}
Get-Content $out
