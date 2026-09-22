# Start the experiment's own CrystalPilot workbench instance from the frozen engine copy.
# Detached and hidden like scripts/restart_server.ps1 in the product, but with its own kernel home,
# its own port and its own log directory, so the user's production instance (port 8010) is untouched.
param([int]$Port = 8021, [int]$Cores = 4)
$exp = 'H:\CrystalPilot-campaigns\scxrd-agent-eval-20260921'
$root = Join-Path $exp 'engine\CrystalPilot'
$logdir = Join-Path $exp 'control\environment\server'
if (-not (Test-Path $logdir)) { New-Item -ItemType Directory -Path $logdir | Out-Null }
$pattern = 'uvicorn server\.app:app.*--port ' + $Port + '(\s|$)'
$existing = @(Get-CimInstance Win32_Process -Filter "Name='python.exe'" | Where-Object { $_.CommandLine -match $pattern })
if ($existing.Count -gt 0) { Write-Output ("already running: pids " + (($existing | ForEach-Object { $_.ProcessId }) -join ',')); exit 0 }
Set-Location $root
$env:PYTHONUTF8 = '1'
$env:CRYSTALPILOT_CPU_CORES = "$Cores"
$env:CRYSTALPILOT_CODEX_HOME = Join-Path $exp 'engine\codex-home-wb'
$python = Join-Path $root '.venv\Scripts\python.exe'
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$p = Start-Process -FilePath $python `
  -ArgumentList '-X','utf8','-m','uvicorn','server.app:app','--loop','asyncio:SelectorEventLoop','--host','127.0.0.1','--port',"$Port" `
  -WorkingDirectory $root `
  -RedirectStandardOutput (Join-Path $logdir "uvicorn_$stamp.log") `
  -RedirectStandardError (Join-Path $logdir "uvicorn_$stamp.err.log") `
  -WindowStyle Hidden -PassThru
Set-Content -Path (Join-Path $logdir 'uvicorn.pid') -Value ("launcher_pid={0} started={1} port={2} codex_home={3}" -f $p.Id, $stamp, $Port, $env:CRYSTALPILOT_CODEX_HOME)
$ok = $false
for ($i = 0; $i -lt 60; $i++) {
  Start-Sleep -Seconds 2
  try { $h = Invoke-WebRequest -UseBasicParsing ("http://127.0.0.1:{0}/api/health" -f $Port) -TimeoutSec 5; if ($h.StatusCode -eq 200) { $ok = $true; break } } catch { }
}
if ($ok) { Write-Output ("server up after ~{0}s" -f (($i + 1) * 2)); Write-Output $h.Content } else { Write-Output 'server NOT up after 120 s'; exit 1 }
