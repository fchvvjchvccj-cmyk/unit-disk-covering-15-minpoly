# 多代理共用计算租约，保证构建、测试和性能试算串行执行。
param([Parameter(Mandatory=$true, Position=0)][string]$Program,
      [Parameter(ValueFromRemainingArguments=$true)][string[]]$ProgramArguments)
$ErrorActionPreference = 'Stop'
$computeMutex = [System.Threading.Mutex]::new($false, 'Local\CodexSxypComputeLease')
$computeAcquired = $false
try {
    try { $computeAcquired = $computeMutex.WaitOne([TimeSpan]::FromMinutes(40)) }
    catch [System.Threading.AbandonedMutexException] { $computeAcquired = $true }
    if (-not $computeAcquired) { throw '等待计算租约超过40分钟' }
    & $Program @ProgramArguments
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} finally {
    if ($computeAcquired) { $computeMutex.ReleaseMutex() }
    $computeMutex.Dispose()
}
