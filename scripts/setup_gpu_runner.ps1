# GitHub Self-Hosted Runner Quick Setup Script for Home Windows PC
# Instructions:
# 1. Open GitHub repo Settings -> Actions -> Runners -> New self-hosted runner (Windows)
# 2. Extract runner package into C:\actions-runner
# 3. Run config.cmd with your token, adding label "gpu"
# 4. Run `.\run.cmd` or install as service (`.\svc.sh install; .\svc.sh start`)

Write-Host "GitHub Self-Hosted Runner configured for GPU Benchmark." -ForegroundColor Green
Write-Host "Label required: gpu, windows, self-hosted" -ForegroundColor Yellow
