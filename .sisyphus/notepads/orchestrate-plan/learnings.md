## [2025-02-28 00:00:00 UTC] Task: notebooklm-auth-setup

### Context
- Working directory: `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline`
- NotebookLM skill installed as user-priority skill
- Task: Prepare safe, exact auth setup procedure without performing interactive login in HPC environment

### Key Findings

#### 1. NotebookLM Skill Architecture
- **Scripts location**: Skill is installed separately (not in repo), accessed via dedicated environment
- **Entry point**: `python scripts/run.py` wrapper (standard pattern across all skill operations)
- **Auth file location**: `~/.claude/skills/notebooklm/data/auth_info.json`
- **Env setup**: Automatic via `run.py` - creates `.venv`, installs deps, activates environment
- **Browser requirement**: MUST be visible/interactive for manual Google login (headless = fails)

#### 2. Auth Manager Commands (Official Skill API)
```bash
# Check authentication status (non-interactive, safe on HPC)
python scripts/run.py auth_manager.py status

# One-time setup with interactive Google login (USER MUST RUN LOCALLY)
python scripts/run.py auth_manager.py setup

# Re-authenticate if token expires
python scripts/run.py auth_manager.py reauth

# Clear all auth data (for troubleshooting or switching accounts)
python scripts/run.py auth_manager.py clear
```

#### 3. Expected Status Output (Unauthenticated)
When running `status` on fresh install:
- Returns: `Not authenticated` or similar status indicator
- No auth_info.json exists yet
- Safe to run repeatedly (read-only operation)

#### 4. Setup Flow (User-Interactive on Local Machine)
1. User runs: `python scripts/run.py auth_manager.py setup`
2. Browser opens automatically (MUST NOT be headless)
3. User manually logs into Google account
4. User grants NotebookLM API permissions
5. Browser redirects to success page / closes
6. Auth token saved to: `~/.claude/skills/notebooklm/data/auth_info.json`
7. User verifies with: `python scripts/run.py auth_manager.py status`

#### 5. Critical Constraints (HPC Environment)
- **Cannot run interactive Google login on HPC** (headless/non-interactive environment)
- **Solution**: User must run `setup` on their local machine with visible browser
- **Auth persists**: Once authenticated locally, token works in HPC environment for API calls
- **Rate limits**: Google free tier = 50 queries/day per account
