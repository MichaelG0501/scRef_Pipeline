## [2025-02-28 00:00:00 UTC] Task: notebooklm-auth-setup

### Recommended Procedure

#### Phase 1: Verification (Can Run on HPC or Local)
Use this non-interactive command to check if authentication is needed:
```bash
python scripts/run.py auth_manager.py status
```

**Expected outputs:**
- If authenticated: Reports account email and token validity
- If unauthenticated: Reports `Not authenticated` or empty state

---

#### Phase 2: One-Time Setup (MUST RUN ON LOCAL MACHINE)
**Prerequisites:**
- ✅ Local machine with visible desktop/browser
- ✅ Python environment available locally
- ✅ NotebookLM skill installed in your Claude environment
- ⚠️  NOT on HPC node (headless environment will fail)

**Step-by-step:**

1. **On your local machine**, open terminal/command prompt
2. **Run setup command:**
   ```bash
   python scripts/run.py auth_manager.py setup
   ```
3. **Browser window opens automatically** (do NOT close it)
4. **Log in to Google** with the account you want to use for NotebookLM
5. **Grant permissions** when NotebookLM requests API access
6. **Wait for success page** - browser should show confirmation
7. **Return to terminal** - setup should complete with success message

**Expected output:**
```
✓ Authentication successful
✓ Token saved to ~/.claude/skills/notebooklm/data/auth_info.json
```

---

#### Phase 3: Verification on HPC
After setup on local machine, test from HPC:
```bash
python scripts/run.py auth_manager.py status
```

**Expected output:**
```
✓ Authenticated as: your.email@gmail.com
✓ Token valid until: [expiry date]
```

---

### Troubleshooting Checklist

| Problem | Solution |
|---------|----------|
| Browser doesn't open in setup | Environment is headless. Run `setup` on local machine instead. |
| "ModuleNotFoundError" in run.py | Use local machine where Python environment is configured. |
| Rate limit error (50/day) | Each Google account has 50 query limit. Wait 24h or switch account. |
| Want to use different Google account | Run: `python scripts/run.py auth_manager.py clear` then re-run setup |
| Token expired | Run: `python scripts/run.py auth_manager.py reauth` (interactive on local machine) |

---

### Security Notes
- **Auth file location**: `~/.claude/skills/notebooklm/data/auth_info.json` (home directory, protected)
- **Never commit auth files** to repository
- **.gitignore protection**: Skill directory is git-ignored by design
- **Clear auth if sharing machine**: `python scripts/run.py auth_manager.py clear`
