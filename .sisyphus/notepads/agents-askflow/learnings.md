# Learnings - Agents Ask Flow
Timestamp: 2026-03-03

Updated AGENTS.md to mandate continuous `ask_questions` invocation. Agents must now call `ask_questions` after EVERY assistant message or reply, not just at task completion. This creates a persistent loop that only stops if the user explicitly selects a "Stop automatic asks" option (tracked via `persistent_stop: true`). Existing rules for headless/non-interactive operation (audit logs, auto-selection, `followup_required`) are preserved and applied to this continuous flow.

Lines added/modified in AGENTS.md:
- Line 11: Updated to "invoke `ask_questions` after every assistant message or reply finishes".
- Line 35: Added sentence about the loop continuing until "Stop automatic asks" is selected.
- Lines 93-125: Added "Continuous ask_questions policy" subsection with mandatory loop details, persistent stop control info, and an example JSON payload.
